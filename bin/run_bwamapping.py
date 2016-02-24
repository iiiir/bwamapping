#!/bin/env python
import sys
import os
import re
import argparse
import glob
import subprocess
import sjm
import util

p = argparse.ArgumentParser(description='Generating the job file for the BWA mapping')
p.add_argument('-r1', '--read1', metavar='FILE', nargs="+", required=True, help='The FASTQ file(s) for reads 1')
p.add_argument('-r2', '--read2', metavar='FILE', nargs="+", help='The FASTQ file(s) for reads 2, if paired-end')
p.add_argument('-o', '--output', metavar='DIR', required=True, help='The output directory')
p.add_argument('-T', '--tmp', metavar='DIR', help='The TMP directory for storing intermediate files (default=output directory')
p.add_argument('--id',metavar='STR', required=True, help="sample id")
p.add_argument('--lb',metavar='STR', required=True, help="library")
p.add_argument('--sm',metavar='STR', required=True, help="sample name")
p.add_argument('--pl',metavar='STR', required=True, help="sequencing platform")
p.add_argument('-j', '--jobfile', metavar='FILE', help='The jobfile name (default: stdout)')
p.add_argument('-m', '--memory', metavar='SIZE', type=int, default=12, help='Memory size (GB) per job (default: 12)')
p.add_argument('-p', '--project', metavar='NAME', default="CompBio", help='Project name for jobs (default: CompBio)')
p.add_argument('-q', '--queue', metavar='NAME', default="normal", help='Queue for jobs (default: normal)')
p.add_argument('-t', '--threads', metavar='COUNT', type=int, default=4, help='Number of threads for BWA alignment, only works for SGE (default: 4)')
p.add_argument('--split', metavar='SIZE', help='Split fastq files before aignment. Value in millions (e.g. 10 implies 10 million reads).')
p.add_argument('--account', metavar='STR', default="swang", help='Accounting string for the purpose of cluster accounting.')
p.add_argument('--submit', action='store_true', help='Submit the jobs')
args = p.parse_args()

outdir=util.Dir(args.output)
outdir.mkdirs()

tmpdir=outdir
if (args.tmp is not None):
    tmpdir=util.Dir(args.tmp)
tmpdir.mkdirs()

sample=args.sample

readgroup = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:%s' % (args.id, args.lb, args.sm, args.pl)

sjm.Job.name_prefix="BWA-mapping"+"."
sjm.Job.memory="%sG"%args.memory
sjm.Job.queue="normal"
sjm.Job.project="CompBio"

tmpdir = getattr(__builtins__, 'str')(tmpdir)

def splitfq(fqname, out_prefix, split_lines):
    ''' 
    Do not support gz file
    This step must finish first to produce file list for mapping
    '''
    fq   = util.File(fqname)
    job  = sjm.Job('Split-fastq-%s' % fq.prefix)
    job.output = os.path.join(tmpdir, out_prefix)
    job.append('split.py %s %s %d' % (fq.path, out_prefix, split_lines) )
    return job

def align(pjob1,pjob2):
    jobs=[]
    reads1 = glob.glob('%s*' % pjob1.output)
    reads2 = glob.glob('%s*' % pjob2.output)
    for i in range(0, len(reads1)):
        job = sjm.Job('bwa_mem_pe-%s' % readfile1.prefix)
        job.output = re.sub(r'[._][Rr]1', '', read1[i].rstrip('.fastq') )+ '.sorted'
        job.append('bwa_mem_pe.sh %s %s %s %s'%(job.output, reads1, reads2, readgroup))
        job.depend(pjob1).denpend(pjob2)
        jobs.append(job)
    return jobs

def dedup(pjobs):
    jobs=[]
    for pjob in pjobs:
        bam=pjob.output
    job=Job(jname)
    job.memory = "12G"
    job.sge_options="-A %s"%args.account
    job.append('%s %s %s'%(cmd, input, output))
    job.append('sam_index.sh %s' % output)
    job.output=output
    return job

job1 = splitfq(args.read1)
job2 = splitfq(args.read2)
jobs = align(job1, job2)

descout = sys.stdout if jobfile is None else open(jobfile.path, "w")
descout.write(Job().depend(*jobs).desc())
descout.flush()

if args.submit:
    print >> sys.stderr, "Submitting jobs (%s) through SJM"%jobfile
    os.system("sjm %s &" %jobfile)
