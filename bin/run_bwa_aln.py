#!/bin/env python
import sys
import os
import re
import argparse
import glob
import subprocess
import sjm
import util

p = argparse.ArgumentParser(description='run_bwa_aln.py -r1 R1.fq -r2 R2.fq -o o.bam -O work_dir --tmp /scratch_space -j sample.sjm --id Test --sm Test --lb Test --pl ILN')
p.add_argument('-r1', '--read1', metavar='FILE', nargs="+", required=True, help='The FASTQ file(s) for reads 1')
p.add_argument('-r2', '--read2', metavar='FILE', nargs="+", required=True, help='The FASTQ file(s) for reads 2, if ubam then provide read1')
p.add_argument('-o', '--output_prefix', metavar='NAME', required=True, help='The output prefix for final bam')
p.add_argument('-O', '--outdir', metavar='DIR', required=True, help='The output directory')
p.add_argument('--merge_bams', action='store_true', help='If files need to be merged')
p.add_argument('-T', '--tmp', metavar='DIR', required=True, help='The TMP directory for storing intermediate files (default=output directory')
p.add_argument('--id',metavar='STR', required=True, help="sample id")
p.add_argument('--lb',metavar='STR', required=True, help="library")
p.add_argument('--sm',metavar='STR', required=True, help="sample name (will creat such a folder in TMP dir to hold temp files)")
p.add_argument('--pl',metavar='STR', required=True, help="sequencing platform")
p.add_argument('-j', '--jobfile', metavar='FILE', help='The jobfile name (default: stdout)')
p.add_argument('-m', '--memory', metavar='SIZE', type=int, default=12, help='Memory size (GB) per job (default: 12)')
p.add_argument('-p', '--project', metavar='NAME', default="CompBio", help='Project name for jobs (default: CompBio)')
p.add_argument('-q', '--queue', metavar='NAME', default="normal", help='Queue for jobs (default: normal)')
p.add_argument('-t', '--threads', metavar='COUNT', type=int, default=4, help='Number of threads for BWA alignment, only works for SGE (default: 4)')
p.add_argument('--account', metavar='STR', default="swang", help='Accounting string for the purpose of cluster accounting.')
p.add_argument('--submit', action='store_true', help='Submit the jobs')
args = p.parse_args()

outdir=util.Dir(args.outdir)
outdir.mkdirs()

tmpdir=util.Dir(os.path.join(args.tmp, args.sm))
tmpdir.mkdirs()

if args.jobfile is None:
    jobfile=None
else:
    jobfile=util.File(args.jobfile)

readgroup = "'@RG\\\\tID:%s\\\\tLB:%s\\\\tSM:%s\\\\tPL:%s'" % (args.id, args.lb, args.sm, args.pl)

sjm.Job.name_prefix="BWA-mapping"+"."
sjm.Job.memory="%sG"%args.memory
sjm.Job.queue="pcgp"
sjm.Job.project="CompBio"

tmpdir = getattr(__builtins__, 'str')(tmpdir)
outdir = getattr(__builtins__, 'str')(outdir)

def align_pe(reads1, reads2):
    jobs=[]
    for i in range(0, len(reads1)):
        read1 = reads1[i]
        read2 = reads2[i]
        readfile1 = util.File(read1)
        readfile2 = util.File(read2)
        if readfile1.path.endswith('.gz'):
            bamname = re.sub(r'[._][Rr]1', '', readfile1.prefix.rstrip('.fastq') ) + '.sorted.bam'
        else:
            bamname   = re.sub(r'[._][Rr]1', '', readfile1.prefix ) + '.sorted.bam'
        bam = util.File( os.path.join(tmpdir, bamname) )
        job = sjm.Job('bwa_aln_pe-%s' % readfile1.prefix)
        job.output = bam
        job.memory = "15G"
        job.append('bwa_aln_pe.sh %s %s %s %s'%(job.output, read1, read2, readgroup))
        jobs.append(job)
    return jobs

def dedup_merge(pjobs, outbam):
    jobs = []
    bams = []
    for pjob in pjobs:
        bams.append(pjob.output.path)
    job = sjm.Job('picard_mdup-%s' % outbam )
    job.memory = "20G"
    job.output = util.File( os.path.join(outdir, outbam) )
    job.append('picard_mdup.sh %s %s'%(job.output, ' '.join(bams) ) )
    job.depend(*pjobs)
    jobs.append(job)
    return jobs

def sam_flagstat(pjobs):
    jobs = []
    for pjob in pjobs:
        bam=util.File(pjob.output)
        job=sjm.Job('samtools-flagstat-%s' % bam.prefix)
        job.memory = "10G"
        job.output = bam.chext("flagstat.txt")
        job.append('samtools flagstat %s > %s'%(bam, job.output))
        job.depend(pjob)
        jobs.append(job)
    return jobs

jobs = align_pe(args.read1, args.read2)
jobs = dedup_merge(jobs, args.output_prefix)
jobs = sam_flagstat(jobs)
descout = sys.stdout if jobfile is None else open(jobfile.path, "w")
descout.write(sjm.Job().depend(*jobs).desc())
descout.flush()

if args.submit:
    print >> sys.stderr, "Submitting jobs (%s) through SJM"%jobfile
    os.system("sjm %s &" %jobfile)
