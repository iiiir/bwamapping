#!/bin/env python
import sys
import os
import argparse
import subprocess
import sjm
import util
import random
import string

p = argparse.ArgumentParser(description='run_bwa_aln_pe_ubam.py -b *.bam -o o.bam -O folder --tmp /scratch_space -j jobs.sjm')
p.add_argument('-b', '--ubams', metavar='FILE', nargs="+", required=True, help='The FASTQ file(s) for reads 1')
p.add_argument('-o', '--output_prefix', metavar='NAME', required=True, help='The output prefix for final bam')
p.add_argument('-O', '--outdir', metavar='DIR', required=True, help='The output directory')
p.add_argument('--merge_bams', action='store_true', help='If files need to be merged')
p.add_argument('-T', '--tmp', metavar='DIR', required=True, help='The TMP directory for storing intermediate files (default=output directory')
p.add_argument('-j', '--jobfile', metavar='FILE', help='The jobfile name (default: stdout)')
p.add_argument('-m', '--memory', metavar='SIZE', type=int, default=12, help='Memory size (GB) per job (default: 12)')
p.add_argument('-p', '--project', metavar='NAME', default="CompBio", help='Project name for jobs (default: CompBio)')
p.add_argument('-q', '--queue', metavar='NAME', default="normal", help='Queue for jobs (default: normal)')
p.add_argument('-t', '--threads', metavar='COUNT', type=int, default=4, help='Number of threads for BWA alignment, only works for SGE (default: 4)')
p.add_argument('--account', metavar='STR', default="swang", help='Accounting string for the purpose of cluster accounting.')
p.add_argument('--submit', action='store_true', help='Submit the jobs')
args = p.parse_args()

outdir = util.Dir(args.outdir)
outdir.mkdirs()

tmplett = ''.join(random.sample(string.ascii_letters,6))
tmpdir  = util.Dir(os.path.join(args.tmp, tmplett))
tmpdir.mkdirs()

if args.jobfile is None:
    jobfile=None
else:
    jobfile=util.File(args.jobfile)

sjm.Job.name_prefix="BWA-mapping"+"."
sjm.Job.memory="%sG"%args.memory
sjm.Job.queue="pcgp"
sjm.Job.project="CompBio"

tmpdir = getattr(__builtins__, 'str')(tmpdir)
outdir = getattr(__builtins__, 'str')(outdir)

def sort_ubam(ubams):
    jobs=[]
    for ubam in ubams:
        ubam = util.File(ubam)
        obam = util.File( os.path.join(tmpdir, os.path.basename(ubam.path.rstrip('u.bam') + '.bam') ) )
        job = sjm.Job('picard_sortUbam-%s' % ubam.prefix)
        job.memory = "20G"
        job.input  = ubam
        job.output = obam
        job.append('picard_sortUbam.sh %s %s' % (job.input, job.output) )
        jobs.append(job)
    return jobs

def align_pe(pjobs):
    jobs=[]
    for pjob in pjobs:
        inbam = pjob.output
        obam  = inbam.chext('aln.bam')
        job   = sjm.Job('bwa_aln_pe-%s' % inbam.prefix)
        job.memory = "20G"
        job.input  = inbam
        job.output = obam
        job.append('bwa_aln_pe_qn.sh %s %s %s'%(job.output, inbam, inbam) )
        job.depend(pjob)
        jobs.append(job)
    return jobs

def merge_aln(pjobs):
    jobs = []
    for pjob in pjobs:
        alnbam = pjob.output
        ubam   = pjob.input
        job = sjm.Job('picard_mergeBam-%s' % alnbam.name )
        job.memory = "10G"
        job.output = util.File( alnbam.path.rstrip('.aln.bam') + '.sort.bam' ) 
        job.append('picard_mergeBam.sh %s %s %s' % (job.output, alnbam, ubam) )
        job.depend(pjob)
        jobs.append(job)
    return jobs

def dedup_merge(pjobs, outbam):
    jobs = []
    bams = []
    for pjob in pjobs:
        bams.append(pjob.output.path)
    job = sjm.Job('samtools_merge-%s' % outbam )
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

jobs = sort_ubam(args.ubams)
jobs = align_pe(jobs)
jobs = merge_aln(jobs)
jobs = dedup_merge(jobs, args.output_prefix)
jobs = sam_flagstat(jobs)
descout = sys.stdout if jobfile is None else open(jobfile.path, "w")
descout.write(sjm.Job().depend(*jobs).desc())
descout.flush()

if args.submit:
    print >> sys.stderr, "Submitting jobs (%s) through SJM"%jobfile
    os.system("sjm %s &" %jobfile)
