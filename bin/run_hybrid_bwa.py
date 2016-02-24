#!/bin/env python

# run_bwamapping.py -r1 `ls_ubam.sh /rgs01/resgen/prod/tartan/index/data/Baker/Baker/SJMMNORM017056_C3/CHIPSEQ_INPUT/bam/SJMMNORM017056_C3.bam` \
# -r2 `ls_ubam.sh /rgs01/resgen/prod/tartan/index/data/Baker/Baker/SJMMNORM017056_C3/CHIPSEQ_INPUT/bam/SJMMNORM017056_C3.bam` \
# -o SJMMNORM017056_C3 \
# -O SJMMNORM017056_C3 \
# --id SJMMNORM017056_C3 \
# --sm SJMMNORM017056_C3 \
# --lb CHIPSEQ_INPUT \
# --pl illumina \
# --tmp /rgs01/scratch_space/swang/ \
# -j SJMMNORM017056_C3.sjm

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
p.add_argument('-r2', '--read2', metavar='FILE', nargs="+", required=True, help='The FASTQ file(s) for reads 2, if ubam then provide read1')
p.add_argument('-o', '--output_prefix', metavar='NAME', required=True, help='The output prefix for final bam')
p.add_argument('-O', '--outdir', metavar='DIR', required=True, help='The output directory')
p.add_argument('-T', '--tmp', metavar='DIR', required=True, help='The TMP directory for storing intermediate files (default=output directory')
p.add_argument('--id',metavar='STR', required=True, help="sample id")
p.add_argument('--lb',metavar='STR', required=True, help="library")
p.add_argument('--sm',metavar='STR', required=True, help="sample name")
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

mm9_chrs = '1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 MT X Y NT_166325 NT_166464 NT_166452 NT_166480 NT_166448 NT_166458 NT_166443 NT_166466 NT_166476 NT_166479 NT_166478 NT_166474 NT_166471 NT_166445 NT_166465 NT_166457 NT_166470 NT_166454 NT_166472 NT_166449 NT_166481 NT_166337 NT_166459 NT_166456 NT_166473 NT_166461 NT_166475 NT_166462 NT_166444 NT_166453 NT_166446 NT_166469 NT_072868 NT_166335 NT_166467 NT_166283 NT_166338 NT_166340 NT_166442 NT_166334 NT_166286 NT_166451 NT_166336 NT_166339 NT_166290 NT_053651 NT_166450 NT_166447 NT_166468 NT_166460 NT_166477 NT_166455 NT_166291 NT_166463 NT_166433 NT_166402 NT_166327 NT_166308 NT_166309 NT_109319 NT_166282 NT_166314 NT_166303 NT_112000 NT_110857 NT_166280 NT_166375 NT_166311 NT_166307 NT_166310 NT_166323 NT_166437 NT_166374 NT_166364 NT_166439 NT_166328 NT_166438 NT_166389 NT_162750 NT_166436 NT_166372 NT_166440 NT_166326 NT_166342 NT_166333 NT_166435 NT_166434 NT_166341 NT_166376 NT_166387 NT_166281 NT_166313 NT_166380 NT_166360 NT_166441 NT_166359 NT_166386 NT_166356 NT_166357 NT_166423 NT_166384 NT_161879 NT_161928 NT_166388 NT_161919 NT_166381 NT_166367 NT_166392 NT_166406 NT_166365 NT_166379 NT_166358 NT_161913 NT_166378 NT_166382 NT_161926 NT_166345 NT_166385 NT_165789 NT_166368 NT_166405 NT_166390 NT_166373 NT_166361 NT_166348 NT_166369 NT_161898 NT_166417 NT_166410 NT_166383 NT_166362 NT_165754 NT_166366 NT_166363 NT_161868 NT_166407 NT_165793 NT_166352 NT_161925 NT_166412 NT_165792 NT_161924 NT_166422 NT_165795 NT_166354 NT_166350 NT_165796 NT_161904 NT_166370 NT_165798 NT_165791 NT_161885 NT_166424 NT_166346 NT_165794 NT_166377 NT_166418 NT_161877 NT_166351 NT_166408 NT_166349 NT_161906 NT_166391 NT_161892 NT_166415 NT_165790 NT_166420 NT_166353 NT_166344 NT_166371 NT_161895 NT_166404 NT_166413 NT_166419 NT_161916 NT_166347 NT_161875 NT_161911 NT_161897 NT_161866 NT_166409 NT_161872 NT_166403 NT_161902 NT_166414 NT_166416 NT_166421 NT_161923 NT_161937'

dm3_chrs = '2L_dm3 2LHet_dm3 2R_dm3 2RHet_dm3 3L_dm3 3LHet_dm3 3R_dm3 3RHet_dm3 4_dm3 U_dm3 Uextra_dm3 X_dm3 XHet_dm3 YHet_dm3 dmel_mitochondrion_genome'

sjm.Job.name_prefix="BWA-mapping"+"."
sjm.Job.memory="%sG"%args.memory
sjm.Job.queue="pcgp"
sjm.Job.project="CompBio"

tmpdir = getattr(__builtins__, 'str')(tmpdir)
outdir = getattr(__builtins__, 'str')(outdir)

def align_se(reads1, reads2):
    jobs=[]
    for i in range(0, len(reads1)):
        read1 = reads1[i]
        read2 = reads2[i]
        readfile1 = util.File(read1)
        readfile2 = util.File(read2)
        bamname   = re.sub(r'[._][Rr]1', '', readfile1.prefix ) + '.sorted.bam'
        bam = util.File( os.path.join(tmpdir, bamname) )
        job = sjm.Job('bwa_aln_se-%s' % readfile1.prefix)
        job.output = bam
        job.append('bwa_aln_se.sh %s %s %s %s'%(job.output, read1, read2, readgroup))
        jobs.append(job)
    return jobs

def split_species(pjobs):
    jobs_mm9 = []
    jobs_dm3 = []
    for pjob in pjobs:
        bamfile = pjob.output
        job1 = sjm.Job('samtools_mm9-%s' % bamfile.prefix)
        job1.memory = "12G"
        job1.output = bamfile.chext('mm9.bam')
        job1.append('samtools view -hb %s %s > %s'%(bamfile, mm9_chrs, job1.output))
        job1.depend(*pjobs)
        jobs_mm9.append(job1)
        job2 = sjm.Job('samtools_dm3-%s' % bamfile.prefix)
        job2.memory = "12G"
        job2.output = bamfile.chext('dm3.bam')
        job2.append('samtools view -hb %s %s > %s'%(bamfile, dm3_chrs, job2.output))
        job2.depend(*pjobs)
        jobs_dm3.append(job2)
    return jobs_mm9, jobs_dm3

def merge_bam(pjobs, out_prefix, suffix=None):
    '''
    Caveat: If output bam exists, needs to apply "-f" to overwrite or task will abort.
    '''
    bams = []
    for pjob in pjobs:
        bams.append(pjob.output.path)
    job = sjm.Job('samtools_merge-%s' % suffix)
    job.memory = "5G"
    outname    = os.path.join(tmpdir, '%s.%s.bam' % (out_prefix, suffix))
    job.output = util.File( outname )
    job.append('samtools merge %s %s && samtools index %s'%(job.output, ' '.join(bams), job.output))
    job.depend(*pjobs)
    return job

def dedup_bam(pjobs):
    jobs = []
    for pjob in pjobs:
        bamfile=pjob.output
        job=sjm.Job('picard_mdup-%s' % bamfile.prefix)
        job.memory = "12G"
        job.output = os.path.join(outdir, bamfile.chext("mdup.bam").name)
        job.append('picard_mdup.sh %s %s'%(job.output, bamfile))
        job.depend(pjob)
        jobs.append(job)
    return jobs

def sam_flagstat(pjobs):
    jobs = []
    for pjob in pjobs:
        bam=util.File(pjob.output)
        job=sjm.Job('samtools-flagstat-%s' % bam.prefix)
        job.memory = "5G"
        job.output = bam.chext("flagstat.txt")
        job.append('samtools flagstat %s > %s'%(bam, job.output))
        job.depend(pjob)
        jobs.append(job)
    return jobs

jobs = align_se(args.read1, args.read2)
jobs_mm9, jobs_dm3 = split_species(jobs)
job_all = merge_bam(jobs, args.output_prefix, 'all')
job_mm9 = merge_bam(jobs_mm9, args.output_prefix, 'mm9')
job_dm3 = merge_bam(jobs_dm3, args.output_prefix, 'dm3')
jobs = [job_all, job_mm9, job_dm3]
jobs = dedup_bam(jobs)
jobs = sam_flagstat(jobs)
descout = sys.stdout if jobfile is None else open(jobfile.path, "w")
descout.write(sjm.Job().depend(*jobs).desc())
descout.flush()

if args.submit:
    print >> sys.stderr, "Submitting jobs (%s) through SJM"%jobfile
    os.system("sjm %s &" %jobfile)
