#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from collections import OrderedDict

import bcbio_postproc

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename, splitext
from optparse import OptionParser

import source
from source.config import Config, defaults
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_genome_resources, determine_run_cnf, \
    determine_sys_cnf
from source import logger
from source.logger import info, err, warn, critical, send_email
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, remove_quotes, \
    file_exists, isfile, splitext_plus
from source.main import set_up_dirs
from source.targetcov.bam_and_bed_utils import prepare_beds, extract_gene_names_and_filter_exons
from source.targetcov.submit_jobs import run_targqc
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed
from source.targetcov.summarize_targetcov import get_bed_targqc_inputs
from source.tools_from_cnf import get_system_path


def proc_args(argv):
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')
    # parser.add_option('--dedup', dest='dedup', action='store_true', default=False, help='count duplicates in coverage metrics')
    # parser.add_option('--no-dedup', dest='dedup', action='store_false', default=False, help='not counting duplicates in coverage metrics')
    parser.add_option('-e', '--extended', dest='extended', action='store_true', default=False, help='count missed variants')
    parser.add_option('--deep-seq', dest='deep_seq', action='store_true', default=False, help='deep targeted sequencing')
    parser.add_option('--no-qualimap', dest='qualimap', action='store_false', default=True, help='do not run qualimap')
    parser.add_option('--bed', dest='bed', help='BED file to run detailed coverage analysis.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('--downsample-to', dest='downsample_to', type='int')

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) == 0:
        critical('No input BAMs/FastQ provided. Usage:\n'
                 '1. ' + sys.argv[0] + ' *.fq.gz -o output_dir -g hg19 [--bed target.bed]\n'
                 '2. ' + sys.argv[0] + ' *.bam -o output_dir [--bed target.bed]\n')

    fastqs_by_sample, bam_by_sample = read_samples(args)

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'), is_targetseq=opts.deep_seq)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    cnf.output_dir = adjust_path(cnf.output_dir)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    cnf.proc_name = 'TargQC'
    set_up_dirs(cnf)
    # cnf.name = 'TargQC_' + cnf.project_name
    info(' '.join(sys.argv))

    check_genome_resources(cnf)

    target_bed, features_bed, genes_fpath = get_bed_targqc_inputs(cnf, cnf.bed)
    if not target_bed:
        info('No input BED is specified, using exons instead from ' + str(features_bed))

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        verify_file(cnf.qsub_runner, is_critical=True)

    return cnf, fastqs_by_sample, bam_by_sample, target_bed, features_bed, genes_fpath


def main():
    cnf, fastqs_by_sample, bam_by_sample, target_bed, exons_bed, genes_fpath = proc_args(sys.argv)

    samples = []
    for sname, (l, r) in fastqs_by_sample.items():
        s = source.TargQC_Sample(sname, join(cnf.output_dir, sname))
        s.l_fpath = l
        s.r_fpath = r
        samples.append(s)
    for sname, bam_fpath in bam_by_sample.items():
        s = source.TargQC_Sample(sname, join(cnf.output_dir, sname), bam=bam_fpath)
        samples.append(s)
    samples.sort(key=lambda _s: _s.key_to_sort())

    targqc_html_fpath = run_targqc(cnf, cnf.output_dir, samples, target_bed, exons_bed, genes_fpath)

    if targqc_html_fpath:
        send_email('TargQC report for ' + cnf.project_name + ':\n  ' + targqc_html_fpath)


def read_samples(args):
    bam_by_sample = find_bams(args)
    if bam_by_sample:
        info('Found ' + str(len(bam_by_sample)) + ' BAMs')

    input_not_bam = [verify_file(fpath) for fpath in args if adjust_path(fpath) not in bam_by_sample]
    input_not_bam = [fpath for fpath in input_not_bam if fpath]
    fastqs_by_sample = dict()
    if not input_not_bam and not bam_by_sample:
        critical('No correct input files')
    if input_not_bam:
        info(str(len(input_not_bam)) + ' correct input not-bam files')
        fastqs_by_sample = find_fastq_pairs(input_not_bam)
        if fastqs_by_sample:
            info('Found FastQ pairs: ' + str(len(fastqs_by_sample)))
        intersection = set(fastqs_by_sample.keys()) & set(bam_by_sample.keys())
        if intersection:
            critical('The following samples both had input BAMs and FastQ: ' + ', '.join(list(intersection)))

    return fastqs_by_sample, bam_by_sample


def find_bams(args):
    bam_by_sample = OrderedDict()
    bad_bam_fpaths = []

    good_args = []
    for arg in args:
        # /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio,Kudos159 /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0160_BHKWMNADXX/bcbio,Kudos160
        fpath = arg.split(',')[0]
        fname, ext = splitext(fpath)
        if ext == '.bam':
            bam_fpath = verify_bam(fpath)
            if not bam_fpath:
                bad_bam_fpaths.append(bam_fpath)
            else:
                if len(arg.split(',')) > 1:
                    sname = arg.split(',')[1]
                else:
                    sname = basename(splitext(bam_fpath)[0])
                bam_by_sample[sname] = bam_fpath
                good_args.append(arg)
    if bad_bam_fpaths:
        critical('BAM files cannot be found, empty or not BAMs:' + ', '.join(bad_bam_fpaths))
    for arg in good_args:
        args.remove(arg)

    return bam_by_sample


def find_fastq_pairs(fpaths):
    info('Finding fastq pairs.')
    fastqs_by_sample_name = dict()
    for fpath in fpaths:
        fn, ext = splitext_plus(basename(fpath))
        if ext in ['.fq', '.fq.gz', '.fastq', '.fastq.gz']:
            sname, l_fpath, r_fpath = None, None, None
            if fn.endswith('_1'):
                sname = fn[:-2]
                l_fpath = fpath
            if fn.endswith('_R1'):
                sname = fn[:-3]
                l_fpath = fpath
            if fn.endswith('_2'):
                sname = fn[:-2]
                r_fpath = fpath
            if fn.endswith('_R2'):
                sname = fn[:-3]
                r_fpath = fpath
            if not sname:
                sname = fn
                info('Cannot detect file for ' + sname)

            l, r = fastqs_by_sample_name.get(sname, (None, None))
            if l and l_fpath:
                critical('Duplicated left FastQ files for ' + sname + ': ' + l + ' and ' + l_fpath)
            if r and r_fpath:
                critical('Duplicated right FastQ files for ' + sname + ': ' + r + ' and ' + r_fpath)
            fastqs_by_sample_name[sname] = l or l_fpath, r or r_fpath

    fixed_fastqs_by_sample_name = dict()
    for sname, (l, r) in fastqs_by_sample_name.items():
        if not l:
            err('ERROR: for sample ' + sname + ', left reads not found')
        if not r:
            err('ERROR: for sample ' + sname + ', left reads not found')
        fixed_fastqs_by_sample_name[sname] = l, r

    return fixed_fastqs_by_sample_name


if __name__ == '__main__':
    main()