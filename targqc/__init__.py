import os
from os.path import join, splitext, dirname

import Utils.reference_data as ref
from Utils.file_utils import safe_mkdir, can_reuse
from Utils.sambamba import index_bam
from Utils.parallel import parallel_view
from Utils.logger import info, critical, debug

from targqc import config
from targqc.fastq import proc_fastq
from targqc.region_coverage import make_region_reports
from targqc.general_report import make_general_reports
from targqc.Target import Target
from targqc.summarize import summarize_targqc

targqc_repr              = 'TargQC'
targqc_name              = 'targqc'

qualimap_name                   = 'qualimap'
qualimap_report_fname           = 'qualimapReport.html'
qualimap_genome_results_fname   = 'genome_results.txt'
qualimap_raw_data_dirname       = 'raw_data_qualimapReport'

qualimap_ishist_fname           = 'insert_size_histogram.txt'
qualimap_covhist_fname          = 'coverage_histogram.txt'
qualimap_gchist_fname           = 'mapped_reads_gc-content_distribution.txt'

picard_name              = 'picard'
picard_ishist_fname      = 'picard_ishist.txt'

fastqc_name              = 'fastqc'
dedup_bam                = 'dedup'  # work/post_processing/dedup and -ready.dedup.bam

fastqc_repr              = 'FastQC'
fastqc_report_fname      = 'fastqc_report.html'


def get_version():
    from targqc import version
    return version.__version__


def get_description():
    from targqc import version
    description = 'TargQC, target coverage evaluation tool. Version ' + version.__version__
    if version.__git_revision__:
        description += ', revision ' + version.__git_revision__
    return description


def start_targqc(work_dir, output_dir, samples, target_bed_fpath, parallel_cfg, bwa_prefix,
                 fai_fpath=None,
                 genome=config.genome,
                 depth_thresholds=config.depth_thresholds,
                 downsample_pairs_num=config.downsample_pairs_num,
                 padding=config.padding,
                 dedup=config.dedup,
                 reuse=config.reuse_intermediate,
                 is_debug=config.is_debug,
                 num_pairs_by_sample=None,
                 reannotate=config.reannotate,
                 ):
    d = get_description()
    info('*'*len(d))
    info(d)
    info('*'*len(d))
    info()

    fai_fpath = fai_fpath or ref.get_fai(genome)
    target = Target(work_dir, output_dir, fai_fpath, padding=padding, bed_fpath=target_bed_fpath,
         reannotate=reannotate, genome=genome, is_debug=is_debug)

    fastq_samples = [s for s in samples if not s.bam and s.l_fpath and s.r_fpath]
    if fastq_samples:
        if not bwa_prefix:
            critical('--bwa-prefix is required when running from fastq')
        with parallel_view(len(fastq_samples), parallel_cfg, join(work_dir, 'sge_fastq')) as view:
            num_pairs_by_sample = proc_fastq(fastq_samples, view, work_dir, bwa_prefix,
                 downsample_pairs_num, num_pairs_by_sample, dedup=dedup)

    info()
    for s in samples:
        if s.bam:
            info(s.name + ': using alignment ' + s.bam)

    with parallel_view(len(samples), parallel_cfg, join(work_dir, 'sge_bam')) as view:
        if all(can_reuse(s.bam + '.bai', s.bam) for s in samples):
            debug('BAM indexes exists')
        else:
            info('Indexing BAMs...')
            view.run(index_bam, [[s.bam] for s in samples])

        info('Making general reports...')
        make_general_reports(view, samples, target, genome, depth_thresholds, padding, num_pairs_by_sample,
                             is_debug=is_debug, reannotate=reannotate)

        info()
        info('Making region-level reports...')
        make_region_reports(view, work_dir, samples, target, genome, depth_thresholds)

    info()
    info('*' * 70)
    info('Summarizing TargQC coverage statistics for all samples...')
    summarize_targqc(parallel_cfg.threads, output_dir, work_dir, samples, bed_fpath=target_bed_fpath)

    # for general_report, per_gene_report, sample in zip(general_reports, per_gene_reports, samples):
    #     info('')
    #     info('*' * 70)
    #     if general_report.txt_fpath and verify_file(general_report.txt_fpath):
    #         info('Summary report: ' + general_report.txt_fpath)
    #     if per_gene_report:
    #         path = per_gene_report if isinstance(per_gene_report, basestring) else per_gene_report.txt_fpath
    #         if path and verify_file(path):
    #             info('Region-based report: ' + path)


class Sample:
    def __init__(self, name, dirpath, work_dir, bam=None, l_fpath=None, r_fpath=None,
                 genome=None, qualimap_dirpath=None, normal_match=None, read_pairs_num=None):
        self.name = name
        self.bam = bam
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.dedup_bam = None
        self.is_wgs = False
        self.qualimap_bed = None
        self.dirpath = dirpath
        self.work_dir = work_dir
        self.phenotype = None
        self.gender = None
        self.genome = None
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.avg_depth = None
        self.read_pairs_num = read_pairs_num

        self.targqc_dirpath                  = None
        self.targqc_html_fpath               = None
        self.targqc_json_fpath               = None
        self.targqc_region_txt               = None
        self.targqc_region_tsv               = None
        self.targqc_norm_depth_vcf_txt       = None
        self.targqc_norm_depth_vcf_tsv       = None

        if dirpath:
            self.targqc_dirpath = dirpath
            self.targqc_txt_fpath            = join(self.targqc_dirpath, 'summary.txt')
            self.targqc_html_fpath           = join(self.targqc_dirpath, 'summary.html')
            self.targqc_json_fpath           = join(self.targqc_dirpath, 'summary.json')
            self.targqc_region_txt           = join(self.targqc_dirpath, 'regions.txt')
            self.targqc_region_tsv           = join(self.targqc_dirpath, 'regions.tsv')
            self.targqc_norm_depth_vcf_txt   = None
            self.targqc_norm_depth_vcf_tsv   = None

        self.qualimap_dirpath                = None
        self.qualimap_html_fpath             = None
        self.qualimap_genome_results_fpath   = None
        self.qualimap_ins_size_hist_fpath    = None
        self.qualimap_cov_hist_fpath         = None
        self.qualimap_gc_hist_fpath          = None

        qualimap_dirpath = qualimap_dirpath or join(self.targqc_dirpath, qualimap_name)
        if qualimap_dirpath:
            self.qualimap_dirpath               = qualimap_dirpath
            self.qualimap_html_fpath            = join(self.qualimap_dirpath, qualimap_report_fname)
            self.qualimap_genome_results_fpath  = join(self.qualimap_dirpath, qualimap_report_fname)
            self.qualimap_raw_dirpath           = join(self.qualimap_dirpath, qualimap_raw_data_dirname)

            self.qualimap_ins_size_hist_fpath   = join(self.qualimap_raw_dirpath, qualimap_ishist_fname)
            self.qualimap_cov_hist_fpath        = join(self.qualimap_raw_dirpath, qualimap_covhist_fname)
            self.qualimap_gc_hist_fpath         = join(self.qualimap_raw_dirpath, qualimap_gchist_fname)

        self.picard_dirpath                 = join(self.targqc_dirpath, picard_name)
        self.picard_ins_size_hist_txt_fpath = join(self.picard_dirpath, picard_ishist_fname)
        self.picard_ins_size_hist_pdf_fpath = join(self.picard_dirpath, splitext(picard_ishist_fname)[0] + '.pdf')

    def __cmp__(self, other):
        return cmp(self.key_to_sort(), other.key_to_sort())

    def key_to_sort(self):
        parts = []

        cur_part = []
        prev_was_num = False

        for c in self.name:
            if prev_was_num == c.isdigit() and c not in ['-', '.']:  # same type of symbol, but not - or .
                cur_part.append(c)
            else:
                if cur_part:
                    part = ''.join(cur_part)
                    if prev_was_num:
                        part = int(part)
                    parts.append(part)
                    cur_part = []

                if c in ['-', '.']:
                    pass
                else:
                    if c.isdigit():
                        prev_was_num = True
                    else:
                        prev_was_num = False
                    cur_part.append(c)
        if cur_part:
            part = ''.join(cur_part)
            if prev_was_num:
                part = int(part)
            parts.append(part)

        return tuple(parts)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name