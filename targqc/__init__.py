from collections import OrderedDict
from os.path import join, splitext

from Utils.file_utils import safe_mkdir
from Utils.sambamba import index_bam
from Utils.parallel import parallel_view
from Utils.logger import info, critical

from targqc.fastq import proc_fastq
from targqc.region_coverage import make_region_reports
from targqc.general_report import make_general_reports
import targqc.config as cfg

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


def start_targqc(work_dir, samples, target,
                 parallel_cfg=None,
                 bwa_prefix=None,
                 genome=None,
                 depth_thresholds=None,
                 downsample_pairs_num=None,
                 padding=None,
                 dedup=None,
                 reuse=None,
                 is_debug=None,
                 ):
    cfg.parallel_cfg = parallel_cfg if parallel_cfg is not None else cfg.parallel_cfg
    cfg.genome = genome if genome is not None else cfg.genome
    cfg.depth_thresholds = depth_thresholds if depth_thresholds is not None else cfg.depth_thresholds
    cfg.downsample_pairs_num = downsample_pairs_num if downsample_pairs_num is not None else cfg.downsample_pairs_num
    cfg.padding = padding if padding is not None else cfg.padding
    cfg.dedup = dedup if dedup is not None else cfg.dedup
    cfg.reuse_intermediate = reuse if reuse is not None else cfg.reuse_intermediate
    cfg.is_debug = is_debug if is_debug is not None else cfg.is_debug

    fastq_samples = [s for s in samples if not s.bam and s.l_fpath and s.r_fpath]
    num_reads_by_sample = dict()
    if fastq_samples:
        if not bwa_prefix:
            critical('--bwa-prefix is required when running from fastq')
        with parallel_view(len(fastq_samples), cfg.parallel_cfg, join(work_dir, 'sge_fastq')) as view:
            num_reads_by_sample = proc_fastq(fastq_samples, view, work_dir, bwa_prefix,
                                             downsample_pairs_num, dedup)

    for s in samples:
        if s.bam:
            info(s.name + ': using alignment ' + s.bam)

    info()
    with parallel_view(len(samples), cfg.parallel_cfg, join(work_dir, 'sge_bam')) as view:
        info('Indexing BAMs...')
        view.run(index_bam, safe_mkdir(join(work_dir, 'index_bam')), [[s.bam] for s in samples])

        info('Making general reports...')
        make_general_reports(view, samples, target, num_reads_by_sample)

        info()
        info('Making region-level reports...')
        make_region_reports(view, work_dir, samples, target)

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
                 genome=None, qualimap_dirpath=None, normal_match=None):
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

        self.targqc_dirpath                  = None
        self.targqc_html_fpath               = None
        self.targqc_json_fpath               = None
        self.targqc_region_txt               = None
        self.targqc_region_tsv               = None
        self.targqc_norm_depth_vcf_txt       = None
        self.targqc_norm_depth_vcf_tsv       = None

        self.qualimap_dirpath                = None
        self.qualimap_html_fpath             = None
        self.qualimap_genome_results_fpath   = None
        self.qualimap_ins_size_hist_fpath    = None
        self.qualimap_cov_hist_fpath         = None
        self.qualimap_gc_hist_fpath          = None

        self.fastqc_dirpath                  = None
        self.fastqc_html_fpath               = None

        self.picard_dirpath                  = None
        self.picard_ins_size_hist_txt_fpath  = None
        self.picard_ins_size_hist_pdf_fpath  = None

        if dirpath:
            self.targqc_dirpath = dirpath
            self.targqc_txt_fpath            = join(self.targqc_dirpath, 'summary.txt')
            self.targqc_html_fpath           = join(self.targqc_dirpath, 'summary.html')
            self.targqc_json_fpath           = join(self.targqc_dirpath, 'summary.json')
            self.targqc_region_txt           = join(self.targqc_dirpath, 'regions.txt')
            self.targqc_region_tsv           = join(self.targqc_dirpath, 'regions.tsv')
            self.targqc_norm_depth_vcf_txt   = None
            self.targqc_norm_depth_vcf_tsv   = None

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