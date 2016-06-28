# coding=utf-8

import math
from collections import OrderedDict
from os.path import join, abspath, realpath, dirname, relpath

import targqc
import targqc.config as tc
from targqc.qualimap.report_parser import parse_qualimap_sample_report
from targqc.qualimap.runner import run_qualimap

from Utils.Region import calc_bases_within_threshs, calc_rate_within_normal, Region, GeneInfo
from Utils.bam_bed_utils import index_bam, get_padded_bed_file, number_mapped_reads_on_target, calc_region_number, \
    intersect_bed, calc_sum_of_regions, number_of_mapped_reads, count_bed_cols, \
    sambamba_depth
from Utils.call_process import run
from Utils.file_utils import intermediate_fname, verify_file, safe_mkdir
from Utils.logger import critical, info, err, warn, debug
from Utils.reporting.reporting import ReportSection, Metric, MetricStorage, SampleReport, PerRegionSampleReport, \
    write_tsv_rows, write_txt_rows


def get_header_metric_storage(depth_thresholds, is_wgs=False, padding=None):
    sections = [
        ReportSection('reads', 'Reads', [
            Metric('Reads', short_name='Reads'),
            Metric('Mapped reads', short_name='Mapped', description='samtools view -c -F 4', ok_threshold='Percentage of mapped reads', bottom=0, multiqc=dict(hidden=True)),
            Metric('Percentage of mapped reads', short_name='Mapped', unit='%', ok_threshold=0.98, bottom=0),
            Metric('Properly paired mapped reads percent', short_name='Mapped paired', unit='%', description='Pecent of properly paired mapped reads.', ok_threshold=0.9, bottom=0),
            Metric('Properly paired reads percent', short_name='Paired', unit='%', description='Pecent of properly paired reads (-f 2).', ok_threshold=0.9, bottom=0),
            Metric('Duplication rate', short_name='Dup rate', description='Percent of mapped reads (-F 4), marked as duplicates (-f 1024)', quality='Less is better', unit='%'),
            Metric('Read min length', short_name='Min len', description='Read minimum length', is_hidden=True, multiqc=dict(hidden=True)),
            Metric('Read max length', short_name='Max len', description='Read maximum length'),
            Metric('Read mean length', short_name='Avg read len', description='Read average length'),
            Metric('Sex', short_name='sex', is_hidden=True),
            Metric('Median insert size', short_name='IS', description='Median insert size'),
            Metric('Median GC', short_name='GC', description='Median GC-content', unit='%'),
        ]),
    ]
    if not is_wgs:
        sections.extend([
            ReportSection('target_metrics', 'Target coverage', [
                Metric('Covered bases in target', short_name='Trg covered', multiqc=dict(hidden=True, name='Target covered bases'), description='Target bases covered by at least 1 read', unit='bp'),
                Metric('Percentage of target covered by at least 1 read', short_name='%', multiqc=dict(name='Target covered'), unit='%'),
                Metric('Percentage of reads mapped on target', short_name='reads on trg', unit='%', description='Percentage of unique mapped reads overlapping target at least by 1 base'),
                Metric('Percentage of reads mapped off target', short_name='reads off trg', unit='%', quality='Less is better', description='Percentage of unique mapped reads that don\'t overlap target even by 1 base'),
                Metric('Percentage of reads mapped on padded target', short_name='reads on trg' + (' %dbp pad' % padding if padding is not None else ' w/ padding'), unit='%', description='Percentage of reads that overlap target at least by 1 base. Should be 1-2% higher.'),
                Metric('Percentage of usable reads', short_name='Usable reads', unit='%', description='Share of unique reads mapped on target in the total number of original reads (reported in the very first column Reads'),
                Metric('Read bases mapped on target', short_name='Read bp on trg', unit='bp', multiqc=dict(hidden=True)),
            ]),
        ])
    else:
        sections.extend([
            ReportSection('target_metrics_wgs', 'Genome coverage', [
                Metric('Covered bases in genome', short_name='Genome covered', multiqc=dict(name='Genome bases covered', hidden=True), description='Genome bases covered by at least 1 read', unit='bp'),
                Metric('Percentage of genome covered by at least 1 read', short_name='%', multiqc=dict(name='Target covered', hidden=True), unit='%'),
                Metric('Covered bases in exome', short_name='CDS covered', multiqc=dict(name='CDS bases covered', hidden=True), description='Covered CDS bases. CDS coordinates are taken from RefSeq', unit='bp'),
                Metric('Percentage of exome covered by at least 1 read', short_name='%', multiqc=dict(name='CDS covered'), description='Percentage of CDS covered by at least 1 read. CDS coordinates are taken from RefSeq', unit='%'),
                Metric('Percentage of reads mapped on exome', short_name='Reads on CDS', multiqc=dict(name='On CDS'), unit='%', description='Percentage of reads mapped on CDS. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of reads mapped off exome', short_name='Off CDS', unit='%', quality='Less is better', description='Percentage of reads mapped outside of CDS. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of usable reads', short_name='Usable reads', unit='%', description='Share of mapped unique reads in all reads (reported in the very first column Reads)'),
            ]),
        ])

    trg_name = 'target' if not is_wgs else 'genome'
    depth_section = ReportSection('depth_metrics', ('Target' if not is_wgs else 'Genome') + ' coverage depth', [
        Metric('Median ' + trg_name + ' coverage depth', short_name='Median', multiqc=dict(name='Depth')),
        Metric('Average ' + trg_name + ' coverage depth', short_name='Avg', multiqc=dict(name='Avg')),
        Metric('Std. dev. of ' + trg_name + ' coverage depth', short_name='Std dev', quality='Less is better'),
        # Metric('Minimal ' + trg_name + ' coverage depth', short_name='Min', is_hidden=True),
        # Metric('Maximum ' + trg_name + ' coverage depth', short_name='Max', is_hidden=True),
        Metric('Percentage of ' + trg_name + ' within 20% of mean depth', short_name='&#177;20% avg', unit='%')
    ])
    for depth in depth_thresholds:
        name = 'Part of ' + trg_name + ' covered at least by ' + str(depth) + 'x'
        depth_section.add_metric(Metric(name, short_name=str(depth) + 'x', description=name, unit='%', multiqc=dict(hidden=True)))
    sections.append(depth_section)

    sections.append(
        ReportSection('qualimap', 'Qualimap stats' + ('' if is_wgs else ' within the target'), [
            Metric('Mean Mapping Quality',  'Mean MQ',            'Mean mapping quality, inside of regions'),
            Metric('Mismatches',            'Mismatches',         'Mismatches, inside of regions', quality='Less is better', multiqc=dict(hidden=True)),
            Metric('Insertions',            'Insertions',         'Insertions, inside of regions', quality='Less is better', multiqc=dict(hidden=True)),
            Metric('Deletions',             'Deletions',          'Deletions, inside of regions', quality='Less is better', multiqc=dict(hidden=True)),
            Metric('Homopolymer indels',    'Homopolymer indels', 'Percentage of homopolymer indels, inside of regions', quality='Less is better', multiqc=dict(hidden=True)),
            Metric('Qualimap',              'Qualimap report',    'Qualimap report', multiqc=dict(hidden=True)),
        ])
    )

    ms = MetricStorage(
        general_section=ReportSection('general_section', '', [
            Metric('Target', short_name='Target', common=True),
            Metric('Reference size', short_name='Reference bp', unit='bp', common=True),
            # Metric('Target ready', short_name='Target ready', common=True),
            Metric('Regions in target', short_name='Regions in target', common=True),
            Metric('Bases in target', short_name='Target bp', unit='bp', common=True),
            Metric('Percentage of reference', short_name='Percentage of reference', unit='%', common=True),
            # Metric('Genes', short_name='Genes', common=True),
            Metric('Genes in target', short_name='Genes in target', common=True),
            Metric('Median human GC', short_name='Median human GC content', unit='%', common=True),
        ]),
        sections=sections
    )
    return ms


class TargetInfo:
    def __init__(self, fpath=None, bed=None, original_target_bed=None, regions_num=None,
                 bases_num=None, fraction=None, genes_fpath=None, genes_num=None):
        self.fpath = realpath(fpath) if fpath else None  # raw source file - to demonstrate where we took one
        self.bed = bed                                   # processed (sorted, merged...), to do real calculations
        self.original_target_bed = original_target_bed
        self.regions_num = regions_num
        self.bases_num = bases_num
        self.fraction = fraction
        self.genes_fpath = realpath(genes_fpath) if genes_fpath else None
        self.genes_num = genes_num
        self.type = None  # 'Regional', 'WGS'


# def _dedup_and_flag_stat(cnf, bam_fpath):
#     bam_stats = flag_stat(cnf, bam_fpath)
#     info('Total reads: ' + Metric.format_value(bam_stats['total']))
#     info('Total mapped reads: ' + Metric.format_value(bam_stats['mapped']))
#     info('Total dup reads: ' + Metric.format_value(bam_stats['duplicates']))
#     info('Total properly paired reads: ' + Metric.format_value(bam_stats['properly paired']))
#
#     # dedup_bam_dirpath = join(cnf.work_dir, source.dedup_bam)
#     # safe_mkdir(dedup_bam_dirpath)
#     # dedup_bam_fpath = join(dedup_bam_dirpath, add_suffix(basename(bam_fpath), source.dedup_bam))
#     # remove_dups(cnf, bam_fpath, output_fpath=dedup_bam_fpath, use_grid=False)
#     # dedup_bam_stats = samtools_flag_stat(cnf, dedup_bam_fpath)
#     # info('Total reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['total']))
#     # info('Total mapped reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['mapped']))
#     return bam_stats  # dedup_bam_fpath, bam_stats, dedup_bam_stats


def parse_qualimap_insert_size(qualimap_insert_size_fpath):
    d = dict()
    zero_insertsize = 0
    with open(qualimap_insert_size_fpath) as f:
        for l in f:
            if not l.strip() or l.startswith('#'):
                continue
            insertsize, count = l.split(None, 1)
            insertsize = int(round(float(insertsize)))
            count = float(count) / 1000000
            if insertsize == 0:
                zero_insertsize = count
            else:
                d[insertsize] = count

    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    median_insert_size = None
    for thisins, thiscount in d.items():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_insert_size = thisins
            break
    return median_insert_size


def parse_qualimap_gc_content(qualimap_gc_content_fpath):
    avg_gc = 0
    avg_human_gc = 0
    with open(qualimap_gc_content_fpath) as f:
        for l in f:
            if l.startswith('#'):
                continue
            sections = l.split(None, 3)
            gc = int(round(float(sections[0])))         # %
            content = float(sections[1]) / 100.0        # share of reads with this GC %
            human_content = float(sections[2]) / 100.0  # share of human reference with this GC %
            avg_gc += gc * content
            avg_human_gc += gc * human_content
    return avg_gc, avg_human_gc


def parse_qualimap_coverage_hist(qualimap_coverage_hist_fpath):
    bases_by_depth = OrderedDict()
    with open(qualimap_coverage_hist_fpath) as f:
        for l in f:
            if l.startswith('#'):
                pass
            else:
                cov, bases = map(int, map(float, l.strip().split()))
                bases_by_depth[cov] = bases

    # calculating median coverage
    num_counts = sum(bases_by_depth.values())
    cum_counts = 0
    median_coverage = None
    for thiscov, thiscount in bases_by_depth.items():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_coverage = thiscov
            break

    return bases_by_depth, median_coverage


def parse_qualimap_results(sample, depth_thresholds):
    if not verify_file(sample.qualimap_html_fpath):
        critical('Qualimap report was not found')

    qualimap_value_by_metric = parse_qualimap_sample_report(sample.qualimap_html_fpath)
    bases_by_depth, median_depth = parse_qualimap_coverage_hist(sample.qualimap_cov_hist_fpath)
    median_gc, median_human_gc = parse_qualimap_gc_content(sample.qualimap_gc_hist_fpath)
    median_ins_size = parse_qualimap_insert_size(sample.qualimap_ins_size_hist_fpath)

    def find_rec(name, percent=False, on_target=True):
        if on_target:
            name_on_target = name + ' (on target)'
            if percent:
                name_on_target += ' %'
            res = qualimap_value_by_metric.get(name_on_target)
            if res:
                return res
        if percent:
            name += ' %'
        return qualimap_value_by_metric.get(name)

    depth_stats = dict(
        ave_depth       = find_rec('Coverage Mean'),
        stddev_depth    = find_rec('Coverage Standard Deviation'),
        median_depth    = median_depth,
        bases_by_depth  = bases_by_depth
    )
    target_stats = dict(
        reference_size  = find_rec('Reference size'),
        target_size     = find_rec('Regions size/percentage of reference'),
        target_fraction = find_rec('Regions size/percentage of reference', percent=True),
    )
    reads_stats = dict(
        total                    = find_rec('Number of reads'),
        mapped                   = find_rec('Mapped reads', on_target=False),
        mapped_rate              = find_rec('Mapped reads', percent=True, on_target=False),
        unmapped                 = find_rec('Unmapped reads'),
        unmapped_rate            = find_rec('Unmapped reads', percent=True),
        mapped_on_target         = find_rec('Mapped reads (on target)'),
        mapped_rate_on_target    = find_rec('Mapped reads (on target)', percent=True),
        mapped_paired            = find_rec('Mapped paired reads', on_target=False),
        mapped_paired_rate       = find_rec('Mapped paired reads', percent=True, on_target=False),
        paired                   = find_rec('Paired reads'),
        paired_rate              = find_rec('Paired reads', percent=True),
        dup                      = find_rec('Duplicated reads (flagged)'),
        dup_rate                 = find_rec('Duplicated reads (flagged)', percent=True),
        min_len                  = find_rec('Read min length'),
        max_len                  = find_rec('Read max length'),
        ave_len                  = find_rec('Read mean length'),
        median_gc                = median_gc,
        median_human_gc          = median_human_gc,
        median_ins_size          = median_ins_size,
    )
    indels_stats = dict(
        mean_mq     = find_rec('Mean Mapping Quality'),
        mismatches  = find_rec('Mismatches'),
        insertions  = find_rec('Insertions'),
        deletions   = find_rec('Deletions'),
        homo_indels = find_rec('Homopolymer indels'),
    )
    return depth_stats, reads_stats, indels_stats, target_stats


chry_key_regions_by_genome = {
    'hg19': join(dirname(abspath(__file__)), 'gender', 'chrY.hg19.bed'),
    'hg38': join(dirname(abspath(__file__)), 'gender', 'chrY.hg38.bed'),
}
MALE_TARGET_REGIONS_FACTOR = 0.7
AVE_DEPTH_THRESHOLD_TO_DETERMINE_SEX = 5
FEMALE_Y_COVERAGE_FACTOR = 10.0

def _determine_sex(work_dir, bam_fpath, ave_depth, target_bed=None):
    info()
    info('Determining sex')

    male_genes_bed = None
    for k in chry_key_regions_by_genome:
        if k in tc.genome:
            male_genes_bed = chry_key_regions_by_genome.get(k)
            break
    if not male_genes_bed:
        warn('Warning: no male key regions for ' + tc.genome + ', cannot identify sex')
        return None

    male_area_size = calc_sum_of_regions(male_genes_bed)
    info('Male region total size: ' + str(male_area_size))

    if target_bed:
        male_genes_bed = intersect_bed(work_dir, target_bed, male_genes_bed)
        target_male_area_size = calc_sum_of_regions(male_genes_bed)
        if target_male_area_size < male_area_size * MALE_TARGET_REGIONS_FACTOR:
            info('Target male region total size is ' + str(target_male_area_size) + ', which is less than the ' +
                 'checked male regions size * ' + str(MALE_TARGET_REGIONS_FACTOR) +
                 ' (' + str(male_area_size * MALE_TARGET_REGIONS_FACTOR) + ') - cannot determine sex')
            return None
        else:
            info('Target male region total size is ' + str(target_male_area_size) + ', which is higher than the ' +
                 'checked male regions size * ' + str(MALE_TARGET_REGIONS_FACTOR) +
                 ' (' + str(male_area_size * MALE_TARGET_REGIONS_FACTOR) + '). ' +
                 'Determining sex based on coverage in those regions.')
    else:
        info('WGS, determining sex based on chrY key regions coverage.')

    info('Detecting sex by comparing the Y chromosome key regions coverage and average coverage depth.')
    if not bam_fpath:
        critical('BAM file is required.')
    index_bam(bam_fpath)

    chry_cov_output_fpath = sambamba_depth(work_dir, male_genes_bed, bam_fpath, [], reuse=tc.reuse_intermediate)
    chry_mean_coverage = get_mean_cov(chry_cov_output_fpath)
    info('Y key regions average depth: ' + str(chry_mean_coverage))
    ave_depth = float(ave_depth)
    info('Sample average depth: ' + str(ave_depth))
    if ave_depth < AVE_DEPTH_THRESHOLD_TO_DETERMINE_SEX:
        info('Sample average depth is too low (less then ' + str(AVE_DEPTH_THRESHOLD_TO_DETERMINE_SEX) +
             ') - cannot determine sex')
        return None

    if chry_mean_coverage == 0:
        info('Y depth is 0 - it\s female')
        sex = 'F'
    else:
        factor = ave_depth / chry_mean_coverage
        info('Sample depth / Y depth = ' + str(factor))
        if factor > FEMALE_Y_COVERAGE_FACTOR:  # if mean target coverage much higher than chrY coverage
            info('Sample depth is more than ' + str(FEMALE_Y_COVERAGE_FACTOR) + ' times higher than Y depth - it\s female')
            sex = 'F'
        else:
            info('Sample depth is not more than ' + str(FEMALE_Y_COVERAGE_FACTOR) + ' times higher than Y depth - it\s male')
            sex = 'M'
    info('Sex is ' + sex)
    info()
    return sex


def get_mean_cov(bedcov_output_fpath):
    mean_cov = []
    mean_cov_col = None
    total_len = 0
    with open(bedcov_output_fpath) as bedcov_file:
        for line in bedcov_file:
            if line.startswith('#'):
                mean_cov_col = line.split('\t').index('meanCoverage')
                continue
            line_tokens = line.replace('\n', '').split()
            start, end = map(int, line_tokens[1:3])
            size = end - start
            mean_cov.append(float(line_tokens[mean_cov_col]) * size)
            total_len += size
    mean_cov = sum(mean_cov) / total_len if total_len > 0 else 0
    return mean_cov


def make_general_report(work_dir, sample, target_bed, gene_by_name_and_chrom):
    target_info = TargetInfo(
        fpath=target_bed, bed=target_bed, original_target_bed=tc.original_target_bed or target_bed,
        genes_num=len(gene_by_name_and_chrom) if gene_by_name_and_chrom else None)
    if target_bed:
        target_info.regions_num = calc_region_number(target_bed)

    run_qualimap(sample.qualimap_dirpath, sample.bam, target_bed, threads=tc.threads_one_sample)

    # try:
    #     picard_ins_size_hist(sample, sample.bam)
    # except subprocess.CalledProcessError as e:
    #     err('Picard insert size histogram command exit with error code ' + str(e.returncode) + ':\n ' + str(e.cmd))
    # else:
    #     os.remove(sample.qualimap_ins_size_hist_fpath)

    depth_stats, reads_stats, indels_stats, target_stats = parse_qualimap_results(sample, tc.depth_thresholds)

    reads_stats['gender'] = _determine_sex(work_dir, sample.bam, depth_stats['ave_depth'], target_bed)
    info()

    if 'bases_by_depth' in depth_stats:
        depth_stats['bases_within_threshs'], depth_stats['rates_within_threshs'] = calc_bases_within_threshs(
            depth_stats['bases_by_depth'],
            target_stats['target_size'] if target_bed else target_stats['reference_size'],
            tc.depth_thresholds)

        depth_stats['wn_20_percent'] = calc_rate_within_normal(
            depth_stats['bases_by_depth'],
            depth_stats['ave_depth'],
            target_stats['target_size'] if target_bed else target_stats['reference_size'])

    if target_stats['target_size']:
        target_info.bases_num = target_stats['target_size']
        target_info.fraction  = target_stats['target_fraction']
    else:
        target_info.bases_num = target_stats['reference_size']

    # if sample.dedup_bam:
    #     bam_fpath = sample.dedup_bam
    reads_stats['mapped_dedup'] = number_of_mapped_reads(work_dir, sample.bam, dedup=True, reuse=tc.reuse_intermediate)

    if target_info.bed:
        reads_stats['mapped_dedup_on_target'] = number_mapped_reads_on_target(
            work_dir, target_bed, sample.bam, dedup=True, reuse=tc.reuse_intermediate) or 0

    if target_info.bed:
        padded_bed = get_padded_bed_file(work_dir, target_info.bed, tc.padding, tc.fai_fpath)
        reads_stats['mapped_dedup_on_padded_target'] = number_mapped_reads_on_target(
            work_dir, padded_bed, sample.bam, dedup=True, reuse=tc.reuse_intermediate) or 0
    elif tc.cds_bed_fpath:
        info('Using the CDS reference BED ' + tc.cds_bed_fpath + ' to calc "reads on CDS"')
        reads_stats['mapped_dedup_on_exome'] = number_mapped_reads_on_target(
            work_dir, tc.cds_bed_fpath, sample.bam, dedup=True, reuse=tc.reuse_intermediate) or 0
    # elif features_no_genes_bed:
    #     info('Using ensemble ' + features_no_genes_bed + ' to calc reads on exome')
    #     reads_stats['mapped_dedup_on_exome'] = number_mapped_reads_on_target(cnf, features_no_genes_bed, bam_fpath) or 0

    summary_report = _build_report(work_dir, depth_stats, reads_stats, indels_stats, sample, target_info)

    return summary_report, depth_stats['ave_depth']


def _build_report(cnf, depth_stats, reads_stats, mm_indels_stats, sample, target_info):
    report = SampleReport(sample, metric_storage=get_header_metric_storage(tc.depth_thresholds, is_wgs=target_info.bed is None, padding=tc.padding))
    report.add_record('Qualimap', value='Qualimap', url=relpath(sample.qualimap_html_fpath, sample.dirpath), silent=True)
    if reads_stats.get('gender') is not None:
        report.add_record('Sex', reads_stats['gender'], silent=True)

    info('* General coverage statistics *')
    report.add_record('Reads', reads_stats['total'])
    report.add_record('Mapped reads', reads_stats['mapped'])
    # report.add_record('Unmapped reads', reads_stats['totaAvgl'] - reads_stats['mapped'])
    percent_mapped = 1.0 * (reads_stats['mapped'] or 0) / reads_stats['total'] if reads_stats['total'] else None
    assert percent_mapped <= 1.0 or percent_mapped is None, str(percent_mapped)
    report.add_record('Percentage of mapped reads', percent_mapped)
    # percent_unmapped = 1.0 * (reads_stats['total'] - reads_stats['mapped']) / reads_stats['total'] if reads_stats['total'] else None
    # assert percent_unmapped <= 1.0 or percent_unmapped is None, str(percent_unmapped)
    # report.add_record('Percentage of unmapped reads', percent_unmapped)
    if reads_stats.get('mapped_paired') is not None:
        total_paired_reads_pecent = 1.0 * (reads_stats['mapped_paired'] or 0) / reads_stats['total'] if reads_stats['total'] else None
        assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
        report.add_record('Properly paired mapped reads percent', total_paired_reads_pecent)
    if reads_stats.get('paired') is not None:
        total_paired_reads_pecent = 1.0 * (reads_stats['paired'] or 0) / reads_stats['total'] if reads_stats['total'] else None
        assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
        report.add_record('Properly paired reads percent', total_paired_reads_pecent)
    # if dedup_bam_stats:
    # dup_rate = 1 - (1.0 * dedup_bam_stats['mapped'] / bam_stats['mapped']) if bam_stats['mapped'] else None
    report.add_record('Duplication rate', reads_stats['dup_rate'])
    # report.add_record('Dedupped mapped reads', reads_stats['mapped'] - reads_stats[''])
    report.add_record('Median GC', reads_stats['median_gc'])
    report.add_record('Median human GC', reads_stats['median_human_gc'])
    report.add_record('Median insert size', reads_stats['median_ins_size'])

    info('')

    if target_info.bed:
        info('* Target coverage statistics *')
        if target_info.original_target_bed:
            report.add_record('Target', target_info.original_target_bed)
            if count_bed_cols(target_info.original_target_bed) == 3:
                report.add_record('Ready target (sorted and annotated)', target_info.fpath)
        else:
            report.add_record('Target', target_info.fpath)
        report.add_record('Bases in target', target_info.bases_num)
        report.add_record('Percentage of reference', target_info.fraction)
        report.add_record('Regions in target', target_info.regions_num)
    else:
        info('* Genome coverage statistics *')
        report.add_record('Target', 'whole genome')
        report.add_record('Reference size', target_info.bases_num)

    report.add_record('Genes in target', target_info.genes_num)

    trg_type = 'target' if target_info.bed else 'genome'

    if 'bases_within_threshs' in depth_stats:
        bases_within_threshs = depth_stats['bases_within_threshs']
        v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
        v_percent_covered_bases_in_targ = 1.0 * (v_covered_bases_in_targ or 0) / target_info.bases_num if target_info.bases_num else None
        assert v_percent_covered_bases_in_targ <= 1.0 or v_percent_covered_bases_in_targ is None, str(v_percent_covered_bases_in_targ)

        report.add_record('Covered bases in ' + trg_type, v_covered_bases_in_targ)
        report.add_record('Percentage of ' + trg_type + ' covered by at least 1 read', v_percent_covered_bases_in_targ)

    if target_info.bed:
        info('Getting number of mapped reads on target...')
        # mapped_reads_on_target = number_mapped_reads_on_target(cnf, target_info.bed, bam_fpath)
        if 'mapped_dedup_on_target' in reads_stats:
            # report.add_record('Reads mapped on target', reads_stats['mapped_on_target'])
            info('Unique mapped on target: ' + str(reads_stats['mapped_dedup_on_target']))
            percent_mapped_dedup_on_target = 1.0 * reads_stats['mapped_dedup_on_target'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
            report.add_record('Percentage of reads mapped on target', percent_mapped_dedup_on_target)
            assert percent_mapped_dedup_on_target <= 1.0 or percent_mapped_dedup_on_target is None, str(percent_mapped_dedup_on_target)

            percent_mapped_dedup_off_target = 1.0 * (reads_stats['mapped_dedup'] - reads_stats['mapped_dedup_on_target']) / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
            report.add_record('Percentage of reads mapped off target', percent_mapped_dedup_off_target)
            assert percent_mapped_dedup_off_target <= 1.0 or percent_mapped_dedup_off_target is None, str(percent_mapped_dedup_off_target)

            percent_usable = 1.0 * reads_stats['mapped_dedup_on_target'] / reads_stats['total'] if reads_stats['total'] != 0 else None
            report.add_record('Percentage of usable reads', percent_usable)
            assert percent_usable <= 1.0 or percent_usable is None, str(percent_usable)

        read_bases_on_targ = int(target_info.bases_num * depth_stats['ave_depth'])  # sum of all coverages
        report.add_record('Read bases mapped on target', read_bases_on_targ)

        if 'mapped_dedup_on_padded_target' in reads_stats:
            # report.add_record('Reads mapped on padded target', reads_stats['mapped_reads_on_padded_target'])
            percent_mapped_on_padded_target = 1.0 * reads_stats['mapped_dedup_on_padded_target'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] else None
            report.add_record('Percentage of reads mapped on padded target', percent_mapped_on_padded_target)
            assert percent_mapped_on_padded_target <= 1.0 or percent_mapped_on_padded_target is None, str(percent_mapped_on_padded_target)

    elif 'mapped_dedup_on_exome' in reads_stats:
        # report.add_record('Reads mapped on target', reads_stats['mapped_on_target'])
        percent_mapped_on_exome = 1.0 * reads_stats['mapped_dedup_on_exome'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
        if percent_mapped_on_exome:
            report.add_record('Percentage of reads mapped on exome', percent_mapped_on_exome)
            assert percent_mapped_on_exome <= 1.0 or percent_mapped_on_exome is None, str(percent_mapped_on_exome)
            percent_mapped_off_exome = 1.0 - percent_mapped_on_exome
            report.add_record('Percentage of reads mapped off exome ', percent_mapped_off_exome)

        percent_usable = 1.0 * reads_stats['mapped_dedup'] / reads_stats['total'] if reads_stats['total'] != 0 else None
        report.add_record('Percentage of usable reads', percent_usable)
        assert percent_usable <= 1.0 or percent_usable is None, str(percent_usable)

    info('')
    report.add_record('Average ' + trg_type + ' coverage depth', depth_stats['ave_depth'])
    report.add_record('Median ' + trg_type + ' coverage depth', depth_stats['median_depth'])
    report.add_record('Std. dev. of ' + trg_type + ' coverage depth', depth_stats['stddev_depth'])
    # report.add_record('Minimal ' + trg_type + ' coverage depth', depth_stats['min_depth'])
    # report.add_record('Maximum ' + trg_type + ' coverage depth', depth_stats['max_depth'])
    if 'wn_20_percent' in depth_stats:
        report.add_record('Percentage of ' + trg_type + ' within 20% of mean depth', depth_stats['wn_20_percent'])
        assert depth_stats['wn_20_percent'] <= 1.0 or depth_stats['wn_20_percent'] is None, str( depth_stats['wn_20_percent'])

    if 'bases_within_threshs' in depth_stats:
        for depth, bases in depth_stats['bases_within_threshs'].items():
            fraction_val = 1.0 * (bases or 0) / target_info.bases_num if target_info.bases_num else 0
            if fraction_val > 0:
                report.add_record('Part of ' + trg_type + ' covered at least by ' + str(depth) + 'x', fraction_val)
            assert fraction_val <= 1.0 or fraction_val is None, str(fraction_val)
    info()

    report.add_record('Read mean length', reads_stats['ave_len'])
    report.add_record('Read min length', reads_stats['min_len'])
    report.add_record('Read max length', reads_stats['max_len'])
    report.add_record('Mean Mapping Quality', mm_indels_stats['mean_mq'])
    report.add_record('Mismatches', mm_indels_stats['mismatches'])
    report.add_record('Insertions', mm_indels_stats['insertions'])
    report.add_record('Deletions', mm_indels_stats['deletions'])
    report.add_record('Homopolymer indels', mm_indels_stats['homo_indels'])

    report.save_json(sample.targqc_json_fpath)
    report.save_txt(sample.targqc_txt_fpath)
    report.save_html(sample.targqc_html_fpath, caption='Target coverage statistics for ' + sample.name, debug=tc.debug)
    info()
    info('Saved to ')
    info('  ' + report.txt_fpath)
    return report
