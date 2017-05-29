# coding=utf-8
import os
from collections import OrderedDict
from os.path import join, abspath, realpath, dirname, relpath
from pybedtools import BedTool

from ensembl import get_merged_cds
from ngs_utils import reference_data, logger
from ngs_utils.bed_utils import get_padded_bed_file, intersect_bed, calc_sum_of_regions, count_bed_cols,\
    calc_bases_within_threshs, calc_rate_within_normal
from ngs_utils.sambamba import index_bam, number_mapped_reads_on_target, number_of_mapped_reads, sambamba_depth
from ngs_utils.file_utils import intermediate_fname, verify_file, safe_mkdir, can_reuse
from ngs_utils.logger import critical, info, err, warn, debug
from ngs_utils.reporting.reporting import ReportSection, Metric, MetricStorage, SampleReport

from targqc.qualimap import report_parser, runner


def get_header_metric_storage(depth_threshs, is_wgs=False, padding=None):
    sections = [
        ReportSection('reads', 'Reads', [
            Metric('Original reads',                       short_name='Orig reads',   multiqc=dict(order=1, kind='reads', min=0)),
            Metric('Reads',                                short_name='Reads',        multiqc=dict(order=1, kind='reads', min=0)),
            Metric('Mapped reads',                         short_name='Mapped',       multiqc=dict(title='Mapped reads', hidden=True, order=2, kind='reads', min=0),  ok_threshold='Percentage of mapped reads', bottom=0, description='samtools view -c -F 4'),
            Metric('Percentage of mapped reads',           short_name='%',            multiqc=dict(title='Mapped', order=3, kind='reads'),                  unit='%', ok_threshold=0.98, bottom=0),
            Metric('Properly paired mapped reads percent', short_name='Paired',       multiqc=dict(order=4, kind='reads'),                                  unit='%', ok_threshold=0.9,  bottom=0, description='Pecent of properly paired mapped reads.'),
            Metric('Duplication rate',                     short_name='Dups',         multiqc=dict(order=5, kind='reads'),                                  unit='%', quality='Less is better',    description='Percent of mapped reads (-F 4), marked as duplicates (-f 1024)'),
            Metric('Read min length',                      short_name='Min len',      multiqc=dict(title='Read min len', hidden=True, kind='reads', min=0), unit='bp',                             description='Read minimum length', is_hidden=True),
            Metric('Read mean length',                     short_name='Avg read len', multiqc=dict(title='Read avg len', order=22, kind='reads', min=0),    unit='bp',                             description='Read average length'),
            Metric('Read max length',                      short_name='Max len',      multiqc=dict(title='Read max len', order=23, kind='reads', min=0),    unit='bp',                             description='Read maximum length'),
            Metric('Median insert size',                   short_name='IS',           multiqc=dict(title='Med. IS', order=21, kind='reads', min=0),         unit='bp',                             description='Median insert size'),
            Metric('Median GC',                            short_name='GC',           multiqc=dict(title='Med. GC', order=21, kind='other'),                unit='%',                              description='Med. GC-content', ),
        ]),
    ]
    if not is_wgs:
        sections.extend([
            ReportSection('target_metrics', 'Target coverage', [
                Metric('Covered bases in target',                         short_name='Trg covered',                         multiqc=dict(hidden=True, title='Trg cov bases', kind='cov'),   unit='bp',                          description='Target bases covered by at least 1 read'),
                Metric('Percentage of target covered by at least 1 read', short_name='%',                                   multiqc=dict(title='Trg cov', order=14, kind='cov'),            unit='%'),
                Metric('Percentage of reads mapped on target',            short_name='Reads on trg',                        multiqc=dict(title='On trg', order=10, kind='trg'),             unit='%',                           description='Percentage of unique mapped reads overlapping target at least by 1 base'),
                Metric('Percentage of reads mapped off target',           short_name='Reads off trg',                       multiqc=dict(hidden=True, title='Off trg', order=11, kind='trg'),            unit='%', quality='Less is better', description='Percentage of unique mapped reads that don\'t overlap target even by 1 base'),
                Metric('Percentage of reads mapped on padded target',     short_name='On trg &#177;' + str(padding) + 'bp', multiqc=dict(order=12, kind='trg'), unit='%',                                                       description='Percentage of reads that overlap target at least by 1 base. Should be 1-2% higher.'),
                Metric('Percentage of usable reads',                      short_name='Usable reads',                        multiqc=dict(order=13, title='Usable', kind='trg'),             unit='%',                           description='Share of unique reads mapped on target in the total number of original reads (reported in the very first column Reads'),
                Metric('Read bases mapped on target',                     short_name='Read bp on trg',                      multiqc=dict(hidden=True, kind='cov'),                          unit='bp'),
            ]),
        ])
    else:
        sections.extend([
            ReportSection('target_metrics_wgs', 'Genome coverage', [
                Metric('Covered bases in genome',                         short_name='Genome covered', multiqc=dict(title='Genome cov bases', hidden=True, kind='cov'),     unit='bp',                          description='Genome bases covered by at least 1 read'),
                Metric('Percentage of genome covered by at least 1 read', short_name='%',              multiqc=dict(title='Genome cov', hidden=True, order=13, kind='cov'), unit='%'),
                Metric('Covered bases in exome',                          short_name='CDS covered',    multiqc=dict(title='CDS cov bases', hidden=True, kind='cov'),        unit='bp',                          description='Covered CDS bases. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of exome covered by at least 1 read',  short_name='%',              multiqc=dict(title='CDS cov', order=14, kind='cov'),                 unit='%',                           description='Percentage of CDS covered by at least 1 read. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of reads mapped on exome',             short_name='Reads on CDS',   multiqc=dict(title='On CDS', order=10, kind='trg'),                  unit='%',                           description='Percentage of reads mapped on CDS. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of reads mapped off exome',            short_name='Off CDS',        multiqc=dict(order=11, kind='trg'),                                  unit='%', quality='Less is better', description='Percentage of reads mapped outside of CDS. CDS coordinates are taken from RefSeq'),
                Metric('Percentage of usable reads',                      short_name='Usable reads',   multiqc=dict(order=12, title='Usable', kind='trg'),                  unit='%',                           description='Share of mapped unique reads in all reads (reported in the very first column Reads)'),
            ]),
        ])

    trg_name = 'target' if not is_wgs else 'genome'
    depth_section = ReportSection('section_name', ('Target' if not is_wgs else 'Genome') + ' coverage depth', [
        Metric('Median ' + trg_name + ' coverage depth',                  short_name='Median',         multiqc=dict(title='Med depth', order=6, kind='cov', min=0)),
        Metric('Mean ' + trg_name + ' coverage depth',                    short_name='Mean',           multiqc=dict(title='Mean depth', order=7, kind='cov', min=0)),
        Metric('Std. dev. of ' + trg_name + ' coverage depth',            short_name='Std dev',        multiqc=dict(order=8, kind='cov', min=0),                                quality='Less is better'),
        Metric('Percentage of ' + trg_name + ' within 20% of med depth',  short_name='&#177;20% med',  multiqc=dict(order=9, kind='cov'),                             unit='%')
    ])
    for depth in depth_threshs:
        name = 'Part of ' + trg_name + ' covered at least by ' + str(depth) + 'x'
        depth_section.add_metric(Metric(name,                             short_name=str(depth) + 'x', multiqc=dict(hidden=True, kind='cov'),                         unit='%', description=name))
    depth_section.add_metric(
        Metric('Estimated ' + trg_name + ' full coverage depth',          short_name='Est full mean',  multiqc=dict(title='Est mean depth', order=7, kind='cov', min=0), description='Estimated mean coverage of full dataset. Calculated as (the total number of raw reads * downsampled mapped reads fraction / total downsampled mapped reads) * downsampled average coverage'),
    )
    sections.append(depth_section)

    sections.append(
        ReportSection('other', 'Other stats' + ('' if is_wgs else ' within the target'), [
            Metric('Mean Mapping Quality',  short_name='Mean MQ',            multiqc=dict(kind='other', min=0),                         description='Mean mapping quality, inside of regions'),
            Metric('Mismatches',            short_name='Mismatches',         multiqc=dict(hidden=True, kind='other', min=0),            description='Mismatches, inside of regions',                       quality='Less is better'),
            Metric('Insertions',            short_name='Insertions',         multiqc=dict(hidden=True, kind='other', min=0),            description='Insertions, inside of regions',                       quality='Less is better'),
            Metric('Deletions',             short_name='Deletions',          multiqc=dict(hidden=True, kind='other', min=0),            description='Deletions, inside of regions',                        quality='Less is better'),
            Metric('Homopolymer indels',    short_name='Homopolymer indels', multiqc=dict(hidden=True, kind='other', min=0),            description='Percentage of homopolymer indels, inside of regions', quality='Less is better'),
            Metric('Qualimap',              short_name='Qualimap report',    multiqc=dict(hidden=True, kind='other'),                   description='Qualimap report'),
            Metric('Sex',                   short_name='sex',                multiqc=dict(order=20,    kind='other'), is_hidden=True),
        ])
    )

    ms = MetricStorage(
        general_section=ReportSection('general_section', '', [
            Metric('Target',                                     short_name='Target',                            common=True),
            Metric('Ready target (clean, sorted and annotated)', short_name='Ready target',                      common=True),
            Metric('Reference size',                             short_name='Reference bp', unit='bp',           common=True),
            Metric('Regions in target',                          short_name='Regions in target',                 common=True),
            Metric('Bases in target',                            short_name='Target bp', unit='bp',              common=True),
            Metric('Percentage of reference',                    short_name='Percentage of reference', unit='%', common=True),
            Metric('Genes in target',                            short_name='Genes in target',                   common=True),
            Metric('Scope',                                      short_name='Scope',                             common=True, is_hidden=True),
        ]),
        sections=sections
    )
    return ms


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
            avg_gc += gc * content
            if len(sections) > 2:
                human_content = float(sections[2]) / 100.0  # share of human reference with this GC %
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


def parse_qualimap_results(sample):
    if not verify_file(sample.qualimap_html_fpath):
        critical('QualiMap report was not found')

    qualimap_value_by_metric = report_parser.parse_qualimap_sample_report(sample.qualimap_html_fpath)
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

def determine_sex(work_dir, bam_fpath, ave_depth, genome, target_bed=None):
    info()
    info('Determining sex')

    male_bed = None
    for k in chry_key_regions_by_genome:
        if k in genome:
            male_bed = BedTool(chry_key_regions_by_genome.get(k))
            break
    if not male_bed:
        warn('Warning: no male key regions for ' + genome + ', cannot identify sex')
        return None

    male_area_size = male_bed.count()
    info('Male region total size: ' + str(male_area_size))

    if target_bed:
        male_bed = BedTool(target_bed).intersect(male_bed).merge()
        target_male_area_size = male_bed.count()
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

    chry_cov_output_fpath = sambamba_depth(work_dir, male_bed, bam_fpath, [])
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


def _qualimap_outputs(sample):
    return [v for k, v in sample.__dict__.items() if k.startswith('qualimap_') and k.endswith('_fpath')]


def make_general_reports(view, samples, target, genome, depth_threshs, bed_padding,
                         num_pairs_by_sample=None, reuse=False, is_debug=False, reannotate=False, fai_fpath=None):
    if all(all(can_reuse(fp, [s.bam, target.qualimap_bed_fpath] if target.bed else s.bam)
               for fp in _qualimap_outputs(s))
           for s in samples):
        debug('All QualiMap files for all samples exist and newer than BAMs and BEDs, reusing')
    else:
        info('Running QualiMap...')
        view.run(runner.run_qualimap,
            [[s.work_dir, s.qualimap_dirpath, _qualimap_outputs(s), s.bam, genome, target.qualimap_bed_fpath, view.cores_per_job]
             for s in samples])

        for s in samples:
            for fp in _qualimap_outputs(s):
                verify_file(fp, is_critical=True)

    summary_reports = []

    for sample in samples:
        info('-'*70)
        info(sample.name)
        debug('-'*70)
        debug('Parsing QualiMap results...')
        depth_stats, reads_stats, indels_stats, target_stats = parse_qualimap_results(sample)

        _prep_report_data(sample, depth_stats, reads_stats, indels_stats, target_stats,
                          target, num_pairs_by_sample, genome, depth_threshs, fai_fpath=fai_fpath)

        r = _build_report(depth_stats, reads_stats, indels_stats, sample, target,
                          depth_threshs, bed_padding, sample_num=len(samples), is_debug=is_debug,
                          reannotate=reannotate)
        summary_reports.append(r)

    return summary_reports


def _prep_report_data(sample, depth_stats, reads_stats, indels_stats, target_stats,
                      target, num_pairs_by_sample, genome, depth_threshs, fai_fpath=None):
    sample.avg_depth = depth_stats['ave_depth']

    if num_pairs_by_sample and sample.name in num_pairs_by_sample:
        reads_stats['original_num_reads'] = num_pairs_by_sample[sample.name] * 2

    chrom_lengths = reference_data.get_chrom_lengths(genome=genome, fai_fpath=fai_fpath)
    if 'Y' in chrom_lengths or 'chrY' in chrom_lengths:
        reads_stats['gender'] = determine_sex(sample.work_dir, sample.bam, depth_stats['ave_depth'],
                                              genome, target.get_capture_bed())
        info()

    if 'bases_by_depth' in depth_stats:
        depth_stats['bases_within_threshs'], depth_stats['rates_within_threshs'] = calc_bases_within_threshs(
            depth_stats['bases_by_depth'],
            target_stats['target_size'] if not target.is_wgs else target_stats['reference_size'],
            depth_threshs)

        if depth_stats['median_depth'] > 0:
            depth_stats['wn_20_percent'] = calc_rate_within_normal(
                depth_stats['bases_by_depth'],
                depth_stats['median_depth'],
                target_stats['target_size'] if not target.is_wgs else target_stats['reference_size'])

    if target_stats['target_size']:
        target.bases_num = target_stats['target_size']
        target.fraction  = target_stats['target_fraction']
    else:
        target.bases_num = target_stats['reference_size']

    reads_stats['mapped_dedup'] = number_of_mapped_reads(sample.work_dir, sample.bam, dedup=True)

    if not target.is_wgs:
        reads_stats['mapped_dedup_on_target'] = number_mapped_reads_on_target(
            sample.work_dir, target.get_capture_bed().cut(range(3)).saveas().fn, sample.bam, dedup=True, target_name='target') or 0

        reads_stats['mapped_dedup_on_padded_target'] = number_mapped_reads_on_target(
            sample.work_dir, target.padded_bed_fpath, sample.bam, dedup=True, target_name='padded_target') or 0

    else:
        cds_bed = get_merged_cds(genome)
        info('Using the CDS reference BED to calc "reads on CDS"')
        reads_stats['mapped_dedup_on_exome'] = number_mapped_reads_on_target(
            sample.work_dir, cds_bed, sample.bam, dedup=True, target_name='exome') or 0

    return depth_stats, reads_stats, indels_stats


def _build_report(depth_stats, reads_stats, mm_indels_stats, sample, target,
                  depth_threshs, bed_padding, sample_num, is_debug=False, reannotate=False):
    report = SampleReport(sample, metric_storage=get_header_metric_storage(depth_threshs, is_wgs=target.bed_fpath is None, padding=bed_padding))

    def _add(_metric_name, _val, url=None):
        return report.add_record(_metric_name, _val, silent=(sample_num > 1 and not is_debug), url=url)

    _add('Qualimap', 'Qualimap', url=relpath(sample.qualimap_html_fpath, sample.dirpath))
    if reads_stats.get('gender') is not None:
        _add('Sex', reads_stats['gender'])

    debug('* General coverage statistics *')
    _add('Reads', reads_stats['total'])
    _add('Mapped reads', reads_stats['mapped'])
    # _add('Unmapped reads', reads_stats['totaAvgl'] - reads_stats['mapped'])
    percent_mapped = 1.0 * (reads_stats['mapped'] or 0) / reads_stats['total'] if reads_stats['total'] else None
    assert percent_mapped <= 1.0 or percent_mapped is None, str(percent_mapped)
    _add('Percentage of mapped reads', percent_mapped)
    # percent_unmapped = 1.0 * (reads_stats['total'] - reads_stats['mapped']) / reads_stats['total'] if reads_stats['total'] else None
    # assert percent_unmapped <= 1.0 or percent_unmapped is None, str(percent_unmapped)
    # _add('Percentage of unmapped reads', percent_unmapped)
    if reads_stats.get('mapped_paired') is not None:
        total_paired_reads_pecent = 1.0 * (reads_stats['mapped_paired'] or 0) / reads_stats['total'] if reads_stats['total'] else None
        assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
        _add('Properly paired mapped reads percent', total_paired_reads_pecent)
    # if reads_stats.get('paired') is not None:
    #     total_paired_reads_pecent = 1.0 * (reads_stats['paired'] or 0) / reads_stats['total'] if reads_stats['total'] else None
    #     assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
    #     _add('Properly paired reads percent', total_paired_reads_pecent)
    # if dedup_bam_stats:
    # dup_rate = 1 - (1.0 * dedup_bam_stats['mapped'] / bam_stats['mapped']) if bam_stats['mapped'] else None
    _add('Duplication rate', reads_stats['dup_rate'])
    # _add('Dedupped mapped reads', reads_stats['mapped'] - reads_stats[''])
    _add('Median GC', reads_stats['median_gc'])
    _add('Median insert size', reads_stats['median_ins_size'])

    debug()

    if not target.is_wgs:
        debug('* Target coverage statistics *')
        if target.original_bed_fpath:
            _add('Target', target.original_bed_fpath)
            if count_bed_cols(target.original_bed_fpath) == 3 or reannotate:
                _add('Ready target (clean, sorted and annotated)', target.capture_bed_fpath)
        else:
            _add('Target', target.capture_bed_fpath)
        _add('Bases in target', target.bases_num)
        _add('Percentage of reference', target.fraction)
        _add('Regions in target', target.regions_num)
        _add('Scope', 'targeted')
        _add('Genes in target', len(target.gene_keys_list))
    else:
        debug('* Genome coverage statistics *')
        _add('Target', 'whole genome')
        _add('Reference size', target.bases_num)
        _add('Scope', 'WGS')

    trg_type = 'target' if not target.is_wgs else 'genome'

    if 'bases_within_threshs' in depth_stats:
        bases_within_threshs = depth_stats['bases_within_threshs']
        v_covered_bases_in_targ = list(bases_within_threshs.items())[0][1]
        v_percent_covered_bases_in_targ = 1.0 * (v_covered_bases_in_targ or 0) / target.bases_num if target.bases_num else None
        assert v_percent_covered_bases_in_targ <= 1.0 or v_percent_covered_bases_in_targ is None, str(v_percent_covered_bases_in_targ)

        _add('Covered bases in ' + trg_type, v_covered_bases_in_targ)
        _add('Percentage of ' + trg_type + ' covered by at least 1 read', v_percent_covered_bases_in_targ)

    if not target.is_wgs:
        debug('Getting number of mapped reads on target...')
        # mapped_reads_on_target = number_mapped_reads_on_target(cnf, target_info.bed, bam_fpath)
        if 'mapped_dedup_on_target' in reads_stats:
            # _add('Reads mapped on target', reads_stats['mapped_on_target'])
            debug('Unique mapped reads on target: ' + str(reads_stats['mapped_dedup_on_target']))
            percent_mapped_dedup_on_target = 1.0 * reads_stats['mapped_dedup_on_target'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
            _add('Percentage of reads mapped on target', percent_mapped_dedup_on_target)
            assert percent_mapped_dedup_on_target <= 1.0 or percent_mapped_dedup_on_target is None, str(percent_mapped_dedup_on_target)

            percent_mapped_dedup_off_target = 1.0 * (reads_stats['mapped_dedup'] - reads_stats['mapped_dedup_on_target']) / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
            _add('Percentage of reads mapped off target', percent_mapped_dedup_off_target)
            assert percent_mapped_dedup_off_target <= 1.0 or percent_mapped_dedup_off_target is None, str(percent_mapped_dedup_off_target)

            percent_usable = 1.0 * reads_stats['mapped_dedup_on_target'] / reads_stats['total'] if reads_stats['total'] != 0 else None
            _add('Percentage of usable reads', percent_usable)
            assert percent_usable <= 1.0 or percent_usable is None, str(percent_usable)

        read_bases_on_targ = int(target.bases_num * depth_stats['ave_depth'])  # sum of all coverages
        _add('Read bases mapped on target', read_bases_on_targ)

        if 'mapped_dedup_on_padded_target' in reads_stats:
            # _add('Reads mapped on padded target', reads_stats['mapped_reads_on_padded_target'])
            percent_mapped_on_padded_target = 1.0 * reads_stats['mapped_dedup_on_padded_target'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] else None
            _add('Percentage of reads mapped on padded target', percent_mapped_on_padded_target)
            assert percent_mapped_on_padded_target <= 1.0 or percent_mapped_on_padded_target is None, str(percent_mapped_on_padded_target)

    elif 'mapped_dedup_on_exome' in reads_stats:
        # _add('Reads mapped on target', reads_stats['mapped_on_target'])
        percent_mapped_on_exome = 1.0 * reads_stats['mapped_dedup_on_exome'] / reads_stats['mapped_dedup'] if reads_stats['mapped_dedup'] != 0 else None
        if percent_mapped_on_exome:
            _add('Percentage of reads mapped on exome', percent_mapped_on_exome)
            assert percent_mapped_on_exome <= 1.0 or percent_mapped_on_exome is None, str(percent_mapped_on_exome)
            percent_mapped_off_exome = 1.0 - percent_mapped_on_exome
            _add('Percentage of reads mapped off exome ', percent_mapped_off_exome)

        percent_usable = 1.0 * reads_stats['mapped_dedup'] / reads_stats['total'] if reads_stats['total'] != 0 else None
        _add('Percentage of usable reads', percent_usable)
        assert percent_usable <= 1.0 or percent_usable is None, str(percent_usable)

    debug()
    _add('Mean ' + trg_type + ' coverage depth', depth_stats['ave_depth'])
    if 'original_num_reads' in reads_stats:
        _add('Original reads', reads_stats['original_num_reads'])
        times_downsampled = 1.0 * reads_stats['original_num_reads'] / reads_stats['total']
        est_full_cov = times_downsampled * depth_stats['ave_depth']
        _add('Estimated ' + trg_type + ' full coverage depth', est_full_cov)
    _add('Median ' + trg_type + ' coverage depth', depth_stats['median_depth'])
    if depth_stats['median_depth'] > 0:
        _add('Std. dev. of ' + trg_type + ' coverage depth', depth_stats['stddev_depth'])
    # _add('Minimal ' + trg_type + ' coverage depth', depth_stats['min_depth'])
    # _add('Maximum ' + trg_type + ' coverage depth', depth_stats['max_depth'])
    if 'wn_20_percent' in depth_stats:
        _add('Percentage of ' + trg_type + ' within 20% of med depth', depth_stats['wn_20_percent'])
        assert depth_stats['wn_20_percent'] <= 1.0 or depth_stats['wn_20_percent'] is None, str(depth_stats['wn_20_percent'])

    if 'bases_within_threshs' in depth_stats:
        for depth, bases in depth_stats['bases_within_threshs'].items():
            fraction_val = 1.0 * (bases or 0) / target.bases_num if target.bases_num else 0
            if fraction_val > 0:
                _add('Part of ' + trg_type + ' covered at least by ' + str(depth) + 'x', fraction_val)
            assert fraction_val <= 1.0 or fraction_val is None, str(fraction_val)
    debug()

    _add('Read mean length', reads_stats['ave_len'])
    _add('Read min length', reads_stats['min_len'])
    _add('Read max length', reads_stats['max_len'])
    _add('Mean Mapping Quality', mm_indels_stats['mean_mq'])
    _add('Mismatches', mm_indels_stats['mismatches'])
    _add('Insertions', mm_indels_stats['insertions'])
    _add('Deletions', mm_indels_stats['deletions'])
    _add('Homopolymer indels', mm_indels_stats['homo_indels'])

    debug()
    info('Saving reports...')
    report.save_json(sample.targqc_json_fpath)
    report.save_txt(sample.targqc_txt_fpath)
    report.save_html(sample.targqc_html_fpath, caption='Target coverage statistics for ' + sample.name)
    debug()
    debug('Saved to ' + dirname(report.txt_fpath))
    return report


def get_total_reads_number_from_fastqc(sample, fastqc_dirpath):
    fastqc_txt_fpaths = find_fastqc_txt(sample, fastqc_dirpath)
    if not fastqc_txt_fpaths:
        return
    num_reads = 0
    for fpath in fastqc_txt_fpaths:
        with open(fpath) as f_in:
            for line in f_in:
                if 'total sequences' in line.lower():
                    num_reads += int(line.strip().split('\t')[-1])
                    break
    return num_reads


def find_fastqc_txt(sample_name, fastqc_dirpath):
    l_fastqc_dirpath = join(fastqc_dirpath, sample_name + '_R1_fastqc')
    r_fastqc_dirpath = join(fastqc_dirpath, sample_name + '_R2_fastqc')
    fastqc_txt_fpaths = [join(l_fastqc_dirpath, 'fastqc_data.txt'), join(r_fastqc_dirpath, 'fastqc_data.txt')]
    if all(os.path.isfile(fpath) for fpath in fastqc_txt_fpaths):
        return fastqc_txt_fpaths
    else:
        return None

