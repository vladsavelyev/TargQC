# coding=utf-8

import math
from collections import OrderedDict
from os.path import dirname

import targqc.config as tc
from Utils.Region import Region, GeneInfo
from Utils.sambamba import sambamba_depth
from Utils.call_process import run
from Utils.file_utils import intermediate_fname, verify_file, safe_mkdir
from Utils.logger import info, debug
from Utils.reporting.reporting import ReportSection, Metric, MetricStorage, PerRegionSampleReport, \
    write_tsv_rows, write_txt_rows


def make_per_gene_report(work_dir, sample, target_bed, features_bed, gene_by_name_and_chrom, sample_avg_depth):
    info('-' * 70)
    info('Detailed exon-level report')
    per_gene_report = None
    if tc.reuse_intermediate and verify_file(sample.targqc_region_tsv, silent=True):
        debug(sample.targqc_region_tsv + ' exists, reusing')
        per_gene_report = _parse_report(sample, gene_by_name_and_chrom)
    else:
        if features_bed or target_bed:
            per_gene_report = _generate_regional_report_from_bam(work_dir, sample,
                target_bed, features_bed, gene_by_name_and_chrom, sample_avg_depth)
    return per_gene_report


def _parse_report(sample, gene_by_name_and_chrom):
    rep = PerRegionSampleReport()
    rep.txt_fpath = sample.targqc_region_txt
    rep.tsv_fpath = sample.targqc_region_tsv

    with open(sample.targqc_region_tsv) as f:
        for l in f:
            rep.rows.append(1)
            fs = l.strip().split('\t')

            if len(fs) < 13 or l.startswith('#'):
                continue
            chrom, start, end, size, gene_name, strand, feature, biotype, transcript_id, min_depth, avg_depth, std_dev, wn20ofmean = fs[:13]

            if gene_name not in ('.', '', None):
                region = Region(gene_name=gene_name, transcript_id=transcript_id, exon_num=None,
                                strand=strand, biotype=biotype, feature=feature, extra_fields=list(), chrom=chrom,
                                start=int(start) if start not in ('.', '', None) else None,
                                end=int(end) if end not in ('.', '', None) else None,
                                size=int(size) if size not in ('.', '', None) else None,
                                min_depth=float(min_depth) if min_depth not in ('.', '', None) else None,
                                avg_depth=float(avg_depth) if avg_depth not in ('.', '', None) else None,
                                std_dev=float(std_dev) if std_dev not in ('.', '', None) else None,
                                rate_within_normal=float(wn20ofmean) if wn20ofmean and wn20ofmean not in ('.', '', None) else None, )
                region.sample_name = gene_by_name_and_chrom[(gene_name, chrom)].sample_name
                depth_thresholds = tc.depth_thresholds
                rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)
                rates = fs[-(len(depth_thresholds)):]
                for i, t in enumerate(rates_within_threshs):
                    rates_within_threshs[t] = float(rates[i]) if rates[i] not in ('.', '', None) else None
                region.rates_within_threshs = rates_within_threshs
                if 'Capture' in feature:
                    gene_by_name_and_chrom[(gene_name, chrom)].add_amplicon(region)
                elif 'CDS' in feature or feature == 'Exon':
                    gene_by_name_and_chrom[(gene_name, chrom)].add_exon(region)
                else:
                    gene_by_name_and_chrom[(gene_name, chrom)].chrom = region.chrom
                    gene_by_name_and_chrom[(gene_name, chrom)].strand = region.strand
                    gene_by_name_and_chrom[(gene_name, chrom)].avg_depth = region.avg_depth
                    gene_by_name_and_chrom[(gene_name, chrom)].min_depth = region.min_depth
                    gene_by_name_and_chrom[(gene_name, chrom)].rates_within_threshs = region.rates_within_threshs
    return rep


def get_detailed_metric_storage(depth_threshs):
    return MetricStorage(
        general_section=ReportSection(metrics=[
            Metric('Sample'),
        ]),
        sections=[ReportSection(metrics=[
            Metric('Chr'),
            Metric('Start'),
            Metric('End'),
            Metric('Size'),
            Metric('Gene'),
            Metric('Strand'),
            Metric('Feature'),
            Metric('Biotype'),
            Metric('Transcript'),
            Metric('Min depth'),
            Metric('Ave depth'),
            Metric('Std dev', description='Coverage depth standard deviation'),
            Metric('W/n 20% of ave depth', description='Percentage of the region that lies within 20% of an avarage depth.', unit='%'),
            # Metric('Norm depth', description='Ave region depth devided by median depth of sample'),
        ] + [
            Metric('{}x'.format(thresh), description='Bases covered by at least {} reads'.format(thresh), unit='%')
            for thresh in depth_threshs
        ])]
    )



def _sambamba_depth_to_regions(sambamba_depth_output_fpath, sample_name, target_type, depth_thresholds):
    read_count_col = None
    mean_cov_col = None
    min_depth_col = None
    std_dev_col = None

    regions = []
    #####################################
    #####################################
    info('Reading coverage statistics...')
    with open(sambamba_depth_output_fpath) as sambabma_depth_file:
        total_regions_count = 0
        for line in sambabma_depth_file:
            fs = line.replace('\n', '').split('\t')
            if line.startswith('#'):
                fs = line.split('\t')
                read_count_col = fs.index('readCount')
                mean_cov_col = fs.index('meanCoverage')
                min_depth_col = fs.index('minDepth') if 'minDepth' in fs else None
                std_dev_col = fs.index('stdDev') if 'stdDev' in fs else None
                continue
            chrom = fs[0]
            start, end = map(int, fs[1:3])
            region_size = end - start
            gene_name = fs[3] if read_count_col != 3 else None
            ave_depth = float(fs[mean_cov_col])
            min_depth = int(fs[min_depth_col]) if min_depth_col is not None else '.'
            std_dev = float(fs[std_dev_col]) if std_dev_col is not None else '.'
            rates_within_threshs = fs[(std_dev_col or mean_cov_col) + 1:-1]

            extra_fields = tuple(fs[4:read_count_col]) if read_count_col > 4 else ()

            region = Region(
                sample_name=sample_name, chrom=chrom,
                start=start, end=end, size=region_size,
                avg_depth=ave_depth,
                gene_name=gene_name, extra_fields=extra_fields)
            regions.append(region)

            region.rates_within_threshs = OrderedDict((depth, float(rate) / 100.0) for (depth, rate) in zip(depth_thresholds, rates_within_threshs))
            region.min_depth = min_depth
            region.std_dev = std_dev

            if target_type == 'amplicons':
                region.feature = 'Capture'
            else:
                if extra_fields:
                    region.exon_num = extra_fields[0]
                    if len(extra_fields) >= 2:
                        region.strand = extra_fields[1]
                    if len(extra_fields) >= 3:
                        region.feature = extra_fields[2]
                    else:
                        region.feature = 'CDS'
                    if len(extra_fields) >= 4:
                        region.biotype = extra_fields[3]
                    if len(extra_fields) >= 5:
                        region.transcript_id = extra_fields[4]

            total_regions_count += 1
            if total_regions_count > 0 and total_regions_count % 10000 == 0:
                info('  Processed {0:,} regions'.format(total_regions_count))
    info('Total regions: ' + str(len(regions)))

    # #####################################
    # #####################################
    # info('Second round of sambamba depth - calculating depth within 20% bounds of average depth')
    # sambamba_depth_output_fpath = sambamba_depth(cnf, bed, bam, depth_thresholds=[int()])
    # if not sambamba_depth_output_fpath:
    #     continue
    #
    # info('Adding rates within normal...')
    # with open(sambamba_depth_output_fpath) as sambabma_depth_file:
    #     total_regions_count = 0
    #     for region, line in zip(regions, (l for l in sambabma_depth_file if not l.startswith('#'))):
    #         line_tokens = line.replace('\n', '').split()
    #         rate_within_low_bound = line_tokens[std_dev_col + 1]
    #         rate_within_higher_bound = line_tokens[std_dev_col + 2]
    #         regions.rate_within_normal = rate_within_low_bound - rate_within_higher_bound
    #
    #         total_regions_count += 1
    #         if total_regions_count > 0 and total_regions_count % 10000 == 0:
    #              info('  Processed {0:,} regions'.format(total_regions_count))
    #     info('Processed {0:,} regions'.format(total_regions_count))

    return regions


def _generate_regional_report_from_bam(work_dir, sample, target_bed, features_bed, gene_by_name_and_chrom, avg_depth):
    depth_thresholds = tc.depth_thresholds
    if avg_depth:
        key_gene_cov_threshold = max(1, int(avg_depth / 2))
        depth_thresholds.append(key_gene_cov_threshold)
        depth_thresholds.sort()

    #####################################
    #####################################
    ready_to_report_genes = []
    ready_to_report_set = set()

    debug('Filtering features BED to have only CDS and Exon features')
    exons_and_cds_features = intermediate_fname(work_dir, features_bed, 'nogenes')
    run('egrep -w "Exon|CDS" ' + features_bed, output_fpath=exons_and_cds_features, reuse=tc.reuse_intermediate)

    # TODO:
    # Process both BED files in parallel, and merge by gene (amplicons, exons, gene-summary).
    # Then run sambamba on the merged set.
    # Then interage the sambamba results and re-format it, and write TSV and TXT in parallel.

    # debug('Interating in parallel')
    # target_and_exons_fpath = join(work_dir, 'target_and_exons.bed')
    # amplicons_f = open(target_bed) if target_bed else None
    # exons_f = open(features_bed) if features_bed else None
    #
    # with open(target_and_exons_fpath, 'w') as out:
    #     # TODO: make sure the order of cur_gene. Handle overlapping genes.
    #     for cur_gene in gene_by_name_and_chrom.values():
    #         amplicon_l = None
    #         exon_l = None
    #
    #         if amplicon_l:  # first line after which the previous "while True" was interrupted - this line belongs to the next gene
    #             assert amplicon_l.split('\t')[3] == cur_gene.gene_name
    #             out.write(amplicon_l)
    #
    #         # iterate over regions until meat another gene
    #         if amplicons_f:
    #             while True:
    #                 amplicon_l = next(amplicons_f, None)
    #                 if not amplicon_l: break
    #                 if amplicon_l.startswith('#'): continue
    #                 chrom, _, _, gname = amplicon_l.strip().split('\t')[:4]
    #                 if (gname, chrom) != (cur_gene.gene_name, cur_gene.chrom) and gname != '.': break  # another gene
    #                 out.write(amplicon_l)
    #
    #         if exon_l:  # first line after which the previous "while True" was interrupted - this line belongs to the next gene
    #             assert exon_l.split('\t')[3] == cur_gene.gene_name
    #             out.write(exon_l)
    #
    #         if exons_f and cur_gene.gene_name not in genes_not_in_refseq:  # gene will be found!
    #             while True:
    #                 exon_l = next(exons_f, None)
    #                 if not exon_l: break
    #                 if exon_l.startswith('#'): continue
    #                 chrom, _, _, gname, _, _, feature = exon_l.strip().split('\t')[:7]
    #                 if (gname, chrom) != (cur_gene.gene_name, cur_gene.chrom): break
    #                 if feature in ['CDS', 'Exon']:
    #                     out.write(exon_l)
    #
    # if target_bed: amplicons_f.close()
    # if features_bed: exons_f.close()

    # debug('Saved mixed BED into ' + target_and_exons_fpath)
    # sys.exit(1)

    for (bed_fpath, target_type) in zip([target_bed, exons_and_cds_features], ['amplicons', 'exons']):  # features are canonical
        if not bed_fpath:
            continue
        info()
        info('Calculating coverage statistics for ' + ('CDS and ncRNA exons...' if target_type == 'exons' else 'the regions in the target BED file...'))
        sambamba_depth_output_fpath = sambamba_depth(work_dir, bed_fpath, sample.bam, depth_thresholds, reuse=tc.reuse_intermediate)
        regions = _sambamba_depth_to_regions(sambamba_depth_output_fpath, sample.name, target_type, depth_thresholds)

        #####################################
        #####################################
        info('Filling regions stats into the gene objects...')
        cur_unannotated_gene = None
        total_regions_count = 0
        for region in regions:
            if region.feature == 'Capture':
                if region.gene_name != '.':
                    cur_unannotated_gene = None
                    gene = gene_by_name_and_chrom[(region.gene_name, region.chrom)]
                    if (gene.gene_name, gene.chrom) not in ready_to_report_set:
                        ready_to_report_genes.append(gene)
                        ready_to_report_set.add((gene.gene_name, gene.chrom))
                    gene.add_amplicon(region)
                else:
                    if cur_unannotated_gene is None:
                        cur_unannotated_gene = GeneInfo(sample_name=region.sample_name,
                            gene_name=region.gene_name, chrom=region.chrom, feature='NotAnnotatedSummary')
                        ready_to_report_genes.append(cur_unannotated_gene)
                    cur_unannotated_gene.add_amplicon(region)

            else:
                gene = gene_by_name_and_chrom[(region.gene_name, region.chrom)]
                gene.add_exon(region)
                if not target_bed:  # in case if only reporting based on features_bed
                    if (gene.gene_name, gene.chrom) not in ready_to_report_set:
                        ready_to_report_genes.append(gene)
                        ready_to_report_set.add((gene.gene_name, gene.chrom))

            # row = [region.chrom, region.start, region.end, region.get_size(), region.gene_name, region.strand,
            #        region.feature, region.biotype, region.transcript_id, region.min_depth, region.avg_depth, region.std_dev,
            #        region.rate_within_normal]
            # row = [Metric.format_value(val, human_readable=True) for val in row]
            # rates = [Metric.format_value(val, unit='%', human_readable=True) for val in region.rates_within_threshs.values()]
            # row.extend(rates)
            # col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

            total_regions_count += 1
            if total_regions_count > 0 and total_regions_count % 10000 == 0:
                info('  Processed {0:,} regions'.format(total_regions_count))
        info('Processed {0:,} regions'.format(total_regions_count))

    #####################################
    #####################################
    report = PerRegionSampleReport(sample=sample, metric_storage=get_detailed_metric_storage(depth_thresholds))
    # report.add_record('Sample', sample.name)
    report.txt_fpath = sample.targqc_region_txt
    report.tsv_fpath = sample.targqc_region_tsv

    debug('Arranging regions to report...')
    regions = []
    for g in ready_to_report_genes:
        for a in g.get_amplicons():
            regions.append(a)
        for e in g.get_exons():
            regions.append(e)
        if g.get_exons():
            process_gene(g, depth_thresholds)
            regions.append(g)

    debug('Preparting report rows...')
    for reg in regions:
        r = report.add_row()
        r.add_record('Chr', reg.chrom)
        r.add_record('Start', reg.start)
        r.add_record('End', reg.end)
        r.add_record('Size', reg.get_size())
        r.add_record('Gene', reg.gene_name)
        r.add_record('Strand', reg.strand)
        r.add_record('Feature', reg.feature)
        r.add_record('Biotype', reg.biotype)
        r.add_record('Transcript', reg.transcript_id)
        r.add_record('Min depth', reg.min_depth)
        r.add_record('Ave depth', reg.avg_depth)
        r.add_record('Std dev', reg.std_dev)
        r.add_record('W/n 20% of ave depth', reg.rate_within_normal)
        for ths in depth_thresholds:
            r.add_record('{}x'.format(ths), reg.rates_within_threshs.get(ths) if reg.rates_within_threshs else None)

    safe_mkdir(dirname(report.tsv_fpath))
    debug('Flattening TSV records...')
    header_rows, flat_rows = report.flatten(None, human_readable=False)
    debug('Writing TSV...')
    write_tsv_rows((header_rows, flat_rows), report.tsv_fpath)
    debug('Flattening TXT records...')
    header_rows, flat_rows = report.flatten(None, human_readable=True)
    debug('Writing TXT...')
    write_txt_rows((header_rows, flat_rows), report.txt_fpath)

    # un_annotated_summary_region = next((g for g in gene_by_name_and_chrom.values() if g.gene_name == '.'), None)
    # if un_annotated_summary_region and un_annotated_amplicons:
    #     un_annotated_summary_region.feature = 'NotAnnotatedSummary'
    #     for ampl in un_annotated_amplicons:
    #         ampl.gene_name = un_annotated_summary_region.gene_name
    #         un_annotated_summary_region.add_amplicon(ampl)
    #         add_region_to_report(report, ampl, depth_thresholds)
    #     add_region_to_report(report, un_annotated_summary_region, depth_thresholds)

    # report.save_txt(sample.targqc_region_txt)
    # report.save_tsv(sample.targqc_region_tsv)
    info('Regions (total ' + str(len(report.rows)) + ') saved into:')
    info('  ' + report.txt_fpath)
    return report


def process_gene(gene, depth_thresholds):
    gene.rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)
    if gene.size == 0:
        gene.start = None
        gene.end = None
        return
    exons = gene.get_exons()
    if not exons:
        return

    total_depth = sum(e.avg_depth * e.size for e in exons)
    gene.size = sum(e.size for e in exons)
    gene.avg_depth = total_depth / gene.size
    if exons[0].std_dev not in ['.', '', None]:
        sum_of_sq_var = sum((((e.avg_depth - e.std_dev) - gene.avg_depth) ** 2 + ((e.avg_depth + e.std_dev) - gene.avg_depth) ** 2) * e.size for e in exons)
        gene.std_dev = math.sqrt(sum_of_sq_var / 2 / float(gene.size))
    for t in depth_thresholds:
        total_rate = sum(e.rates_within_threshs[t] * e.size for e in exons)
        rate = total_rate / gene.size
        gene.rates_within_threshs[t] = rate


# def add_region_to_report(report, region, depth_threshs, file_to_write=None, col_widths=None):
#     if file_to_write:
#         report.rows.append(1)
#         rep_region = Row(parent_report=report)
#     else:
#         rep_region = report.add_row()
#     rep_region.add_record('Chr', region.chrom)
#     rep_region.add_record('Start', region.start)
#     rep_region.add_record('End', region.end)
#     rep_region.add_record('Size', region.get_size())
#     rep_region.add_record('Gene', region.gene_name)
#     rep_region.add_record('Strand', region.strand)
#     rep_region.add_record('Feature', region.feature)
#     rep_region.add_record('Biotype', region.biotype)
#     rep_region.add_record('Transcript', region.transcript_id)
#     rep_region.add_record('Min depth', region.min_depth)
#     rep_region.add_record('Ave depth', region.avg_depth)
#     rep_region.add_record('Std dev', region.std_dev)
#     rep_region.add_record('W/n 20% of ave depth', region.rate_within_normal)
#
#     if region.rates_within_threshs is None:
#         warn('Error: no rates_within_threshs for ' + str(region))
#     for thresh in depth_threshs:
#         rep_region.add_record('{}x'.format(thresh), region.rates_within_threshs.get(thresh) if region.rates_within_threshs else None)
#
#     if file_to_write:
#         for fpath in fpaths_to_write:
#             human_readable = fpath.endswith('txt')
#             flat_row = []
#             for m in report.metric_storage.get_metrics(None, skip_general_section=True):
#                 rec = BaseReport.find_record(rep_region.records, m.name)
#                 if rec:
#                     flat_row.append(rec.format(human_readable=human_readable))
#
#             with open(fpath, 'a') as out:
#                 if fpath.endswith('tsv'):
#                     out.write('\t'.join([val for val in flat_row]) + '\n')
#                 else:
#                     col_widths = col_widths or repeat(0)
#                     for val, w in izip(flat_row, col_widths):
#                         out.write(val + (' ' * (w - len(val) + 2)))
#                     out.write('\n')
#
#     return col_widths

# def _bases_by_depth(depth_vals, depth_thresholds):
#     bases_by_min_depth = {depth: 0 for depth in depth_thresholds}
#
#     for depth_value in depth_vals:
#         for threshold in depth_thresholds:
#             if depth_value >= threshold:
#                 bases_by_min_depth[threshold] += 1
#
#         return [1.0 * bases_by_min_depth[thres] / len(depth_vals) if depth_vals else 0
#                 for thres in depth_thresholds]


