# coding=utf-8

from collections import defaultdict
from os.path import isfile, join

from ngs_utils import reference_data
from ngs_utils.bed_utils import count_bed_cols
from ngs_utils.sambamba import sambamba_depth
from ngs_utils.call_process import run
from ngs_utils.file_utils import intermediate_fname, verify_file, file_transaction, can_reuse
from ngs_utils.logger import info, debug
from ngs_utils.utils import OrderedDefaultDict

import ensembl as ebl


def make_region_reports(view, work_dir, samples, target, genome, depth_thresholds):
    bed_fpath = target.bed_fpath or target.wgs_bed_fpath

    if all(can_reuse(s.targqc_region_tsv, [s.bam, bed_fpath]) for s in samples):
        debug('All region reports exist, reusing')
        return [s.targqc_region_tsv for s in samples]

    info('Calculating coverage statistics for CDS and exon regions from RefSeq...')

    depth_thresholds_by_sample = dict()
    for s in samples:
        depth_thresholds_by_sample[s.name] = depth_thresholds

    debug()
    debug('Running sambamba...')
    sambamba_depth_output_fpaths = view.run(sambamba_depth,
        [[s.work_dir, bed_fpath, s.bam, depth_thresholds_by_sample[s.name], None, s.name]
         for s in samples])
    assert len(sambamba_depth_output_fpaths) == len(samples), \
        'Number of sambamba results = ' + str(len(sambamba_depth_output_fpaths)) + \
        ' which is less then the number of samples ' + str(len(samples))

    debug()
    debug('Parsing sambamba results and writing results...')
    view.run(_proc_sambamba_depth,
        [[sambamba_output_fpath, s.targqc_region_tsv, s.name, depth_thresholds_by_sample[s.name]]
         for sambamba_output_fpath, s in zip(sambamba_depth_output_fpaths, samples)])

    info('Done.')
    return [s.targqc_region_tsv for s in samples]


def _proc_sambamba_depth(sambamba_depth_output_fpath, output_fpath, sample_name, depth_thresholds):
    read_count_col = None
    mean_cov_col = None
    median_cov_col = None
    min_depth_col = None
    std_dev_col = None
    wn_20_pcnt_col = None

    regions_by_genekey = defaultdict(list)
    #####################################
    #####################################
    if can_reuse(output_fpath, sambamba_depth_output_fpath):
        return output_fpath

    debug('Reading coverage statistics and writing regions to ' + output_fpath)

    def write_line(f, fields):
        f.write('\t'.join(fields) + '\n')

    with file_transaction(None, output_fpath) as tx:
        with open(sambamba_depth_output_fpath) as sambabma_depth_file, open(tx, 'w') as out:
            total_regions_count = 0
            for line in sambabma_depth_file:
                fs = line.strip('\n').split('\t')
                if line.startswith('#'):
                    fs = line.split('\t')
                    read_count_col = fs.index('readCount')
                    mean_cov_col = fs.index('meanCoverage')
                    median_cov_col = fs.index('medianCoverage') if 'medianCoverage' in fs else None
                    min_depth_col = fs.index('minDepth') if 'minDepth' in fs else None
                    std_dev_col = fs.index('stdDev') if 'stdDev' in fs else None
                    wn_20_pcnt_col = fs.index('percentWithin20PercentOfMedian') if 'percentWithin20PercentOfMedian' in fs else None

                    write_line(out, [
                        '#Chr',
                        'Start',
                        'End',
                        'Size',
                        'Gene',
                        'Exon',
                        'Strand',
                        'Feature',
                        'Biotype',
                        'Transcript',
                        'Tx overlap',
                        'Exome overlap',
                        'Min depth',
                        'Avg depth',
                        'Median depth',
                        'Std dev',
                        'W/n 20% of median',
                    ] + ['{}x'.format(ths) for ths in depth_thresholds])
                    continue

                chrom = fs[0]
                start, end = int(fs[1]), int(fs[2])
                region_size = end - start
                gene_name = fs[ebl.BedCols.GENE] if read_count_col != ebl.BedCols.GENE else '.'
                exon = fs[ebl.BedCols.EXON]
                strand = fs[ebl.BedCols.STRAND]
                feature = fs[ebl.BedCols.FEATURE]
                biotype = fs[ebl.BedCols.BIOTYPE]
                transcript = fs[ebl.BedCols.ENSEMBL_ID]
                transcript_overlap = fs[ebl.BedCols.TX_OVERLAP_PERCENTAGE]
                exome_overlap = fs[ebl.BedCols.EXON_OVERLAPS_PERCENTAGE]
                avg_depth = float(fs[mean_cov_col])
                min_depth = int(fs[min_depth_col]) if min_depth_col is not None else '.'
                std_dev = float(fs[std_dev_col]) if std_dev_col is not None else '.'
                median_depth = int(fs[median_cov_col]) if median_cov_col is not None else '.'
                rate_within_normal = float(fs[wn_20_pcnt_col]) if wn_20_pcnt_col is not None else '.'
                last_cov_col = max(mean_cov_col or 0, median_cov_col or 0, std_dev_col or 0, wn_20_pcnt_col or 0)
                rates_within_threshs = fs[last_cov_col + 1:-1]

                write_line(out, map(str, [
                        chrom,
                        start,
                        end,
                        region_size,
                        gene_name,
                        exon,
                        strand,
                        feature,
                        biotype,
                        transcript,
                        ((transcript_overlap + '%') if transcript_overlap not in ['', None, '.'] else '.'),
                        ((exome_overlap + '%') if exome_overlap not in ['', None, '.'] else '.'),
                        min_depth,
                        avg_depth,
                        median_depth,
                        std_dev,
                        rate_within_normal,
                    ] + rates_within_threshs))

                total_regions_count += 1
                if total_regions_count > 0 and total_regions_count % 10000 == 0:
                    debug('  Processed {0:,} regions'.format(total_regions_count))
        debug('Total regions: ' + str(len(regions_by_genekey)))
    return output_fpath

# def _get_values_from_row(fields, cols):
#     return [fields[col] if col else None for col in cols]


# def _parse_report(sample, gene_by_name_and_chrom):
#     rep = PerRegionSampleReport()
#     rep.txt_fpath = sample.targqc_region_txt
#     rep.tsv_fpath = sample.targqc_region_tsv
#
#     with open(sample.targqc_region_tsv) as f:
#         min_depth_col = None
#         median_depth_col = None
#         avg_depth_col = None
#         std_dev_col = None
#         wn_20_pcnt_col = None
#         for l in f:
#             rep.rows.append(1)
#             fs = l.strip().split('\t')
#
#             if l.startswith('#'):
#                 if len(fs) > 9:
#                     min_depth_col = fs.index('Min depth') if 'Min depth' in fs else min_depth_col
#                     median_depth_col = fs.index('Median depth') if 'Median depth' in fs else median_depth_col
#                     avg_depth_col = fs.index('Ave depth') if 'Ave depth' in fs else avg_depth_col
#                     std_dev_col = fs.index('Std dev') if 'Std dev' in fs else std_dev_col
#                     wn_20_pcnt_col = fs.index('W/n 20% of median depth') if 'W/n 20% of median depth' in fs else wn_20_pcnt_col
#                 continue
#             chrom, start, end, size, gene_name, strand, feature, biotype, transcript_id = fs[:9]
#             min_depth, median_depth, avg_depth, std_dev, wn20ofmedian = \
#                 _get_values_from_row(fs, [min_depth_col, median_depth_col, avg_depth_col, std_dev_col, wn_20_pcnt_col])
#
#             if gene_name not in ('.', '', None):
#                 region = Region(gene_name=gene_name, transcript_id=transcript_id, exon_num=None,
#                                 strand=strand, biotype=biotype, feature=feature, extra_fields=list(), chrom=chrom,
#                                 start=int(start) if start not in ('.', '', None) else None,
#                                 end=int(end) if end not in ('.', '', None) else None,
#                                 size=int(size) if size not in ('.', '', None) else None,
#                                 min_depth=float(min_depth) if min_depth not in ('.', '', None) else None,
#                                 median_depth=float(median_depth) if median_depth not in ('.', '', None) else None,
#                                 avg_depth=float(avg_depth) if avg_depth not in ('.', '', None) else None,
#                                 std_dev=float(std_dev) if std_dev not in ('.', '', None) else None,
#                                 rate_within_normal=float(wn20ofmedian) if wn20ofmedian and wn20ofmedian not in ('.', '', None) else None, )
#                 region.sample_name = gene_by_name_and_chrom[(gene_name, chrom)].sample_name
#                 depth_thresholds = cfg.depth_thresholds
#                 rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)
#                 rates = fs[-(len(depth_thresholds)):]
#                 for i, t in enumerate(rates_within_threshs):
#                     rates_within_threshs[t] = float(rates[i]) if rates[i] not in ('.', '', None) else None
#                 region.rates_within_threshs = rates_within_threshs
#                 if 'Capture' in feature:
#                     gene_by_name_and_chrom[(gene_name, chrom)].add_amplicon(region)
#                 elif 'CDS' in feature or feature == 'Exon':
#                     gene_by_name_and_chrom[(gene_name, chrom)].add_exon(region)
#                 else:
#                     gene_by_name_and_chrom[(gene_name, chrom)].chrom = region.chrom
#                     gene_by_name_and_chrom[(gene_name, chrom)].strand = region.strand
#                     gene_by_name_and_chrom[(gene_name, chrom)].avg_depth = region.avg_depth
#                     gene_by_name_and_chrom[(gene_name, chrom)].min_depth = region.min_depth
#                     gene_by_name_and_chrom[(gene_name, chrom)].median_depth = region.median_depth
#                     gene_by_name_and_chrom[(gene_name, chrom)].rates_within_threshs = region.rates_within_threshs
#     return rep


# def get_detailed_metric_storage(depth_threshs):
#     return MetricStorage(
#         general_section=ReportSection(metrics=[
#             Metric('Sample'),
#         ]),
#         sections=[ReportSection(metrics=[
#             Metric('Chr'),
#             Metric('Start'),
#             Metric('End'),
#             Metric('Size'),
#             Metric('Gene'),
#             Metric('Strand'),
#             Metric('Feature'),
#             Metric('Biotype'),
#             Metric('Transcript'),
#             Metric('Min depth'),
#             Metric('Median depth'),
#             Metric('Ave depth'),
#             Metric('Std dev', description='Coverage depth standard deviation'),
#             Metric('W/n 20% of median depth', description='Percentage of the region that lies within 20% of an median depth.', unit='%'),
#             # Metric('Norm depth', description='Ave region depth devided by median depth of sample'),
#         ] + [
#             Metric('{}x'.format(thresh), description='Bases covered by at least {} reads'.format(thresh), unit='%')
#             for thresh in depth_threshs
#         ])]
#     )


# def _sambamba_depth_to_regions(sambamba_depth_output_fpath, sample_name, target_type, depth_thresholds):
#     read_count_col = None
#     mean_cov_col = None
#     median_cov_col = None
#     min_depth_col = None
#     std_dev_col = None
#     wn_20_pcnt_col = None
#
#     regions_by_genekey = defaultdict(list)
#     #####################################
#     #####################################
#     debug('Reading coverage statistics...')
#     with open(sambamba_depth_output_fpath) as sambabma_depth_file:
#         total_regions_count = 0
#         for line in sambabma_depth_file:
#             fs = line.replace('\n', '').split('\t')
#             if line.startswith('#'):
#                 fs = line.split('\t')
#                 read_count_col = fs.index('readCount')
#                 mean_cov_col = fs.index('meanCoverage')
#                 median_cov_col = fs.index('medianCoverage') if 'medianCoverage' in fs else None
#                 min_depth_col = fs.index('minDepth') if 'minDepth' in fs else None
#                 std_dev_col = fs.index('stdDev') if 'stdDev' in fs else None
#                 wn_20_pcnt_col = fs.index('W/n 20% of median depth') if 'W/n 20% of median depth' in fs else None
#                 continue
#             chrom = fs[0]
#             start, end = map(int, fs[1:3])
#             region_size = end - start
#             gene_name = fs[3] if read_count_col != 3 else None
#             ave_depth = float(fs[mean_cov_col])
#             min_depth = int(fs[min_depth_col]) if min_depth_col is not None else '.'
#             std_dev = float(fs[std_dev_col]) if std_dev_col is not None else '.'
#             median_depth = int(fs[median_cov_col]) if median_cov_col is not None else '.'
#             rate_within_normal = float(fs[wn_20_pcnt_col]) if wn_20_pcnt_col is not None else '.'
#             last_cov_col = max(mean_cov_col, median_cov_col, std_dev_col, wn_20_pcnt_col)
#             rates_within_threshs = fs[last_cov_col + 1:-1]
#
#             extra_fields = tuple(fs[4:read_count_col]) if read_count_col > 4 else ()
#
#             region = Region(
#                 sample_name=sample_name, chrom=chrom,
#                 start=start, end=end, size=region_size,
#                 avg_depth=ave_depth, median_depth=median_depth,
#                 gene_name=gene_name, extra_fields=extra_fields,
#                 rate_within_normal=rate_within_normal)
#             regions_by_genekey[(gene_name, chrom)].append(region)
#
#             region.rates_within_threshs = OrderedDict((depth, float(rate) / 100.0) for (depth, rate) in zip(depth_thresholds, rates_within_threshs))
#             region.min_depth = min_depth
#             region.std_dev = std_dev
#
#             if target_type == 'amplicons':
#                 region.feature = 'Capture'
#             else:
#                 if extra_fields:
#                     region.exon_num = extra_fields[0]
#                     if len(extra_fields) >= 2:
#                         region.strand = extra_fields[1]
#                     if len(extra_fields) >= 3:
#                         region.feature = extra_fields[2]
#                     else:
#                         region.feature = 'CDS'
#                     if len(extra_fields) >= 4:
#                         region.biotype = extra_fields[3]
#                     if len(extra_fields) >= 5:
#                         region.transcript_id = extra_fields[4]
#
#             total_regions_count += 1
#             if total_regions_count > 0 and total_regions_count % 10000 == 0:
#                 debug('  Processed {0:,} regions'.format(total_regions_count))
#     debug('Total regions: ' + str(len(regions_by_genekey)))
#
#     # #####################################
#     # #####################################
#     # info('Second round of sambamba depth - calculating depth within 20% bounds of average depth')
#     # sambamba_depth_output_fpath = sambamba_depth(cnf, bed, bam, depth_thresholds=[int()])
#     # if not sambamba_depth_output_fpath:
#     #     continue
#     #
#     # info('Adding rates within normal...')
#     # with open(sambamba_depth_output_fpath) as sambabma_depth_file:
#     #     total_regions_count = 0
#     #     for region, line in zip(regions, (l for l in sambabma_depth_file if not l.startswith('#'))):
#     #         line_tokens = line.replace('\n', '').split()
#     #         rate_within_low_bound = line_tokens[std_dev_col + 1]
#     #         rate_within_higher_bound = line_tokens[std_dev_col + 2]
#     #         regions.rate_within_normal = rate_within_low_bound - rate_within_higher_bound
#     #
#     #         total_regions_count += 1
#     #         if total_regions_count > 0 and total_regions_count % 10000 == 0:
#     #              info('  Processed {0:,} regions'.format(total_regions_count))
#     #     info('Processed {0:,} regions'.format(total_regions_count))
#
#     return regions_by_genekey


# def _write_regions(refseq_regions_by_genekey, target_regions_by_genekey, genekeys_list,
#                    output_fpath, depth_thresholds, reuse=False):
#     if reuse and isfile(output_fpath) and verify_file(output_fpath):
#         debug('Regional report ' + output_fpath + ' exists, reusing')
#         return
#
#     debug('Writing regions to ' + output_fpath)
#     with file_transaction(None, output_fpath) as tx:
#         with open(tx, 'w') as out:
#             out.write('#' + '\t'.join([
#                 'Chr',
#                 'Start',
#                 'End',
#                 'Size',
#                 'Gene',
#                 'Strand',
#                 'Feature',
#                 'Biotype',
#                 'Transcript',
#                 'Min depth',
#                 'Median depth',
#                 'Std dev',
#                 'W/n 20% of median',
#             ] + ['{}x'.format(ths) for ths in depth_thresholds]) + '\n')
#
#             for (gene_name, chrom) in genekeys_list:
#                 for r in refseq_regions_by_genekey[(gene_name, chrom)] + \
#                          target_regions_by_genekey[(gene_name, chrom)]:
#                     out.write('\t'.join(map(str, [
#                         r.chrom,
#                         r.start,
#                         r.end,
#                         r.get_size(),
#                         r.gene_name,
#                         r.strand or '.',
#                         r.feature,
#                         r.biotype or '.',
#                         r.transcript_id or '.',
#                         r.min_depth,
#                         r.median_depth,
#                         r.avg_depth,
#                         r.std_dev,
#                         r.rate_within_normal,
#                     ] + [r.rates_within_threshs.get(ths) for ths in depth_thresholds])) + '\n')


# def _prev_generate_regional_reports_from_bam(view, work_dir, samples, target, features_bed):
#     # if reuse and verify_file(sample.targqc_region_tsv, silent=True):
#     #     debug(sample.targqc_region_tsv + ' exists, reusing')
#     #     return _parse_report(sample, target)
#
#     #####################################
#     #####################################
#     ready_to_report_genes = []
#     ready_to_report_set = set()
#
#     # debug('Interating in parallel')
#     # target_and_exons_fpath = join(work_dir, 'target_and_exons.bed')
#     # amplicons_f = open(target_bed) if target_bed else None
#     # exons_f = open(features_bed) if features_bed else None
#     #
#     # with open(target_and_exons_fpath, 'w') as out:
#     #     for cur_gene in gene_by_name_and_chrom.values():
#     #         amplicon_l = None
#     #         exon_l = None
#     #
#     #         if amplicon_l:  # first line after which the previous "while True" was interrupted - this line belongs to the next gene
#     #             assert amplicon_l.split('\t')[3] == cur_gene.gene_name
#     #             out.write(amplicon_l)
#     #
#     #         # iterate over regions until meat another gene
#     #         if amplicons_f:
#     #             while True:
#     #                 amplicon_l = next(amplicons_f, None)
#     #                 if not amplicon_l: break
#     #                 if amplicon_l.startswith('#'): continue
#     #                 chrom, _, _, gname = amplicon_l.strip().split('\t')[:4]
#     #                 if (gname, chrom) != (cur_gene.gene_name, cur_gene.chrom) and gname != '.': break  # another gene
#     #                 out.write(amplicon_l)
#     #
#     #         if exon_l:  # first line after which the previous "while True" was interrupted - this line belongs to the next gene
#     #             assert exon_l.split('\t')[3] == cur_gene.gene_name
#     #             out.write(exon_l)
#     #
#     #         if exons_f and cur_gene.gene_name not in genes_not_in_refseq:  # gene will be found!
#     #             while True:
#     #                 exon_l = next(exons_f, None)
#     #                 if not exon_l: break
#     #                 if exon_l.startswith('#'): continue
#     #                 chrom, _, _, gname, _, _, feature = exon_l.strip().split('\t')[:7]
#     #                 if (gname, chrom) != (cur_gene.gene_name, cur_gene.chrom): break
#     #                 if feature in ['CDS', 'Exon']:
#     #                     out.write(exon_l)
#     #
#     # if target_bed: amplicons_f.close()
#     # if features_bed: exons_f.close()
#
#     # debug('Saved mixed BED into ' + target_and_exons_fpath)
#     # sys.exit(1)
#
#     gene_by_name_and_chrom = build_gene_objects_list(features_bed, target.gene_keys_list)
#
#     for (bed_fpath, target_type) in zip([target.bed_fpath, exons_and_cds_features], ['amplicons', 'exons']):  # features are canonical
#         if not bed_fpath:
#             continue
#         info()
#         info('Calculating coverage statistics for ' + ('CDS and ncRNA exons...' if target_type == 'exons' else 'the regions in the target BED file...'))
#         sambamba_depth_output_fpath = sambamba_depth(work_dir, bed_fpath, sample.bam, depth_thresholds, reuse=cfg.reuse_intermediate)
#         regions = _sambamba_depth_to_regions(sambamba_depth_output_fpath, sample.name, target_type, depth_thresholds)
#
#         #####################################
#         #####################################
#         info('Filling regions stats into the gene objects...')
#         cur_unannotated_gene = None
#         total_regions_count = 0
#         for region in regions:
#             if region.feature == 'Capture':
#                 if region.gene_name != '.':
#                     cur_unannotated_gene = None
#                     gene = gene_by_name_and_chrom[(region.gene_name, region.chrom)]
#                     if (gene.gene_name, gene.chrom) not in ready_to_report_set:
#                         ready_to_report_genes.append(gene)
#                         ready_to_report_set.add((gene.gene_name, gene.chrom))
#                     gene.add_amplicon(region)
#                 else:
#                     if cur_unannotated_gene is None:
#                         cur_unannotated_gene = GeneInfo(gene_name=region.gene_name, chrom=region.chrom, feature='NotAnnotatedSummary')
#                         ready_to_report_genes.append(cur_unannotated_gene)
#                     cur_unannotated_gene.add_amplicon(region)
#
#             else:
#                 gene = gene_by_name_and_chrom[(region.gene_name, region.chrom)]
#                 gene.add_exon(region)
#                 if not target.bed_fpath:  # in case if only reporting based on features_bed
#                     if (gene.gene_name, gene.chrom) not in ready_to_report_set:
#                         ready_to_report_genes.append(gene)
#                         ready_to_report_set.add((gene.gene_name, gene.chrom))
#
#             # row = [region.chrom, region.start, region.end, region.get_size(), region.gene_name, region.strand,
#             #        region.feature, region.biotype, region.transcript_id, region.min_depth, region.avg_depth, region.std_dev,
#             #        region.rate_within_normal]
#             # row = [Metric.format_value(val, human_readable=True) for val in row]
#             # rates = [Metric.format_value(val, unit='%', human_readable=True) for val in region.rates_within_threshs.values()]
#             # row.extend(rates)
#             # col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]
#
#             total_regions_count += 1
#             if total_regions_count > 0 and total_regions_count % 10000 == 0:
#                 info('  Processed {0:,} regions'.format(total_regions_count))
#         info('Processed {0:,} regions'.format(total_regions_count))
#
#     #####################################
#     #####################################
#     report = PerRegionSampleReport(sample=sample, metric_storage=get_detailed_metric_storage(depth_thresholds))
#     # report.add_record('Sample', sample.name)
#     report.txt_fpath = sample.targqc_region_txt
#     report.tsv_fpath = sample.targqc_region_tsv
#
#     debug('Arranging regions to report...')
#     regions = []
#     for g in ready_to_report_genes:
#         for a in g.get_amplicons():
#             regions.append(a)
#         for e in g.get_exons():
#             regions.append(e)
#         if g.get_exons():
#             process_gene(g, depth_thresholds)
#             regions.append(g)
#
#     debug('Preparting report rows...')
#     for reg in regions:
#         r = report.add_row()
#         r.add_record('Chr', reg.chrom)
#         r.add_record('Start', reg.start)
#         r.add_record('End', reg.end)
#         r.add_record('Size', reg.get_size())
#         r.add_record('Gene', reg.gene_name)
#         r.add_record('Strand', reg.strand)
#         r.add_record('Feature', reg.feature)
#         r.add_record('Biotype', reg.biotype)
#         r.add_record('Transcript', reg.transcript_id)
#         r.add_record('Min depth', reg.min_depth)
#         r.add_record('Median depth', reg.median_depth)
#         r.add_record('Ave depth', reg.avg_depth)
#         r.add_record('Std dev', reg.std_dev)
#         r.add_record('W/n 20% of median depth', reg.rate_within_normal)
#         for ths in depth_thresholds:
#             r.add_record('{}x'.format(ths), reg.rates_within_threshs.get(ths) if reg.rates_within_threshs else None)
#
#     safe_mkdir(dirname(report.tsv_fpath))
#     debug('Flattening TSV records...')
#     header_rows, flat_rows = report.flatten(None, human_readable=False)
#     debug('Writing TSV...')
#     write_tsv_rows((header_rows, flat_rows), report.tsv_fpath)
#     debug('Flattening TXT records...')
#     header_rows, flat_rows = report.flatten(None, human_readable=True)
#     debug('Writing TXT...')
#     write_txt_rows((header_rows, flat_rows), report.txt_fpath)
#
#     # un_annotated_summary_region = next((g for g in gene_by_name_and_chrom.values() if g.gene_name == '.'), None)
#     # if un_annotated_summary_region and un_annotated_amplicons:
#     #     un_annotated_summary_region.feature = 'NotAnnotatedSummary'
#     #     for ampl in un_annotated_amplicons:
#     #         ampl.gene_name = un_annotated_summary_region.gene_name
#     #         un_annotated_summary_region.add_amplicon(ampl)
#     #         add_region_to_report(report, ampl, depth_thresholds)
#     #     add_region_to_report(report, un_annotated_summary_region, depth_thresholds)
#
#     # report.save_txt(sample.targqc_region_txt)
#     # report.save_tsv(sample.targqc_region_tsv)
#     info('Regions (total ' + str(len(report.rows)) + ') saved into:')
#     info('  ' + report.txt_fpath)
#     return report


# def process_gene(gene, depth_thresholds):
#     gene.rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)
#     if gene.size == 0:
#         gene.start = None
#         gene.end = None
#         return
#     exons = gene.get_exons()
#     if not exons:
#         return
#
#     total_depth = sum(e.avg_depth * e.size for e in exons)
#     gene.size = sum(e.size for e in exons)
#     gene.avg_depth = total_depth / gene.size
#     if exons[0].std_dev not in ['.', '', None]:
#         sum_of_sq_var = sum((((e.avg_depth - e.std_dev) - gene.avg_depth) ** 2 + ((e.avg_depth + e.std_dev) - gene.avg_depth) ** 2) * e.size for e in exons)
#         gene.std_dev = math.sqrt(sum_of_sq_var / 2 / float(gene.size))
#     for t in depth_thresholds:
#         total_rate = sum(e.rates_within_threshs[t] * e.size for e in exons)
#         rate = total_rate / gene.size
#         gene.rates_within_threshs[t] = rate


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


