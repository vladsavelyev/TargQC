import shutil
import os
from os import listdir
from os.path import relpath, join, exists, dirname, basename, abspath, splitext
from collections import OrderedDict, defaultdict

import targqc
import targqc.config as tc
from ngs_utils.file_utils import verify_dir, verify_file, adjust_path, symlink_plus, file_transaction, add_suffix
from ngs_utils.logger import info, err, debug
from ngs_utils.reporting.reporting import PerRegionSampleReport, BaseReport, Metric, ReportSection, MetricStorage, \
    FullReport
from targqc.general_report import get_header_metric_storage
from targqc.qualimap.runner import run_multisample_qualimap


def _make_targetcov_symlinks(samples):
    for sample in samples:
        new_link = join(
            dirname(dirname(sample.targetcov_detailed_txt)),
            basename(sample.targetcov_detailed_txt))
        if exists(new_link):
            os.unlink(new_link)
        symlink_plus(sample.targetcov_detailed_txt, new_link)
        info('TargetCov TXT symlink saved to ' + new_link)


def make_tarqc_html_report(output_dir, work_dir, samples, bed_fpath=None, tag_by_sample=None):
    # header_storage = get_header_metric_storage(tc.depth_thresholds,
    #                                            is_wgs=bed_fpath is not None,
    #                                            padding=tc.padding)

    jsons_by_sample = {s.name: s.targqc_json_fpath for s in samples if verify_file(s.targqc_json_fpath)}
    # htmls_by_sample = {s.name: s.targqc_html_fpath for s in samples if verify_file(s.targqc_html_fpath)}
    htmls_by_sample = dict()

    if not jsons_by_sample:
        return None, None, None

    targqc_full_report = FullReport.construct_from_sample_report_jsons(samples, output_dir, jsons_by_sample, htmls_by_sample)

    for sample_report in targqc_full_report.sample_reports:
        if tag_by_sample:
            sample_report.set_project_tag(tag_by_sample[sample_report.sample.name])
        if verify_file(sample_report.sample.qualimap_html_fpath):
            url = relpath(sample_report.sample.qualimap_html_fpath, output_dir)
            r = sample_report.find_record(sample_report.records, 'Qualimap')
            if r:
                r.url = url
            else:
                sample_report.add_record(metric_name='Qualimap', value='Qualimap', url=url, silent=True)

    if len(samples) > 1:
        run_multisample_qualimap(output_dir, work_dir, samples, targqc_full_report)

    fn = splitext(basename(samples[0].targqc_txt_fpath))[0]
    tsv_fpath = targqc_full_report.save_tsv(join(output_dir, fn + '.tsv'))
    html_fpath = targqc_full_report.save_html(join(output_dir, fn + '.html'), 'TargQC')

    return tsv_fpath, html_fpath


def combined_regional_reports(work_dir, output_dir, samples):
    if not any(verify_file(s.targqc_region_tsv, silent=True) for s in samples):
        return None, None

    tsv_region_rep_fpath = join(output_dir, basename(samples[0].targqc_region_tsv))
    debug('Combining regional reports, writing to ' + tsv_region_rep_fpath)
    with file_transaction(work_dir, tsv_region_rep_fpath) as tx_tsv:
        with open(tx_tsv, 'w') as tsv_out:
            # sample_i = 0
            # for s in samples:
            #     if s.targqc_region_txt and verify_file(s.targqc_region_txt):
            #         with open(s.targqc_region_txt) as txt_in:
            #             for l in txt_in:
            #                 if l.startswith('#'):
            #                     if not l.startswith('##') and sample_i == 0:
            #                         txt_out.write('#Sample' + ' '*(max(len('#Sample'), len(s.name)) - len('#Sample')) + ' ' + l.replace('#Chr', 'Chr '))
            #                 else:
            #                     txt_out.write(s.name + ' '*(max(len('#Sample'), len(s.name)) - len(s.name)) + ' ' + l)
            #         sample_i += 1
            sample_i = 0
            for s in samples:
                if s.targqc_region_tsv and verify_file(s.targqc_region_tsv):
                    with open(s.targqc_region_tsv) as tsv_in:
                        for l in tsv_in:
                            if l.startswith('#'):
                                if not l.startswith('##') and sample_i == 0:
                                    tsv_out.write('#Sample\t' + l[1:])
                            else:
                                tsv_out.write(s.name + '\t' + l)
                    sample_i += 1

    return tsv_region_rep_fpath


# def summarize_targqc(summary_threads, output_dir, work_dir, samples, bed_fpath=None, tag_by_sample=None):
    # best_for_regions_fpath = None
    # if any(verify_file(s.targqc_region_tsv, silent=True) for s in samples):
    #     best_for_regions_fpath = _save_best_details_for_each_gene(tc.depth_thresholds, samples, output_dir)
    # ''' 1. best_regions = get_best_regions()
    #     2. best_for_regions_fpath = save_per_region_report()
    #     3. calc median coverage across best regions
    #     4. flagged_regions_report_fpath = _generate_flagged_regions_report(
    #          output_dir, 'Best', average_coverage, genes, depth_threshs)
    # '''


def _prep_best_report(metric_storage, samples):
    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)

    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    m = metric_storage.find_metric('Average sample depth')
    ave_sample_depth = max(BaseReport.find_record(s.report.records, m.name).value for s in samples)
    report.add_record('Average sample depth', ave_sample_depth)

    return report


def _get_targqc_metric_storage(metric_storages_by_report_type):
    class SectionId:
        def __init__(self, name, title):
            self.name = name
            self.title = title

        def __hash__(self):
            #return hash((self.name, self.title))
            return hash(self.name)  # use title from the first metric_storage

        def __eq__(self, other):
            #return (self.name, self.title) == (other.name, other.title)
            return self.name == other.name  # use title from the first metric_storage

    metrics_by_sections = OrderedDict()
    general_section_id = None
    general_section_metric_list = []

    for report_type, metric_storage in metric_storages_by_report_type:
        for section in metric_storage.sections:
            section_id = SectionId(section.name, section.title)
            if section_id not in metrics_by_sections.keys():
                metrics_by_sections[section_id] = []

            metrics_by_sections[section_id] += [metric
                for metric in metric_storage.get_metrics(sections=[section], skip_general_section=True)
                if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]

        # specific behaviour for general section
        general_section_metric_list += [metric
            for metric in metric_storage.general_section.metrics
            if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]
        if not general_section_id:
            general_section_id = SectionId(metric_storage.general_section.name, metric_storage.general_section.title)

    sections = []
    for section_id, metric_list in metrics_by_sections.items():
        sections.append(ReportSection(section_id.name, section_id.title, metric_list))

    return MetricStorage(
        general_section=ReportSection(
            general_section_id.name, general_section_id.title, general_section_metric_list),
        sections=sections)


def _correct_qualimap_genome_results(samples):
    """ fixing java.lang.Double.parseDouble error on entries like "6,082.49"
    """
    for s in samples:
        if verify_file(s.qualimap_genome_results_fpath):
            correction_is_needed = False
            with open(s.qualimap_genome_results_fpath, 'r') as f:
                content = f.readlines()
                metrics_started = False
                for line in content:
                    if ">> Reference" in line:
                        metrics_started = True
                    if metrics_started:
                        if line.find(',') != -1:
                            correction_is_needed = True
                            break
            if correction_is_needed:
                with open(s.qualimap_genome_results_fpath, 'w') as f:
                    metrics_started = False
                    for line in content:
                        if ">> Reference" in line:
                            metrics_started = True
                        if metrics_started:
                            if line.find(',') != -1:
                                line = line.replace(',', '')
                        f.write(line)


def _correct_qualimap_insert_size_histogram(work_dir, samples):
    """ replacing Qualimap insert size histogram with Picard one.
    """
    for s in samples:
        qualimap1_dirname = dirname(s.qualimap_ins_size_hist_fpath).replace('raw_data_qualimapReport', 'raw_data')
        qualimap2_dirname = dirname(s.qualimap_ins_size_hist_fpath)
        if exists(qualimap1_dirname):
            if not exists(qualimap2_dirname):
                shutil.move(qualimap1_dirname, qualimap2_dirname)
            else:
                shutil.rmtree(qualimap1_dirname)
        elif not exists(qualimap2_dirname):
            continue  # no data from both Qualimap v.1 and Qualimap v.2

        # if qualimap histogram exits and reuse_intermediate, skip
        if verify_file(s.qualimap_ins_size_hist_fpath, silent=True) and tc.reuse_intermediate:
            pass
        else:
            if verify_file(s.picard_ins_size_hist_txt_fpath):
                with open(s.picard_ins_size_hist_txt_fpath, 'r') as picard_f:
                    one_line_to_stop = False
                    for line in picard_f:
                        if one_line_to_stop:
                            break
                        if line.startswith('## HISTOGRAM'):
                            one_line_to_stop = True

                    with file_transaction(work_dir, s.qualimap_ins_size_hist_fpath) as tx:
                        with open(tx, 'w') as qualimap_f:
                            for line in picard_f:
                                qualimap_f.write(line)


def select_best(values, fn=max):
    vs = [v for v in values if v is not None]
    return fn(vs) if len(vs) > 0 else None


def get_int_val(v):
    v = _get_num(v)
    return int(v) if v else None

def get_float_val(v):
    v = _get_num(v)
    return float(v) if v else None

def _get_num(v):
    v = get_val(v)
    return ''.join(c for c in v if c.isdigit() or c == '.') if v else None

def get_val(v):
    return v.strip() if v.strip() not in ['.', '-', ''] else None


def _save_best_details_for_each_gene(depth_threshs, samples, output_dir):
    metric_storage = get_detailed_metric_storage(depth_threshs)

    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)
    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    total_regions = 0
    fpaths = [s.targqc_region_tsv for s in samples if verify_file(s.targqc_region_tsv)]
    if not fpaths:
        err('No targetcov detailed per-gene report was generated; skipping.')
        return None

    open_tsv_files = [open(fpath) for fpath in fpaths]

    first_col = 0
    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break
        l = lines_for_each_sample[0]
        if l.startswith('##'):
            continue
        elif l.startswith('#'):
            if l.startswith('#Sample'):
                first_col = 1
            break

    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break

        if all([not l.startswith('#') and ('Whole-Gene' in l or 'Gene-Exon' in l) for l in lines_for_each_sample]):
            shared_fields = lines_for_each_sample[0].split('\t')[first_col:first_col+9]
            reg = report.add_row()
            reg.add_record('Chr', get_val(shared_fields[0]))
            reg.add_record('Start', get_int_val(shared_fields[1]))
            reg.add_record('End', get_int_val(shared_fields[2]))
            reg.add_record('Size', get_int_val(shared_fields[3]))
            reg.add_record('Gene', get_val(shared_fields[4]))
            reg.add_record('Strand', get_val(shared_fields[5]))
            reg.add_record('Feature', get_val(shared_fields[6]))
            reg.add_record('Biotype', get_val(shared_fields[7]))
            reg.add_record('Transcript', get_val(shared_fields[8]))

            min_depths, ave_depths, stddevs, withins = ([], [], [], [])
            percents_by_threshs = {t: [] for t in depth_threshs}

            for l in lines_for_each_sample:
                fs = l.split('\t')

                min_depths.append(get_int_val(fs[first_col+9]))
                ave_depths.append(get_float_val(fs[first_col+10]))
                stddevs.append(get_float_val(fs[first_col+11]))
                withins.append(get_float_val(fs[first_col+12]))
                for t, f in zip(depth_threshs, fs[first_col+13:]):
                    percents_by_threshs[t].append(get_float_val(f))

            # counting bests
            reg.add_record('Min depth', select_best(min_depths))
            reg.add_record('Ave depth', select_best(ave_depths))
            reg.add_record('Std dev', select_best(stddevs, max))
            reg.add_record('W/n 20% of median depth', select_best(withins))
            for t in depth_threshs:
                reg.add_record('{}x'.format(t), select_best(percents_by_threshs[t]))

            total_regions += 1

    for f in open_tsv_files:
        f.close()

    gene_report_basename = add_suffix(samples[0].targqc_region_tsv, 'best')
    txt_rep_fpath = report.save_txt(join(output_dir, gene_report_basename + '.txt'))
    tsv_rep_fpath = report.save_tsv(join(output_dir, gene_report_basename + '.tsv'))
    info('')
    info('Best values for the regions (total ' + str(total_regions) + ') saved into:')
    info('  ' + txt_rep_fpath)

    return txt_rep_fpath




