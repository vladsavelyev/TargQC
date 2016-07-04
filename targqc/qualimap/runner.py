import subprocess
from distutils.version import LooseVersion
from os.path import getsize, dirname, join, abspath, relpath, isdir, exists
from os import listdir
import shutil

from Utils.call_process import run
from Utils.file_utils import safe_mkdir, verify_file, verify_dir, file_transaction, file_exists
from Utils.logger import info, warn, err, critical
from Utils.reporting.reporting import write_tsv_rows

import targqc.config as tc


def get_qualimap_max_mem(bam):
    mem_m = getsize(bam) / 3 / 1024 / 1024
    mem_m = min(max(mem_m, 1200), 16000)
    return mem_m


def find_executable():
    root_dirpath = abspath(join(dirname(__file__)))
    for dname in listdir(abspath(join(dirname(__file__)))):
        if isdir(join(root_dirpath, dname)) and dname.startswith('qualimap') and file_exists(join(root_dirpath, dname, 'qualimap')):
            return join(root_dirpath, dname, 'qualimap')
    critical('Error: could not find Qualimap executable')


def fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath):
    with open(bed_fpath) as inn, open(qualimap_bed_fpath, 'w') as out:
        for l in inn:
            fields = l.strip().split('\t')

            if len(fields) < 3:
                continue
            try:
                int(fields[1]), int(fields[2])
            except ValueError:
                continue

            if len(fields) < 4:
                fields.append('.')

            if len(fields) < 5:
                fields.append('0')

            if len(fields) < 6:
                fields.append('+')

            out.write('\t'.join(fields) + '\n')


def run_qualimap(output_dir, bam_fpath, bed_fpath=None, threads=1, reuse_intermediate=False):
    safe_mkdir(dirname(output_dir))
    safe_mkdir(output_dir)

    bed_cmd = ''
    if bed_fpath:
        qualimap_bed_fpath = join(output_dir, 'tmp_qualimap.bed')
        fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)
        bed_cmd = ' -gff ' + qualimap_bed_fpath + ' '
        info('Using amplicons/capture panel ' + qualimap_bed_fpath)

    mem_cmdl = ''
    mem_m = get_qualimap_max_mem(bam_fpath)
    mem = str(int(mem_m)) + 'M'
    mem_cmdl = '--java-mem-size=' + mem

    cmdline = (find_executable() + ' bamqc --skip-duplicated -nt {threads} {mem_cmdl} -nr 5000 '
        '-bam {bam_fpath} -outdir {output_dir} {bed_cmd} -c -gd HUMAN').format(**locals())
    report_fpath = join(output_dir, 'qualimapReport.html')

    run(cmdline, output_fpath=report_fpath, stdout_to_outputfile=False, env_vars=dict(DISPLAY=None),
        checks=[lambda _1, _2: verify_file(report_fpath)], reuse=reuse_intermediate)


def run_multisample_qualimap(output_dir, work_dir, samples, targqc_full_report):
    """ 1. Generates Qualimap2 plots and put into plots_dirpath
        2. Adds records to targqc_full_report.plots
    """
    plots_dirpath = join(output_dir, 'plots')
    if tc.reuse_intermediate and isdir(plots_dirpath) and [f for f in listdir(plots_dirpath) if not f.startswith('.')]:
        info('Qualimap miltisample plots exist - ' + plots_dirpath + ', reusing...')
    else:
        # Qualimap2 run for multi-sample plots
        if len([s.qualimap_html_fpath for s in samples if s.qualimap_html_fpath]) > 0:
            if find_executable() is not None and get_qualimap_type(find_executable()) == 'full':
                qualimap_output_dir = join(work_dir, 'qualimap_multi_bamqc')

                _correct_qualimap_genome_results(samples)
                _correct_qualimap_insert_size_histogram(samples)

                safe_mkdir(qualimap_output_dir)
                rows = []
                for sample in samples:
                    if sample.qualimap_html_fpath:
                        rows += [[sample.name, sample.qualimap_html_fpath]]

                data_fpath = write_tsv_rows(([], rows), join(qualimap_output_dir, 'qualimap_results_by_sample.tsv'))
                qualimap_plots_dirpath = join(qualimap_output_dir, 'images_multisampleBamQcReport')
                cmdline = find_executable() + ' multi-bamqc --data {data_fpath} -outdir {qualimap_output_dir}'.format(**locals())
                run(cmdline, env_vars=dict(DISPLAY=None),
                    checks=[lambda _1, _2: verify_dir(qualimap_output_dir)], reuse=tc.reuse_intermediate)

                if not verify_dir(qualimap_plots_dirpath):
                    warn('Warning: Qualimap for multi-sample analysis failed to finish. TargQC will not contain plots.')
                    return None
                else:
                    if exists(plots_dirpath):
                        shutil.rmtree(plots_dirpath)
                    shutil.move(qualimap_plots_dirpath, plots_dirpath)
            else:
                warn('Warning: Qualimap for multi-sample analysis was not found. TargQC will not contain plots.')
                return None

    targqc_full_report.plots = []
    for plot_fpath in listdir(plots_dirpath):
        plot_fpath = join(plots_dirpath, plot_fpath)
        if verify_file(plot_fpath) and plot_fpath.endswith('.png'):
            targqc_full_report.plots.append(relpath(plot_fpath, output_dir))


def get_qualimap_type(tool_cmdline):
    """Qualimap supports multi-bamqc functionality only starting from v.2.0
    """
    if LooseVersion(_get_qualimap_version(tool_cmdline)) >= LooseVersion("2.0"):
        return "full"
    else:
        return "limited"


def _get_qualimap_version(tool_cmdline):
    cmdline = tool_cmdline + ' -version'  # actually, Qualimap doesn't have -version option

    version = None
    with subprocess.Popen(cmdline,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True).stdout as stdout:
        out = stdout.read().strip()
        flag = "QualiMap v."
        if out.startswith(flag) >= 0:
            version = out.split(flag)[-1].strip()
    if not version:
        info('WARNING: could not determine Qualimap version, using 1.0')
        return '1.0'
    if version.split('.') > 2:  # only major version
        version = '.'.join(version.split('.')[:2])
    return version


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


def _correct_qualimap_insert_size_histogram(samples):
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

                    with file_transaction(None, s.qualimap_ins_size_hist_fpath) as tx:
                        with open(tx, 'w') as qualimap_f:
                            for line in picard_f:
                                qualimap_f.write(line)
