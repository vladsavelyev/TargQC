from os.path import getsize, dirname, join, abspath

from Utils.call_process import run
from Utils.file_utils import safe_mkdir, verify_file
from Utils.logger import info

import targqc.config as tc


def get_qualimap_max_mem(bam):
    mem_m = getsize(bam) / 3 / 1024 / 1024
    mem_m = min(max(mem_m, 1200), 16000)
    return mem_m


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


def run_qualimap(output_dir, bam_fpath, bed_fpath=None, threads=1):
    info('Running QualiMap')
    safe_mkdir(dirname(output_dir))
    safe_mkdir(output_dir)

    bed_cmd = ''
    if bed_fpath:
        qualimap_bed_fpath = join(output_dir, 'tmp_qualimap.bed')
        fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)
        bed_cmd = ' -gff ' + qualimap_bed_fpath + ' '
        info('Using amplicons/capture panel ' + qualimap_bed_fpath)

    qualimap = abspath(join(dirname(__file__), 'qualimap', 'qualimap'))

    mem_cmdl = ''
    mem_m = get_qualimap_max_mem(bam_fpath)
    mem = str(int(mem_m)) + 'M'
    mem_cmdl = '--java-mem-size=' + mem

    cmdline = ('{qualimap} bamqc --skip-duplicated -nt {threads} {mem_cmdl} -nr 5000 '
        '-bam {bam_fpath} -outdir {output_dir} {bed_cmd} -c -gd HUMAN').format(**locals())
    report_fpath = join(output_dir, 'qualimapReport.html')

    run(cmdline, output_fpath=report_fpath, stdout_to_outputfile=False, env_vars=dict(DISPLAY=None),
        checks=[lambda _1, _2: verify_file(report_fpath)], reuse=tc.reuse_intermediate)
