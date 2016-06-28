import subprocess
from os.path import dirname, basename, join, abspath, isdir

import shutil

from Utils.call_process import run
from Utils.file_utils import safe_mkdir, verify_file, file_transaction
from Utils.logger import info, err


picard = 'java -jar ' + abspath(join(dirname(__file__), 'picard', 'picard.jar'))


def picard_ins_size_hist(sample, bam_fpath):
    safe_mkdir(dirname(sample.picard_ins_size_hist_txt_fpath))
    safe_mkdir(dirname(sample.picard_ins_size_hist_pdf_fpath))
    info('Picard ins size hist for "' + basename(bam_fpath) + '"')
    cmdline = picard + ' CollectInsertSizeMetrics' \
              ' I={bam_fpath}' \
              ' O={sample.picard_ins_size_hist_txt_fpath}' \
              ' H={sample.picard_ins_size_hist_pdf_fpath}' \
              ' VALIDATION_STRINGENCY=LENIENT'
    cmdline = cmdline.format(**locals())
    run(cmdline, output_fpath=sample.picard_ins_size_hist_txt_fpath, stdout_to_outputfile=False)


def replace_qualimap_insert_size_histogram(cnf, sample):
    """ replacing Qualimap insert size histogram with Picard one.
    """
    qualimap1_dirname = dirname(sample.qualimap_ins_size_hist_fpath).replace('raw_data_qualimapReport', 'raw_data')
    qualimap2_dirname = dirname(sample.qualimap_ins_size_hist_fpath)
    if isdir(qualimap1_dirname):
        if not dir(qualimap2_dirname):
            shutil.move(qualimap1_dirname, qualimap2_dirname)
        else:
            shutil.rmtree(qualimap1_dirname)
    elif not isdir(qualimap2_dirname):
        return  # no data from both Qualimap v.1 and Qualimap v.2

    # if qualimap histogram exits and reuse_intermediate, skip
    if verify_file(sample.qualimap_ins_size_hist_fpath, silent=True) and cnf.reuse_intermediate:
        pass
    else:
        if verify_file(sample.picard_ins_size_hist_txt_fpath):
            with open(sample.picard_ins_size_hist_txt_fpath, 'r') as picard_f:
                one_line_to_stop = False
                for line in picard_f:
                    if one_line_to_stop:
                        break
                    if line.startswith('## HISTOGRAM'):
                        one_line_to_stop = True

                with file_transaction(cnf.work_dir, sample.qualimap_ins_size_hist_fpath) as tx:
                    with open(tx, 'w') as qualimap_f:
                        for line in picard_f:
                            qualimap_f.write(line)
