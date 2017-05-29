import os
import sys
import shutil
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getsize, getmtime
from datetime import datetime
from collections import namedtuple

from ngs_utils.testing import BaseTestCase, info, check_call, swap_output


class BaseTargQC(BaseTestCase):
    script = 'targqc'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir)
    gold_standard_dir = join(join(dirname(__file__), BaseTestCase.gold_standard_dir))

    syn3_url = 'http://quast.bioinf.spbau.ru/static/chr21.tar.gz'
    bwa_url = 'http://quast.bioinf.spbau.ru/static/bwa.tar.gz'

    Sample = namedtuple('Sample', 'name bam l_fastq r_fastq')
    samples = [
        Sample('syn3-tumor', 'syn3-tumor.bam', 'syn3-tumor_R1.fq.gz', 'syn3-tumor_R2.fq.gz'),
        Sample('syn3-normal', 'syn3-normal.bam', 'syn3-normal_R1.fq.gz', 'syn3-normal_R2.fq.gz'),
    ]

    bwa_dir = join(data_dir, 'bwa')
    bwa_path = join(bwa_dir, 'hg19-chr21.fa')

    def setUp(self):
        BaseTestCase.setUp(self)

        self.syn3_dir = join(BaseTargQC.data_dir, 'chr21')
        self.bed3 = join(self.syn3_dir, 'NGv3.chr21.3col.bed')
        self.bed4 = join(self.syn3_dir, 'NGv3.chr21.4col.bed')
        self.bams = [join(self.syn3_dir, s.bam) for s in BaseTargQC.samples]
        self.fastqs = [join(self.syn3_dir, s.l_fastq) for s in BaseTargQC.samples] + \
                      [join(self.syn3_dir, s.r_fastq) for s in BaseTargQC.samples]

        if not isdir(self.syn3_dir):
            info(self.syn3_dir + ' does not exist, downloading test data')
            cur_dir = os.getcwd()
            os.chdir(self.data_dir)
            check_call(['wget', self.syn3_url])
            check_call(['tar', '-xzvf', basename(self.syn3_url)])
            os.chdir(cur_dir)
        if not isdir(self.bwa_dir):
            info(self.bwa_dir + ' does not exist, downloading bwa reference data')
            cur_dir = os.getcwd()
            os.chdir(self.data_dir)
            check_call(['wget', self.bwa_url])
            check_call(['tar', '-xzvf', basename(self.bwa_url)])
            os.chdir(cur_dir)

    def _test(self, output_dirname=None, used_samples=samples, bams=None, fastq=None, bed=None,
              debug=False, reuse_intermediate=False, reuse_output_dir=False, reannotate=False,
              genome='hg19-chr21', bwa=None, threads=None, ipython=None, keep_work_dir=True):
        os.chdir(self.results_dir)
        cmdl = [self.script]
        output_dir = None
        if output_dirname:
            output_dir = join(self.results_dir, output_dirname)
            cmdl.extend(['-o', output_dir])
        if bams: cmdl.extend(bams)
        if fastq: cmdl.extend(fastq)
        if bed: cmdl.extend(['--bed', bed])
        if debug: cmdl.append('--debug')
        if reuse_intermediate:
            reuse_output_dir = True
            cmdl.append('--reuse')
        if reannotate: cmdl.append('--reannotate')
        if genome: 
            cmdl.extend(['-g', genome])
        if bwa:
            cmdl.extend(['--bwa', bwa])
        if threads: cmdl.extend(['-t', str(threads)])
        if ipython: cmdl.extend('-s sge -q queue -r pename=smp -r --local'.split())
        if keep_work_dir: cmdl.append('--keep-work-dir')

        output_dir = output_dir or self._default_output_dir()
        if reuse_output_dir is False:
            swap_output(output_dir)

        info('-' * 100)
        check_call(cmdl)
        info('-' * 100)
        info('')

        self._check_results(output_dir, used_samples)

        if self.remove_work_dir_on_success and not reuse_intermediate and not reuse_output_dir:
            work_dir = join(output_dir, 'work')
            if not isdir(work_dir):
                info('Work dir for run ' + output_dirname + ' does not exist under ' + work_dir)
            else:
                shutil.rmtree(work_dir)
        info('')

    @staticmethod
    def _default_output_dir():
        return join(os.getcwd(), 'targqc')

    def _check_results(self, output_dir, used_samples):
        assert isdir(output_dir)
        self._check_file_throws(join(output_dir, 'regions.tsv'))
        self._check_file_throws(join(output_dir, 'summary.tsv'))
        self._check_file_throws(join(output_dir, 'summary.html'), ignore_matching_lines='report_date')
        for s in used_samples:
            s_dir = join(output_dir, s.name)
            assert isdir(s_dir)
            self._check_file_throws(join(s_dir, 'regions.tsv'))
            self._check_file_throws(join(s_dir, 'summary.txt'))
            self._check_file_throws(join(s_dir, 'summary.html'), ignore_matching_lines='report_date')
            self._check_file_throws(join(s_dir, 'summary.json'), ignore_matching_lines='work_dir')
        # TODO: check line numbers and some values isntead of diff?
