import os
from genericpath import getmtime
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath

from datetime import datetime
from collections import namedtuple
from nose import SkipTest

from ngs_utils.file_utils import add_suffix
from ngs_utils.testing import BaseTestCase, swap_output, check_call
from tests import BaseTargQC, info


class AnnotateBedTests(BaseTestCase):
    script = 'annotate_bed.py'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir, 'bed')
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, 'annotate_bed')
    gold_standard_dir = join(join(dirname(__file__), BaseTestCase.gold_standard_dir), 'annotate_bed')

    def test_hg19(self):
        self._test('hg19', 'hg19')

    def test_hg38(self):
        self._test('hg38', 'hg38')

    def test_hg38_seq2c(self):
        self._test('hg38_seq2c', 'hg38', ['--seq2c'])

    def test_hg38_ext(self):
        self._test('hg38_ext', 'hg38', ['--extended'])

    def test_hg38_ext_features(self):
        self._test('hg38_ext_features', 'hg38', ['--extended', '--output-features'])

    def test_hg38_ext_features_highconf(self):
        self._test('hg38_ext_features_high_conf', 'hg38', ['--extended', '--output-features', '--high-confidence'])

    def test_hg38_short(self):
        self._test('hg38_short', 'hg38', ['--short'])

    def test_mm10(self):
        self._test('mm10', 'mm10', ['--extended'])

    def _test(self, name, genome, opts=None):
        os.chdir(self.results_dir)
        input_fname = genome + '.bed'
        output_fname = add_suffix(input_fname, 'anno')
        input_fpath = join(self.data_dir, input_fname)
        output_fpath = join(self.results_dir, output_fname)

        cmdl = [self.script,
                input_fpath,
                '-o', output_fpath] + \
               (opts or []) + \
               ['--debug'] + \
               ['-g', genome]

        swap_output(output_fpath)

        info('-' * 100)
        check_call(cmdl)
        info('-' * 100)
        info('')

        self._check_file(output_fpath)

