from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath

from nose import SkipTest

from . import BaseTargQC, info


class UnitTests(BaseTargQC):
    def test_01_simple(self):
        self._test('simple', [self.samples[0]], bams=[self.bams[0]], bed=self.bed4)

    def test_02_bed3(self):
        self._test('bed3', bams=self.bams, bed=self.bed3)

    def test_03_bed4(self):
        self._test('bed4', bams=self.bams, bed=self.bed4)

    def test_04_bed4_reannotate(self):
        self._test('bed4_reannotate', bams=self.bams, bed=self.bed4, reannotate=True)

    def test_05_wgs(self):
        self._test('wgs', bams=self.bams)

    def test_06_threads(self):
        self._test('threads', bams=self.bams, bed=self.bed4, threads='2')

    def test_07_ipython(self):
        self._test('ipython', bams=self.bams, bed=self.bed4, ipython=True, threads='2')

    def test_08_fastq(self):
        self._test('fastq', fastq=self.fastqs, bwa=self.bwa_path, bed=self.bed4)

    def test_09_debug_and_reuse(self):
        self._test('debug_and_reuse', bams=self.bams, bed=self.bed4, debug=True, reuse_intermediate=True)

    def test_10_reuse_output(self):
        self._test('reuse_output_dir', bams=self.bams, bed=self.bed4, reuse_output_dir=True)

    def test_11_full_hg19(self):
        raise SkipTest

    def test_12_full_hg38(self):
        raise SkipTest

    # def test_13_api(self):
    #     import targqc
    #     import ngs_utils.reference_data as ref
    #     from ngs_utils.file_utils import safe_mkdir
    #     from ngs_utils.parallel import ParallelCfg
    #
    #     genome = 'hg19-chr21'
    #     fai_fpath = ref.get_fai(genome)
    #     output_dirname = 'api'
    #     output_dir = join(self.results_dir, output_dirname)
    #     work_dir = join(output_dir, 'work')
    #     samples = sorted([targqc.Sample(s.name,
    #          dirpath=safe_mkdir(join(output_dir, s.name)),
    #          work_dir=safe_mkdir(join(work_dir, s.name)),
    #          bam=join(self.syn3_dir, s.bam))
    #                for s in self.samples], key=lambda _s: _s.key_to_sort())
    #     parallel_cfg = ParallelCfg(None, None, None, 1, None)
    #     info('-' * 100)
    #     targqc.start_targqc(work_dir, output_dir, samples, self.bed4,
    #                         parallel_cfg, self.bwa_path,
    #                         fai_fpath=fai_fpath,
    #                         genome=genome)
    #     info('-' * 100)
    #     info('')
    #
    #     info()
    #     self._check_results(output_dir, self.samples)
