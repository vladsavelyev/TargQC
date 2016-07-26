from collections import defaultdict
from os.path import isfile, join
from pybedtools import BedTool
import GeneAnnotation as ga
from GeneAnnotation.annotate_bed import annotate, get_sort_key, tx_sort_key
from Utils import reference_data
from Utils.bed_utils import filter_bed_with_gene_set, get_gene_keys_from_bed, calc_region_number, merge_overlaps
from Utils.bed_utils import remove_comments, sort_bed, count_bed_cols, cut, verify_bed
from Utils.file_utils import iterate_file, add_suffix, intermediate_fname, file_transaction, verify_file
from Utils.logger import debug, info, warn
from Utils.utils import OrderedDefaultDict
from targqc import config as cfg


class Target:
    def __init__(self, work_dir, fai_fpath, reuse, bed_fpath=None):
        self.bed = None
        self.original_bed_fpath = None
        self.bed_fpath = None
        self.qualimap_bed_fpath = None
        self.padded_bed_fpath = None

        self.gene_keys_set = set()  # set of pairs (gene_name, chrom)
        self.gene_keys_list = list()  # list of pairs (gene_name, chrom)
        self.regions_num = None

        self.bases_num = None
        self.fraction = None

        if bed_fpath:
            self.is_wgs = False
            verify_bed(bed_fpath, is_critical=True)
            self.original_bed_fpath = bed_fpath
            self._make_target_bed(bed_fpath, work_dir, fai_fpath, reuse)
        else:
            self.is_wgs = True
            self._make_wgs_regions_file(work_dir)

    def get_capture_bed(self):
        if not self.is_wgs:
            return self.bed.filter(lambda x: x[ga.BedCols.FEATURE] == 'capture')
        else:
            return None

    def _make_target_bed(self, bed_fpath, work_dir, fai_fpath=None, reuse=False):
        debug()
        info('Cleaning target BED file...')
        clean_target_bed_fpath = intermediate_fname(work_dir, bed_fpath, 'clean')
        if isfile(clean_target_bed_fpath) and verify_file(clean_target_bed_fpath, cmp_date_fpath=bed_fpath):
            pass
        else:
            bed = BedTool(bed_fpath)\
                .filter(lambda x: x.chrom and not any(x.chrom.startswith(e) for e in ['#', ' ', 'track', 'browser']))\
                .remove_invalid()\
                .cut(range(3))
            with file_transaction(work_dir, clean_target_bed_fpath) as tx:
                bed.saveas(tx)
            debug('Saved to ' + clean_target_bed_fpath)

        debug()
        info('Sorting target BED file...')
        sort_target_bed_fpath = sort_bed(clean_target_bed_fpath, work_dir=work_dir, fai_fpath=fai_fpath, reuse=reuse)
        debug('Saved to ' + sort_target_bed_fpath)

        debug()
        info('Annotating target BED file and collecting overlapping genome features...')
        ann_target_bed_fpath = add_suffix(sort_target_bed_fpath, 'ann')
        annotate(sort_target_bed_fpath, ann_target_bed_fpath, reuse=reuse, genome=cfg.genome,
                 work_dir=work_dir, is_debug=cfg.debug, extended=True, output_features=True)
        # TODO prepare BEDs: annotate with --report-features, then use that whole BED to find coverage, but only "cature" for summary reports
        debug('Saved to ' + ann_target_bed_fpath)

        self.bed_fpath = intermediate_fname(work_dir, ann_target_bed_fpath, 'clean')
        if isfile(self.bed_fpath) and verify_file(self.bed_fpath, cmp_date_fpath=ann_target_bed_fpath):
            pass
        else:
            self.bed = BedTool(ann_target_bed_fpath).remove_invalid()
            with file_transaction(work_dir, self.bed_fpath) as tx:
                self.bed.saveas(tx)
        self.bed = BedTool(self.bed_fpath)

        gene_key_set, gene_key_list = get_gene_keys_from_bed(bed_fpath)
        self.gene_keys_set = gene_key_set
        self.gene_keys_list = gene_key_list
        self.regions_num = self.get_capture_bed().count()

        self._make_qualimap_bed(work_dir)
        self._make_padded_bed(work_dir)

    def _make_padded_bed(self, work_dir):
        if self.is_wgs:
            return None

        self.padded_bed_fpath = intermediate_fname(work_dir, self.bed_fpath, 'padded')
        if isfile(self.padded_bed_fpath) and verify_file(self.padded_bed_fpath, cmp_date_fpath=self.bed_fpath):
            debug('Padded BED file ' + self.padded_bed_fpath + ' is ready, reusing')
            return BedTool(self.padded_bed_fpath)

        padded_bed = self.bed.slop(b=cfg.padding, g=cfg.fai_fpath).sort().merge()
        with file_transaction(work_dir, self.padded_bed_fpath) as tx:
            padded_bed.saveas(tx)
        return BedTool(self.padded_bed_fpath)

    def _make_qualimap_bed(self, work_dir):
        if self.is_wgs:
            return None

        self.qualimap_bed_fpath = intermediate_fname(work_dir, self.bed_fpath, 'qualimap_ready')
        if isfile(self.qualimap_bed_fpath) and verify_file(self.qualimap_bed_fpath, cmp_date_fpath=self.bed_fpath):
            debug('Qualimap-ready BED file ' + self.qualimap_bed_fpath + ' is ready, reusing')
            return self.qualimap_bed_fpath

        debug('Merging and saving BED into required bed6 format for Qualimap')
        bed = self.bed.sort().merge()
        with file_transaction(work_dir, self.qualimap_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for i, region in enumerate(x for x in bed):
                    region = [x for x in list(region) if x]
                    fillers = [str(i), "1.0", "+"]
                    full = region + fillers[:6 - len(region)]
                    out.write("\t".join(full) + "\n")
        return self.qualimap_bed_fpath

    def _make_wgs_regions_file(self, work_dir):
        self.wgs_bed_fpath = join(work_dir, 'targqc_features_to_report.bed')
        if isfile(self.wgs_bed_fpath) and cfg.reuse_intermediate and verify_file(self.wgs_bed_fpath):
            debug(self.wgs_bed_fpath + ' exists, reusing')
            return self.wgs_bed_fpath

        chr_order = reference_data.get_chrom_order(cfg.genome)

        r_by_tx_by_gene = OrderedDefaultDict(lambda: defaultdict(list))
        all_features = ga.get_all_features(cfg.genome, high_confidence=True)
        all_features = all_features.saveas('/Users/vlad/vagrant/Pre_Alignment_Suite/tmp.bed')

        info('Select best transcript to report')
        for r in all_features:
            if r[ga.BedCols.FEATURE] != 'gene':
                gene = r[ga.BedCols.HUGO]
                tx = r[ga.BedCols.ENSEMBL_ID]
                r_by_tx_by_gene[gene][tx].append(r)

        with file_transaction(work_dir, self.wgs_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for gname, r_by_tx in r_by_tx_by_gene.iteritems():
                    all_tx = (x for xx in r_by_tx.itervalues() for x in xx if x[ga.BedCols.FEATURE] == 'transcript')
                    tx_sorted_list = [x[ga.BedCols.ENSEMBL_ID] for x in sorted(all_tx, key=tx_sort_key)]
                    if not tx_sorted_list:
                        continue
                    tx_id = tx_sorted_list[0]
                    for r in sorted(r_by_tx[tx_id], key=get_sort_key(chr_order)):
                        out.write('\t'.join(str(f) for f in r) + '\n')
        return self.wgs_bed_fpath