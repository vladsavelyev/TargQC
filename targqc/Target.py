from collections import defaultdict

from ngs_utils import reference_data
from ngs_utils.bed_utils import sort_bed, verify_bed, get_genes_from_bed
from ngs_utils.file_utils import add_suffix, intermediate_fname, file_transaction, verify_file, can_reuse
from ngs_utils.logger import debug
from ngs_utils.utils import OrderedDefaultDict
from os.path import join, basename
from pybedtools import BedTool

import ensembl as ebl
from ensembl.bed_annotation import overlap_with_features, get_sort_key, tx_priority_sort_key
from targqc import config as cfg


class Target:
    def __init__(self, work_dir, output_dir, fai_fpath, bed_fpath=None,
                 padding=None, reannotate=False, genome=None, is_debug=False):
        self.bed = None
        self.original_bed_fpath = None
        self.bed_fpath = None  # with genomic features
        self.capture_bed_fpath = None  # w/o genomic features
        self.qualimap_bed_fpath = None
        self.padded_bed_fpath = None

        self.gene_keys_set = set()  # set of pairs (gene_name, chrom)
        self.gene_keys_list = list()  # list of pairs (gene_name, chrom)
        self.regions_num = None

        self.bases_num = None
        self.fraction = None

        if bed_fpath:
            debug('Using target BED file ' + bed_fpath)
            self.is_wgs = False
            verify_bed(bed_fpath, is_critical=True)
            self.original_bed_fpath = bed_fpath
            self._make_target_bed(bed_fpath, work_dir, output_dir, padding=padding,
                is_debug=is_debug, fai_fpath=fai_fpath, genome=genome, reannotate=reannotate)
        else:
            debug('No input BED. Assuming whole genome. For region-based reports, analysing RefSeq CDS.')
            self.is_wgs = True
            self._make_wgs_regions_file(work_dir, genome=genome)

    def get_capture_bed(self):
        if not self.is_wgs:
            if self.bed.field_count() > ebl.BedCols.FEATURE:
                return self.bed.filter(lambda x: x[ebl.BedCols.FEATURE] == 'capture')
            else:
                return self.bed
        else:
            return None

    def _make_target_bed(self, bed_fpath, work_dir, output_dir, is_debug,
                         padding=None, fai_fpath=None, genome=None, reannotate=False):
        clean_target_bed_fpath = intermediate_fname(work_dir, bed_fpath, 'clean')
        if not can_reuse(clean_target_bed_fpath, bed_fpath):
            debug()
            debug('Cleaning target BED file...')
            bed = BedTool(bed_fpath)
            if bed.field_count() > 4:
                bed = bed.cut(range(4))
            bed = bed\
                .filter(lambda x: x.chrom and not any(x.chrom.startswith(e) for e in ['#', ' ', 'track', 'browser']))\
                .remove_invalid()
            with file_transaction(work_dir, clean_target_bed_fpath) as tx:
                bed.saveas(tx)
            debug('Saved to ' + clean_target_bed_fpath)
            verify_file(clean_target_bed_fpath, is_critical=True)

        sort_target_bed_fpath = intermediate_fname(work_dir, clean_target_bed_fpath, 'sorted')
        if not can_reuse(sort_target_bed_fpath, clean_target_bed_fpath):
            debug()
            debug('Sorting target BED file...')
            sort_target_bed_fpath = sort_bed(clean_target_bed_fpath, output_bed_fpath=sort_target_bed_fpath, fai_fpath=fai_fpath)
            debug('Saved to ' + sort_target_bed_fpath)
            verify_file(sort_target_bed_fpath, is_critical=True)

        if genome in ebl.SUPPORTED_GENOMES:
            ann_target_bed_fpath = intermediate_fname(work_dir, sort_target_bed_fpath, 'ann_plus_features')
            if not can_reuse(ann_target_bed_fpath, sort_target_bed_fpath):
                debug()
                if BedTool(sort_target_bed_fpath).field_count() == 3 or reannotate:
                        debug('Annotating target BED file and collecting overlapping genome features')
                        overlap_with_features(sort_target_bed_fpath, ann_target_bed_fpath, work_dir=work_dir,
                             genome=genome, extended=True, reannotate=reannotate, only_canonical=True)
                else:
                    debug('Overlapping with genomic features:')
                    overlap_with_features(sort_target_bed_fpath, ann_target_bed_fpath, work_dir=work_dir,
                         genome=genome, extended=True, only_canonical=True)
                debug('Saved to ' + ann_target_bed_fpath)
                verify_file(ann_target_bed_fpath, is_critical=True)
        else:
            ann_target_bed_fpath = sort_target_bed_fpath

        final_clean_target_bed_fpath = intermediate_fname(work_dir, ann_target_bed_fpath, 'clean')
        if not can_reuse(final_clean_target_bed_fpath, ann_target_bed_fpath):
            bed = BedTool(ann_target_bed_fpath).remove_invalid()
            with file_transaction(work_dir, final_clean_target_bed_fpath) as tx:
                bed.saveas(tx)
            verify_file(final_clean_target_bed_fpath, is_critical=True)

        self.bed_fpath = final_clean_target_bed_fpath
        self.bed = BedTool(self.bed_fpath)
        
        self.capture_bed_fpath = add_suffix(join(output_dir, basename(bed_fpath)), 'clean_sorted_ann')
        if not can_reuse(self.capture_bed_fpath, self.bed_fpath):
            with file_transaction(work_dir, self.capture_bed_fpath) as tx:
                self.get_capture_bed().saveas(tx)

        gene_key_set, gene_key_list = get_genes_from_bed(bed_fpath)
        self.gene_keys_set = gene_key_set
        self.gene_keys_list = gene_key_list
        self.regions_num = self.get_capture_bed().count()

        self._make_qualimap_bed(work_dir)
        if padding:
            self._make_padded_bed(work_dir, fai_fpath, padding)

    def _make_padded_bed(self, work_dir, fai_fpath, padding):
        if self.is_wgs:
            return None

        self.padded_bed_fpath = intermediate_fname(work_dir, self.capture_bed_fpath, 'padded')
        if can_reuse(self.padded_bed_fpath, self.capture_bed_fpath):
            return BedTool(self.padded_bed_fpath)

        padded_bed = self.bed.slop(b=padding, g=fai_fpath).sort().merge()
        with file_transaction(work_dir, self.padded_bed_fpath) as tx:
            padded_bed.saveas(tx)
        verify_file(self.padded_bed_fpath, is_critical=True)
        return BedTool(self.padded_bed_fpath)

    def _make_qualimap_bed(self, work_dir):
        if self.is_wgs:
            return None

        self.qualimap_bed_fpath = intermediate_fname(work_dir, self.capture_bed_fpath, 'qualimap_ready')
        if can_reuse(self.qualimap_bed_fpath, self.capture_bed_fpath):
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
        verify_file(self.qualimap_bed_fpath, is_critical=True)
        return self.qualimap_bed_fpath

    def _make_wgs_regions_file(self, work_dir, genome=None):
        self.wgs_bed_fpath = join(work_dir, 'targqc_features_to_report.bed')
        if can_reuse(self.wgs_bed_fpath, ebl.ensembl_gtf_fpath(genome)):
            return self.wgs_bed_fpath

        chr_order = reference_data.get_chrom_order(genome or cfg.genome)

        r_by_tx_by_gene = OrderedDefaultDict(lambda: defaultdict(list))
        all_features = ebl.get_all_features(genome or cfg.genome, high_confidence=True)

        debug('Select best transcript to report')
        for r in all_features:
            if r[ebl.BedCols.FEATURE] != 'gene':
                gene = r[ebl.BedCols.HUGO]
                tx = r[ebl.BedCols.ENSEMBL_ID]
                r_by_tx_by_gene[gene][tx].append(r.fields)

        with file_transaction(work_dir, self.wgs_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for gname, r_by_tx in r_by_tx_by_gene.items():
                    all_tx = (x for xx in r_by_tx.values() for x in xx if x[ebl.BedCols.FEATURE] == 'transcript')
                    tx_sorted_list = [x[ebl.BedCols.ENSEMBL_ID] for x in sorted(all_tx, key=tx_priority_sort_key)]
                    if not tx_sorted_list:
                        continue
                    tx_id = tx_sorted_list[0]
                    for r in sorted(r_by_tx[tx_id], key=get_sort_key(chr_order)):
                        out.write('\t'.join(str(f) for f in r) + '\n')
        return self.wgs_bed_fpath