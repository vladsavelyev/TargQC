from os.path import abspath, join

from GeneAnnotation.annotate_bed import annotate
from Utils.bed_utils import filter_bed_with_gene_set, get_gene_keys_from_bed, calc_region_number
from Utils.file_utils import add_suffix
from Utils.logger import debug, warn, critical
from targqc import config as cfg


class Target:
    def __init__(self, bed_fpath=None, genes_fpath=None):
        self.bed_fpath = bed_fpath
        self.original_bed_fpath = bed_fpath
        self.type = 'Regional' if bed_fpath else 'WGS'

        self.genes_fpath = abspath(genes_fpath) if genes_fpath else None
        self.gene_keys_set = set()   # set of pairs (gene_name, chrom)
        self.gene_keys_list = []    # list of pairs (gene_name, chrom)
        self.genes_not_in_refseq = []

        self.regions_num = None
        self.bases_num = None
        self.fraction = None

    def extract_gene_names_and_filter_exons(self, work_dir, features_bed_fpath, reuse=False):
        debug()
        debug('Getting gene list')

        # if genes_fpath:
        #     with open(genes_fpath) as f:
        #         gene_key_list = [g.strip() for g in f.read().split('\n') if g]
        #         gene_key_set = set(gene_key_list)
        #     info('Using genes from ' + genes_fpath + ', filtering exons and amplicons with this genes.')
        #     if target_bed:
        #         target_bed = filter_bed_with_gene_set(cnf, target_bed, gene_key_set)
        #     if exons_bed:
        #         exons_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_key_set)
        #         exons_no_genes_bed = filter_bed_with_gene_set(cnf, exons_no_genes_bed, gene_key_set)
        # else:

        gene_key_set, gene_key_list = get_gene_keys_from_bed(self.bed_fpath)

        debug('Using genes from the target ' + self.bed_fpath)

        debug('Trying filtering exons with these ' + str(len(gene_key_list)) + ' genes.')
        features_filt_bed_fpath, genes_in_refseq = filter_bed_with_gene_set(work_dir, features_bed_fpath, gene_key_set, suffix='target_genes_1st_round', reuse=reuse)
        if not genes_in_refseq:
            debug()
            warn('No gene symbols from the target BED file was found in the RefSeq features. Re-annotating target...')
            self.bed_fpath = annotate(self.bed_fpath, features_bed_fpath, add_suffix(self.bed_fpath, 'ann'), reuse=reuse, genome=cfg.genome)
            #info('Merging regions within genes...')
            #target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)
            # debug('Sorting amplicons_bed by (chrom, gene_name, start)')
            # target_bed = sort_bed(work_dir, target_bed)
            debug('Getting gene names again...')
            gene_key_set, gene_key_list = get_gene_keys_from_bed(self.bed_fpath)
            debug()
            debug('Using genes from the new amplicons list, filtering features with this genes again.')
            features_filt_bed_fpath, genes_in_refseq = filter_bed_with_gene_set(work_dir, features_bed_fpath, gene_key_set, suffix='target_genes_2nd_round', reuse=reuse)
            if not genes_in_refseq:
                critical('No gene symbols from the target BED file was found in the RefSeq features.')

        features_bed_fpath = features_filt_bed_fpath

        self.gene_keys_set = gene_key_set
        self.gene_keys_list = gene_key_list
        self.genes_not_in_refseq = gene_key_set - genes_in_refseq

        self.regions_num = calc_region_number(self.bed_fpath)

        return features_bed_fpath

