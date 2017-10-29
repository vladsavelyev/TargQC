from os.path import isfile

from utilz.bed_utils import get_genes_from_bed, get_total_bed_size


def get_bed_genes(bed_fpath):
    gene_set, gene_list = get_genes_from_bed(bed_fpath)
    return [gn for gn in gene_list if (gn and gn != '.')]


def get_genes_from_file(genes_fpath):
    with open(genes_fpath) as f:
        gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])
    return gene_names


def get_key_genes(genome, genes_file=None, get_key_genes_file=None):
    key_gene_names = get_genes_from_file(genes_file or get_key_genes_file())
    if genome in ['mm10', 'rn6']:
        return [gn[0] + gn[1:].lower() for gn in key_gene_names]
    else:
        return key_gene_names


def get_target_genes(genome, bed_file=None, get_key_genes_file=None):
    if is_small_target(bed_file):
        return get_bed_genes(bed_file)
    else:
        return get_key_genes(genome, get_key_genes_file=get_key_genes_file)


def is_small_target(bed_file=None):
    return bed_file and isfile(bed_file) and get_total_bed_size(bed_file) < 10 * 1000 * 1000
