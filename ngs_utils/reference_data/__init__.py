from os.path import dirname, join, abspath, splitext, isfile

from ngs_utils.file_utils import verify_file, adjust_path, verify_dir
from ngs_utils.logger import critical, debug


SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr21', 'GRCh37', 'mm10']


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + str(genome) + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))

def _get(relative_path, genome=None):
    if genome:
        check_genome(genome)
    else:
        genome = ''
    relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    return path


######################
######## FAI #########
def get_fai(genome):
    return _get(join('fai', '{genome}.fa.fai'), genome)

def get_chrom_lengths(genome=None, fai_fpath=None):
    assert genome or fai_fpath, 'One of genome or fai_fpath should be not None: ' \
                                'genome=' + str(genome) + ' fai_fpath=' + str(fai_fpath)

    if not fai_fpath:
        check_genome(genome)
        fai_fpath = get_fai(genome)
    else:
        fai_fpath = verify_file(fai_fpath, is_critical=True)
        if not fai_fpath.endswith('.fai') and not fai_fpath.endswith('.fa'):
            critical('Error: .fai or .fa is accepted.')

    chr_lengths = []

    if fai_fpath.endswith('.fa'):
        debug('Reading genome sequence (.fa) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append((chrom, len(record.seq)))

    else:
        debug('Reading genome index file (.fai) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], int(line.split()[1])
                    chr_lengths.append((chrom, length))

    return chr_lengths

def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order

def ucsc_to_ensembl(genome):
    """ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.ucscToEnsembl.tsv
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.ucscToEnsembl.tsv
    """
    return _get(join('fai', '{genome}.ucscToEnsembl.tsv'), genome)
