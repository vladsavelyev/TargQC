import pybedtools

from GeneAnnotation.annotate_bed import annotate
from Utils.bam_bed_utils import remove_comments, sort_bed, count_bed_cols, cut, verify_bed
from Utils.file_utils import iterate_file, add_suffix, intermediate_fname, file_transaction
from Utils.logger import debug, info, warn


def prepare_beds(work_dir, fai_fpath=None, features_bed=None, target_bed_fpath=None, cds_bed_fpath=None, reuse=False):
    if features_bed is None and target_bed_fpath is None:
        warn('No input target BED, and no features BED in the system config specified. Not making detailed per-gene reports.')

    if target_bed_fpath:
        target_bed_fpath = verify_bed(target_bed_fpath, is_critical=True)

    if features_bed:
        features_bed = verify_bed(features_bed, is_critical=True)

    # # Features
    # features_no_genes_bed = None
    # if features_bed:
    #     # info()
    #     # info('Merging regions within genes...')
    #     # exons_bed = group_and_merge_regions_by_gene(cnf, exons_bed, keep_genes=True)
    #     #
    #     # info()
    #     # info('Sorting exons by (chrom, gene name, start)')
    #     # exons_bed = sort_bed(cnf, exons_bed)
    #
    #     debug('Filtering the features bed file to have only non-gene and no-transcript records...')
    #     features_no_genes_bed = intermediate_fname(work_dir, features_bed, 'no_genes')
    #     call_process.run('grep -vw Gene ' + features_bed + ' | grep -vw Transcript', output_fpath=features_no_genes_bed, reuse=reuse)

    if target_bed_fpath:
        debug()
        debug('Cleaning target...')
        clean_target_bed_fpath = intermediate_fname(work_dir, target_bed_fpath, 'clean')
        target_bed = pybedtools.BedTool(target_bed_fpath)
        target_bed = target_bed.filter(lambda x: x.chrom and
             not any(x.chrom.startswith(e) for e in ['#', ' ', 'track', 'browser']))
        target_bed = target_bed.remove_invalid()
        target_bed = target_bed.cut(range(3))
        with file_transaction(work_dir, clean_target_bed_fpath) as tx:
            target_bed.saveas(tx)
        debug('Saved to ' + clean_target_bed_fpath)

        debug()
        debug('Sorting target...')
        sort_target_bed_fpath = sort_bed(clean_target_bed_fpath, work_dir=work_dir, fai_fpath=fai_fpath, reuse=reuse)
        debug('Saved to ' + sort_target_bed_fpath)

        debug()
        info('Annotating target...')
        ann_target_bed_fpath = add_suffix(sort_target_bed_fpath, 'ann')
        annotate(sort_target_bed_fpath, features_bed, ann_target_bed_fpath, reuse=reuse)
        debug('Saved to ' + ann_target_bed_fpath)

        target_bed_fpath = ann_target_bed_fpath

    return target_bed_fpath, features_bed