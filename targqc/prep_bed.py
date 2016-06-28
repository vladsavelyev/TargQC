from GeneAnnotation.annotate_bed import annotate
from Utils.bam_bed_utils import remove_comments, sort_bed, count_bed_cols, cut, verify_bed
from Utils.file_utils import iterate_file, add_suffix
from Utils.logger import debug, info, warn


def prepare_beds(work_dir, fai_fpath=None, features_bed=None, target_bed=None, seq2c_bed=None, cds_bed_fpath=None, reuse=False):
    if features_bed is None and target_bed is None:
        warn('No input target BED, and no features BED in the system config specified. Not making detailed per-gene reports.')

    if target_bed:
        target_bed = verify_bed(target_bed, is_critical=True)

    if seq2c_bed:
        seq2c_bed = verify_bed(seq2c_bed, is_critical=True)

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

    ori_target_bed_path = target_bed
    if target_bed:
        debug()
        info('Remove comments in target...')
        target_bed = remove_comments(work_dir, target_bed, reuse=reuse)

        debug()
        debug('Cutting target...')
        target_bed = cut(target_bed, 4, reuse=reuse)

        debug()
        info('Sorting target...')
        target_bed = sort_bed(target_bed, work_dir=work_dir, fai_fpath=fai_fpath, reuse=reuse)

        cols = count_bed_cols(target_bed)
        if not features_bed:
            warn('Cannot re-annotate BED file without features')
        else:
            info('Annotating target...')
            target_bed = annotate(target_bed, features_bed, add_suffix(target_bed, 'ann'), reuse=reuse)


    return features_bed, target_bed, seq2c_bed