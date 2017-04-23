import os
import subprocess
import traceback
from os.path import join, dirname, abspath, basename, isfile, getmtime
import sys

from pybedtools import BedTool

from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file, splitext_plus, which, can_reuse
from ngs_utils.logger import debug, warn, err, critical


def get_executable():
    sys_path = which('sambamba')
    if not sys_path:
        critical('Error: sambamba executable is not found')
    return sys_path


def index_bam(bam_fpath, sambamba=None, samtools=None):
    sambamba = sambamba or get_executable()
    indexed_bam = bam_fpath + '.bai'
    if not can_reuse(indexed_bam, cmp_f=bam_fpath, silent=True):
        # info('Indexing BAM, writing ' + indexed_bam + '...')
        cmdline = '{sambamba} index {bam_fpath}'.format(**locals())
        res = run(cmdline, output_fpath=indexed_bam, stdout_to_outputfile=False, stdout_tx=False)
        # if not isfile(indexed_bam) or getmtime(indexed_bam) < getmtime(bam_fpath):
        #     samtools = samtools or get_system_path(cnf, 'samtools')
        #     cmdline = '{samtools} index {bam_fpath}'.format(**locals())
        #     call(cnf, cmdline)
    # else:
    #     debug('Actual "bai" index exists.')


def call_sambamba(cmdl, bam_fpath, output_fpath=None, command_name='', no_index=False):
    if not no_index:
        index_bam(bam_fpath)
    sambamba = get_executable()
    run(sambamba + ' ' + cmdl, output_fpath=output_fpath)
    return output_fpath


def sambamba_depth(work_dir, bed, bam, depth_thresholds=None,
                   output_fpath=None, sample_name=None, threads=1):
    sample_name = sample_name or splitext_plus(basename(bam))[0]
    depth_thresholds = depth_thresholds or []
    
    if isinstance(bed, BedTool):
        bed = bed.saveas().fn
    if not output_fpath:
        output_fpath = join(work_dir,
            splitext_plus(basename(bed))[0] + '_' + sample_name + '_sambamba_depth.txt')

    if can_reuse(output_fpath, [bam, bed]):
        return output_fpath

    thresholds_str = ''.join([' -T' + str(d) for d in depth_thresholds])
    cmdline = ('depth region -F "not duplicate and not failed_quality_control" '
               '-t {threads} -L {bed} {thresholds_str} {bam}').format(**locals())

    call_sambamba(cmdline, bam_fpath=bam, output_fpath=output_fpath)
    return output_fpath


def remove_dups(bam, output_fpath):
    cmdline = 'view --format=bam -F "not duplicate" {bam}'.format(**locals())  # -F (=not) 1024 (=duplicate)
    return call_sambamba(cmdline, bam_fpath=bam, output_fpath=output_fpath, command_name='not_duplicate')


def count_in_bam(work_dir, bam, query, dedup=False, bed=None, use_grid=False, sample_name=None, target_name=None):
    if dedup:
        query += ' and not duplicate'
    name = 'num_' + (query.replace(' ', '_') or 'reads')
    if bed is not None and isinstance(bed, BedTool):
        bed = bed.saveas().fn
    if bed is not None:
        target_name = target_name or ('target_' + basename(bed))
        name += '_on_' + target_name

    sample_name = sample_name or basename(bam)
    output_fpath = join(work_dir, sample_name + '_' + name)

    if can_reuse(output_fpath, cmp_f=bam):
        pass
    else:
        cmdline = 'view -c -F "{query}" {bam}'.format(**locals())
        if bed is not None:
            cmdline += ' -L ' + bed

        call_sambamba(cmdline, bam_fpath=bam, output_fpath=output_fpath, command_name=name)

    with open(output_fpath) as f:
        return int(f.read().strip())


def number_of_reads(work_dir, bam, dedup=False, use_grid=False, sample_name=None, reuse=False):
    return count_in_bam(work_dir, bam, '', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_mapped_reads(work_dir, bam, dedup=False, use_grid=False, sample_name=None, reuse=False):
    return count_in_bam(work_dir, bam, 'not unmapped', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_properly_paired_reads(work_dir, bam, dedup=False, use_grid=False, sample_name=None, reuse=False):
    return count_in_bam(work_dir, bam, 'proper_pair', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_dup_reads(work_dir, bam, use_grid=False, sample_name=None, reuse=False):
    return count_in_bam(work_dir, bam, 'duplicate', use_grid=use_grid, sample_name=sample_name)


def number_mapped_reads_on_target(work_dir, bed, bam, dedup=False, use_grid=False, sample_name=None, target_name=None):
    return count_in_bam(work_dir, bam, 'not unmapped', dedup, bed=bed, use_grid=use_grid, sample_name=sample_name, target_name=target_name)


# def flag_stat(cnf, bam):
#     output_fpath = join(cnf.work_dir, basename(bam) + '_flag_stats')
#     cmdline = 'flagstat {bam}'.format(**locals())
#     call_sambamba(cmdline, output_fpath=output_fpath, bam_fpath=bam, command_name='flagstat')
#     stats = dict()
#     with open(output_fpath) as f:
#         lines = f.readlines()
#         for stat, fun in [('total', number_of_reads),
#                           ('duplicates', number_of_dup_reads),  # '-f 1024'
#                           ('mapped', number_of_mapped_reads),   # '-F 4'
#                           ('properly paired', number_of_properly_paired_reads)]:  # '-f 2'
#             try:
#                 val = next(l.split()[0] for l in lines if stat in l)
#             except StopIteration:
#                 warn('Cannot extract ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
#                 val = None
#             else:
#                 try:
#                     val = int(val)
#                 except ValueError:
#                     warn('Cannot parse value ' + str(val) + ' from ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
#                     val = None
#             if val is not None:
#                 stats[stat] = val
#             else:
#                 stats[stat] = fun(cnf, bam)
#     return stats
