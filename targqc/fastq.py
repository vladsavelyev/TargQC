import os
import random
import gzip
from os.path import splitext, dirname, join, basename, isfile, abspath

from ngs_utils import sambamba
from ngs_utils.bam_utils import verify_bam
from ngs_utils.call_process import run
from ngs_utils.file_utils import open_gzipsafe, file_transaction, verify_file, add_suffix, safe_mkdir, which, can_reuse
from ngs_utils.logger import critical, debug, info, warn, err

try:
    from itertools import izip as zip
except ImportError:
    pass


def make_downsampled_fpath(work_dir, fastq_fpath):
    return join(work_dir, add_suffix(basename(fastq_fpath), 'subset'))

def make_pair_counts_fpath(work_dir):
    return join(work_dir, 'full_pairs_count.txt')

def make_bam_fpath(work_dir):
    return join(work_dir, 'downsampled.bam')


def proc_fastq(samples, parall_view, work_dir, bwa_prefix, downsample_to, num_pairs_by_sample=None, dedup=True):
    num_pairs_by_sample = num_pairs_by_sample or dict()
    if downsample_to:
        # Read pairs counts
        debug()
        if all(s.name in num_pairs_by_sample for s in samples):
            debug('Using read pairs counts extracted from FastQC reports')
        elif all(can_reuse(make_pair_counts_fpath(s.work_dir), s.l_fpath) for s in samples):
            debug('Reusing pairs counts, reading from files')
            num_pairs_by_sample = {s.name: int(open(make_pair_counts_fpath(s.work_dir)).read().strip()) for s in samples}
        else:
            info('Counting read pairs')
            num_pairs = parall_view.run(count_read_pairs, [[s.name, s.work_dir, s.l_fpath] for s in samples])
            num_pairs_by_sample = {s.name: pairs_count for s, pairs_count in zip(samples, num_pairs)}

        # Downsampling
        debug()
        if all(can_reuse(make_downsampled_fpath(s.work_dir, s.l_fpath), s.l_fpath) and
               can_reuse(make_downsampled_fpath(s.work_dir, s.r_fpath), s.r_fpath) for s in samples):
            debug('Reusing downsampled FastQ')
            for s in samples:
                s.l_fpath = make_downsampled_fpath(s.work_dir, s.l_fpath)
                s.r_fpath = make_downsampled_fpath(s.work_dir, s.r_fpath)
        else:
            if isinstance(downsample_to, float):
                info('Downsampling FastQ to ' + str(float(downsample_to)) + ' fraction of reads')
            else:
                info('Downsampling FastQ to ' + str(int(downsample_to)) + ' read pairs')
            fastq_pairs = parall_view.run(downsample,
                [[s.work_dir, s.name, s.l_fpath, s.r_fpath, downsample_to, num_pairs_by_sample.get(s.name)]
                 for s in samples])
            for s, (l_r, r_r) in zip(samples, fastq_pairs):
                s.l_fpath = l_r
                s.r_fpath = r_r
    else:
        info('Skipping downsampling')

    debug()
    if all(can_reuse(make_bam_fpath(s.work_dir), [s.l_fpath, s.r_fpath]) for s in samples):
        debug('All downsampled BAM exists, reusing')
        for s in samples:
            s.bam = make_bam_fpath(s.work_dir)
    else:
        bwa = which('bwa')
        if not isfile(bwa):
            critical('BWA not found under ' + bwa)
        smb = sambamba.get_executable()
        if not (bwa and smb):
            if not bwa:         err('Error: bwa is required for the alignment pipeline')
            if not smb:         err('Error: sambamba is required for the alignment pipeline')
            critical('Tools required for alignment not found')
        info('Aligning reads to the reference')
        bam_fpaths = parall_view.run(align,
            [[s.work_dir, s.name, s.l_fpath, s.r_fpath, bwa, smb, bwa_prefix, dedup, parall_view.cores_per_job]
             for s in samples])

        bam_fpaths = map(verify_bam, bam_fpaths)
        if len(bam_fpaths) < len(samples):
            critical('Some samples were not aligned successfully.')
        for bam, s in zip(bam_fpaths, samples):
            s.bam = bam

    return num_pairs_by_sample


def count_read_pairs(s_name, work_dir, fastq_fpath):
    from ngs_utils.logger import info

    pairs_counts_fpath = make_pair_counts_fpath(work_dir)
    if can_reuse(pairs_counts_fpath, fastq_fpath):
        with open(pairs_counts_fpath) as f:
            return int(f.read().strip())
    else:
        info('Counting read pairs in ' + s_name + ', writing to ' + pairs_counts_fpath)
        pairs_number = _count_records_in_fastq(fastq_fpath)
        with open(pairs_counts_fpath, 'w') as out:
            out.write(str(pairs_number))
        return pairs_number


def _count_records_in_fastq(fastq_fpath):
    return sum(1 for _ in open_gzipsafe(fastq_fpath)) / 4


def downsample(work_dir, sample_name, fastq_left_fpath, fastq_right_fpath, downsample_to, num_pairs=None):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    """
    sample_name = sample_name or splitext(''.join(lc if lc == rc else '' for lc, rc in zip(fastq_left_fpath, fastq_right_fpath)))[0]

    l_out_fpath = make_downsampled_fpath(work_dir, fastq_left_fpath)
    r_out_fpath = make_downsampled_fpath(work_dir, fastq_right_fpath)
    if can_reuse(l_out_fpath, [fastq_left_fpath, fastq_right_fpath]):
        return l_out_fpath, r_out_fpath

    info('Processing ' + sample_name)
    if num_pairs is None:
        info(sample_name + ': counting number of reads in fastq...')
        num_pairs = _count_records_in_fastq(fastq_left_fpath)
    if num_pairs > LIMIT:
        info(sample_name + ' the number of reads is higher than ' + str(LIMIT) +
             ', sampling from only first ' + str(LIMIT))
        num_pairs = LIMIT
    info(sample_name + ': ' + str(num_pairs) + ' reads')
    num_downsample_pairs = int(downsample_to * num_pairs) if isinstance(downsample_to, float) else downsample_to
    if num_pairs <= num_downsample_pairs:
        info(sample_name + ': and it is less than ' + str(num_downsample_pairs) + ', so no downsampling.')
        return fastq_left_fpath, fastq_right_fpath
    else:
        info(sample_name + ': downsampling to ' + str(num_downsample_pairs))
        rand_records = sorted(random.sample(range(num_pairs), num_downsample_pairs))

    info('Opening ' + fastq_left_fpath)
    fh1 = open_gzipsafe(fastq_left_fpath)
    info('Opening ' + fastq_right_fpath)
    fh2 = open_gzipsafe(fastq_right_fpath) if fastq_right_fpath else None

    out_files = (l_out_fpath, r_out_fpath) if r_out_fpath else (l_out_fpath,)

    written_records = 0
    with file_transaction(work_dir, out_files) as tx_out_files:
        if isinstance(tx_out_files, str):
            tx_out_f1 = tx_out_files
        else:
            tx_out_f1, tx_out_f2 = tx_out_files
        info('Opening ' + str(tx_out_f1) + ' to write')
        sub1 = open_gzipsafe(tx_out_f1, "w")
        info('Opening ' + str(tx_out_f2) + ' to write')
        sub2 = open_gzipsafe(tx_out_f2, "w") if r_out_fpath else None
        rec_no = -1
        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fh1.readline()
                if fh2:
                    for i in range(4): fh2.readline()
            for i in range(4):
                sub1.write(fh1.readline())
                if sub2:
                    sub2.write(fh2.readline())
            written_records += 1
            if written_records % 10000 == 0:
                info(sample_name + ': written ' + str(written_records) + ', rec_no ' + str(rec_no + 1))
            if rec_no > num_pairs:
                info(sample_name + ' reached the limit of ' + str(num_pairs), ' read lines, stopping.')
                break
        info(sample_name + ': done, written ' + str(written_records) + ', rec_no ' + str(rec_no))
        fh1.close()
        sub1.close()
        if fastq_right_fpath:
            fh2.close()
            sub2.close()

    info(sample_name + ': done downsampling, saved to ' + l_out_fpath + ' and ' + r_out_fpath + ', total ' + str(written_records) + ' paired reads written')
    return l_out_fpath, r_out_fpath


def align(work_dir, sample_name, l_fpath, r_fpath, bwa, smb, bwa_prefix, dedup=True, threads=1):
    info('Running bwa to align reads...')
    bam_fpath = make_bam_fpath(work_dir)
    if can_reuse(bam_fpath, [l_fpath, r_fpath]):
        return bam_fpath

    tmp_dirpath = join(work_dir, 'sambamba_tmp_dir')
    safe_mkdir(tmp_dirpath)

    bwa_cmdline = ('{bwa} mem -t {threads} -v 2 {bwa_prefix} {l_fpath} {r_fpath} | ' +
                   '{smb} view /dev/stdin -t {threads} -f bam -S -o - | ' +
                   '{smb} sort /dev/stdin -t {threads} --tmpdir {tmp_dirpath} -o {bam_fpath}').format(**locals())
    run(bwa_cmdline, output_fpath=bam_fpath, stdout_to_outputfile=False)

    if dedup:
        dedup_bam_fpath = add_suffix(bam_fpath, 'dedup')
        dedup_cmdl = '{smb} markdup -t {threads} {bam_fpath} {dedup_bam_fpath}'.format(**locals())
        run(dedup_cmdl, output_fpath=dedup_bam_fpath, stdout_to_outputfile=False)
        verify_bam(dedup_bam_fpath)
        os.rename(dedup_bam_fpath, bam_fpath)

    sambamba.index_bam(bam_fpath)

# samtools view -b -S -u - |
# sambamba sort -N -t 8 -m 682M --tmpdir /Molly/saveliev/cancer-dream-syn3/work/align/syn3-normal/split/tx/tmpwdXndE/syn3-normal-sort-1_20000000-sorttmp-full
    # -o /Molly/saveliev/cancer-dream-syn3/work/align/syn3-normal/split/tx/tmpwdXndE/syn3-normal-sort-1_20000000.bam
    # /dev/stdin

    # if dedup:
    #     info()
    #     info('Calling SamBlaster to mark duplicates')
    #     markdup_sam_fpath = markdup_sam(sam_fpath, samblaster)
    #     if markdup_sam_fpath:
    #         sam_fpath = markdup_sam_fpath
    # info()

    # info('Converting to BAM')
    # cmdline = sambamba.get_executable() + ' view -t {threads} -S -f bam {sam_fpath}'.format(**locals())
    # run(cmdline, output_fpath=bam_fpath, reuse=cfg.reuse_intermediate)
    #
    # info()
    # info('Sorting BAM')
    # prefix = splitext(sorted_bam_fpath)[0]
    # cmdline = sambamba.get_executable() + ' sort -t {threads} {bam_fpath} -o {sorted_bam_fpath}'.format(**locals())
    # run(cmdline, output_fpath=sorted_bam_fpath, stdout_to_outputfile=False, reuse=cfg.reuse_intermediate)

    return bam_fpath


LIMIT = 500*1000*1000


# def markdup_sam(in_sam_fpath, samblaster=None, reuse=False):
#     """Perform non-stream based deduplication of SAM input files using samblaster.
#     """
#     if not samblaster:
#         samblaster = 'samblaster'
#         if which(samblaster) is None:
#             warn('No samblaster, can\'t mark duplicates.')
#             return None
#
#     out_sam_fpath = add_suffix(in_sam_fpath, 'markdup')
#     cmdline = '{samblaster} -i {in_sam_fpath} -o {out_sam_fpath}'.format(**locals())
#     run(cmdline, output_fpath=out_sam_fpath, stdout_to_outputfile=False, reuse=reuse)
#     return out_sam_fpath


