import os
import random
import gzip
from itertools import izip, product
from os.path import splitext, dirname, join, basename, isfile

from Utils import sambamba
from Utils.bam_utils import verify_bam
from Utils.call_process import run
from Utils.file_utils import open_gzipsafe, file_transaction, file_exists, intermediate_fname, verify_file, add_suffix, \
    splitext_plus, safe_mkdir, which
from Utils.logger import critical, debug, info, warn, err


def proc_fastq(samples, parall_view, work_dir, bwa_prefix, downsample_pairs_num, num_reads_by_sample=None,
               dedup=True, reuse=False):
    num_reads_by_sample = num_reads_by_sample or dict()
    if downsample_pairs_num:
        if all(s.name in num_reads_by_sample for s in samples):
            pass
        else:
            info('Counting read numbers')
            read_counts = parall_view.run(count_reads, [[s.name, s.work_dir, s.l_fpath, reuse] for s in samples])
            for s, read_count in zip(samples, read_counts):
                num_reads_by_sample[s.name] = read_count

        info('Downsampling the reads to ' + str(int(downsample_pairs_num)))
        fastq_pairs = parall_view.run(downsample,
            [[s.work_dir, s.name, s.work_dir, s.l_fpath, s.r_fpath, downsample_pairs_num, 'subset']
             for s in samples])
        for s, (l_r, r_r) in zip(samples, fastq_pairs):
            s.l_fpath = l_r
            s.r_fpath = r_r
    else:
        info('Skipping downsampling')

    bwa = which('bwa')
    samtools = which('samtools')
    sb = sambamba.get_executable()
    if not (bwa and samtools and sb):
        if not bwa:         err('Error: bwa is required for alignment')
        if not samtools:    err('Error: samtools is required for alignment')
        if not sb:          err('Error: sambamba is required')
        critical()
    info()
    info('Aligning reads to the reference')
    bam_fpaths = parall_view.run(align,
        [[s.work_dir, s.name, s.l_fpath, s.r_fpath, bwa, samtools, sb, bwa_prefix, dedup, parall_view.cores_per_job, reuse]
         for s in samples])

    bam_fpaths = map(verify_bam, bam_fpaths)
    if len(bam_fpaths) < len(samples):
        critical('Some samples were not aligned successfully.')
    for bam, s in zip(bam_fpaths, samples):
        s.bam = bam

    return num_reads_by_sample


def count_reads(s_name, work_dir, fastq_fpath, reuse=False):
    from os.path import join, isfile
    from Utils.file_utils import verify_file
    from Utils.logger import info

    read_counts_fpath = join(work_dir, 'original_read_count.txt')
    if reuse and isfile(read_counts_fpath) and verify_file(read_counts_fpath, cmp_date_fpath=fastq_fpath):
        with open(read_counts_fpath) as f:
            return int(f.read().strip())
    else:
        info('Counting read numbers in ' + s_name + ', writing to ' + read_counts_fpath)
        pairs_number = count_records_in_fastq(fastq_fpath)
        read_number = pairs_number * 2
        with open(read_counts_fpath, 'w') as out:
            out.write(str(read_number))
        return read_number


def count_records_in_fastq(fastq_fpath):
    return sum(1 for _ in open_gzipsafe(fastq_fpath)) / 4


def downsample(work_dir, sample_name, output_dir, fastq_left_fpath, fastq_right_fpath, N_pairs,
               suffix=None, num_read_pairs=None, reuse=False):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    """
    sample_name = sample_name or splitext(''.join(lc if lc == rc else '' for lc, rc in izip(fastq_left_fpath, fastq_right_fpath)))[0]

    l_out_fpath = join(output_dir, add_suffix(basename(fastq_left_fpath), suffix or 'subset'))
    r_out_fpath = join(output_dir, add_suffix(basename(fastq_right_fpath), suffix or 'subset'))
    if reuse and verify_file(l_out_fpath, silent=True) and verify_file(r_out_fpath, silent=True):
        debug(l_out_fpath + ' and ' + r_out_fpath + ' exist, reusing.')
        return l_out_fpath, r_out_fpath

    info('Processing ' + sample_name)
    N_pairs = int(N_pairs)
    if num_read_pairs is None:
        info(sample_name + ': counting number of reads in fastq...')
        num_read_pairs = count_records_in_fastq(fastq_left_fpath)
    if num_read_pairs > LIMIT:
        info(sample_name + ' the number of reads is higher than ' + str(LIMIT) +
             ', sampling from only first ' + str(LIMIT))
        num_read_pairs = LIMIT
    info(sample_name + ': ' + str(num_read_pairs) + ' reads')
    if num_read_pairs <= N_pairs:
        info(sample_name + ': and it is less than ' + str(N_pairs) + ', so no downsampling.')
        return fastq_left_fpath, fastq_right_fpath
    else:
        info(sample_name + ': downsampling to ' + str(N_pairs))
        rand_records = sorted(random.sample(xrange(num_read_pairs), N_pairs))

    info('Opening ' + fastq_left_fpath)
    fh1 = open_gzipsafe(fastq_left_fpath)
    info('Opening ' + fastq_right_fpath)
    fh2 = open_gzipsafe(fastq_right_fpath) if fastq_right_fpath else None

    out_files = (l_out_fpath, r_out_fpath) if r_out_fpath else (l_out_fpath,)

    written_records = 0
    with file_transaction(work_dir, out_files) as tx_out_files:
        if isinstance(tx_out_files, basestring):
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
            if rec_no > num_read_pairs:
                info(sample_name + ' reached the limit of ' + str(num_read_pairs), ' read lines, stopping.')
                break
        info(sample_name + ': done, written ' + str(written_records) + ', rec_no ' + str(rec_no))
        fh1.close()
        sub1.close()
        if fastq_right_fpath:
            fh2.close()
            sub2.close()

    info(sample_name + ': done downsampling, saved to ' + l_out_fpath + ' and ' + r_out_fpath + ', total ' + str(written_records) + ' paired reads written')
    return l_out_fpath, r_out_fpath


def align(work_dir, sample_name, l_fpath, r_fpath, bwa, samtools, smb, bwa_prefix, dedup=True, threads=1, reuse=False):
    info('Running bwa to align reads...')
    bam_fpath = join(work_dir, sample_name + '_downsampled.bam')
    if isfile(bam_fpath) and verify_bam(bam_fpath) and reuse:
        debug(bam_fpath + ' exists, reusing')
        return bam_fpath

    tmp_dirpath = join(work_dir, 'sambamba_tmp_dir')
    safe_mkdir(tmp_dirpath)

    bwa_cmdline = ('{bwa} mem -t {threads} -v 2 {bwa_prefix} {l_fpath} {r_fpath} | ' +
                   '{smb} view /dev/stdin -t {threads} -f bam -S -o - | ' +
                   '{smb} sort /dev/stdin -t {threads} --tmpdir {tmp_dirpath} -o {bam_fpath}').format(**locals())
    run(bwa_cmdline, output_fpath=bam_fpath, stdout_to_outputfile=False, reuse=reuse)

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


