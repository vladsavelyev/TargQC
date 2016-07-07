import os
import random
import gzip
from itertools import izip, product
from os.path import splitext, dirname, join, basename

from Utils import sambamba
from Utils.call_process import run
from Utils.file_utils import open_gzipsafe, file_transaction, file_exists, intermediate_fname, verify_file, add_suffix, \
    splitext_plus, safe_mkdir, which
from Utils.logger import critical, debug, info, warn

import targqc.config as cfg


def count_records_in_fastq(fastq_fpath):
    return sum(1 for _ in open_gzipsafe(fastq_fpath)) / 4


def downsample(work_dir, sample_name, output_dir, fastq_left_fpath, fastq_right_fpath, N_pairs,
               suffix=None, num_read_pairs=None):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    """
    sample_name = sample_name or splitext(''.join(lc if lc == rc else '' for lc, rc in izip(fastq_left_fpath, fastq_right_fpath)))[0]

    l_out_fpath = join(output_dir, add_suffix(basename(fastq_left_fpath), suffix or 'subset'))
    r_out_fpath = join(output_dir, add_suffix(basename(fastq_right_fpath), suffix or 'subset'))
    if cfg.reuse_intermediate and verify_file(l_out_fpath, silent=True) and verify_file(r_out_fpath, silent=True):
        debug(l_out_fpath + ' and ' + r_out_fpath + ' exist, reusing.')
        return l_out_fpath, r_out_fpath

    info('Processing ' + sample_name)
    N_pairs = int(N_pairs)
    info(sample_name + ': getting number of reads in fastq...')
    num_read_pairs = num_read_pairs or count_records_in_fastq(fastq_left_fpath)
    if num_read_pairs > LIMIT:
        info(sample_name + ' the number of reads is higher than ' + str(LIMIT) +
             ', sampling from only first ' + str(LIMIT))
        num_read_pairs = LIMIT
    info(sample_name + ': ' + str(num_read_pairs) + ' reads')
    if num_read_pairs < N_pairs:
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
            rec_no += 1
            if written_records % 10000 == 0:
                info(sample_name + ': written ' + str(written_records) + ', rec_no ' + str(rec_no))
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


def align(work_dir, sample_name, l_fpath, r_fpath, bwa, samblaster, samtools, sambabma, bwa_prefix, dedup=True, threads=1):
    info('Aligning reads')
    bam_fpath = join(work_dir, sample_name + '_downsampled.bam')
    tmp_dirpath = join(work_dir, 'sambamba_tmp_dir')
    safe_mkdir(tmp_dirpath)

    bwa_cmdline = ('{bwa} mem -t {threads} -v 2 {bwa_prefix} {l_fpath} {r_fpath} | ' +
                   ('{samblaster} | ' if dedup else '') +
                   '{samtools} view --threads {threads} -b -S -u - | '
                   '{sambabma} sort -t {threads} --tmpdir {tmp_dirpath} '
                        '-o {bam_fpath} /dev/stdin').format(**locals())
    run(bwa_cmdline, output_fpath=bam_fpath, stdout_to_outputfile=False, reuse=cfg.reuse_intermediate)
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


def markdup_sam(in_sam_fpath, samblaster=None):
    """Perform non-stream based deduplication of SAM input files using samblaster.
    """
    if not samblaster:
        samblaster = 'samblaster'
        if which(samblaster) is None:
            warn('No samblaster, can\'t mark duplicates.')
            return None

    out_sam_fpath = add_suffix(in_sam_fpath, 'markdup')
    cmdline = '{samblaster} -i {in_sam_fpath} -o {out_sam_fpath}'.format(**locals())
    run(cmdline, output_fpath=out_sam_fpath, stdout_to_outputfile=False, reuse=cfg.reuse_intermediate)
    return out_sam_fpath


