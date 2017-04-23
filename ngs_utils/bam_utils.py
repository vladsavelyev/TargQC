from ngs_utils.file_utils import verify_file, adjust_path
from ngs_utils.logger import critical, err


def verify_bam(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    logfn = critical if is_critical else err
    if not fpath.endswith('.bam'):
        logfn('The file ' + fpath + ' is supposed to be BAM but does not have the .bam '
            'extension. Please, make sure you pass proper file.')
        return None

    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda baitiki: bool(baitiki.translate(None, textchars))
    if not is_binary_string(open(fpath).read(3)):
        logfn('The BAM file ' + fpath + ' must be a binary file.')
        return None

    return fpath


