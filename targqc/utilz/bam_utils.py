from targqc.utilz.logger import critical, err
from targqc.utilz.file_utils import verify_file, adjust_path


def verify_bam(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    logfn = critical if is_critical else err
    if not fpath.endswith('.bam'):
        logfn('The file ' + fpath + ' is supposed to be BAM but does not have the .bam '
            'extension. Please, make sure you pass proper file.')
        return None

    # TODO: check if binary

    return fpath


