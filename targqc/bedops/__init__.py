from os.path import abspath, join, dirname, isfile
from sys import platform as _platform

from Utils.file_utils import which
from Utils.logger import warn


def get_executable():
    if 'darwin' in _platform:
        path = abspath(join(dirname(__file__), 'bedops_osx'))
    else:
        path = abspath(join(dirname(__file__), 'bedops_lnx'))
    if isfile(path):
        return path
    else:
        sys_path = which('bedops')
        warn('BedOps executable not found in ' + path + ', using system BedOps at ' + sys_path)
        return sys_path
