import os
import shutil
import six
import subprocess
import sys
import unittest
from datetime import datetime
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, \
    expanduser, \
    realpath
from targqc.utilz.utils import is_travis
from utilz.file_utils import verify_dir, verify_file


def info(msg=''):
    sys.stdout.write(msg + '\n')

def call(cmdl, suppress_output=False):
    info(cmdl if isinstance(cmdl, six.string_types) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, six.string_types):
        return subprocess.call(cmdl, shell=True, executable='/bin/bash',
                               stdout=subprocess.DEVNULL if suppress_output else None)
    else:
        return subprocess.call(cmdl,
                               stdout=subprocess.DEVNULL if suppress_output else None)

def check_call(cmdl):
    info(cmdl if isinstance(cmdl, six.string_types) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, six.string_types):
        subprocess.check_call(cmdl, shell=True, executable='/bin/bash')
    else:
        subprocess.check_call(cmdl)

def check_output(cmdl):
    info(cmdl if isinstance(cmdl, six.string_types) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, six.string_types):
        return subprocess.check_output(cmdl, shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
    else:
        return subprocess.check_output(cmdl, stderr=subprocess.STDOUT)

def swap_output(output_path):
    if not exists(output_path):
        return

    last_changed = datetime.fromtimestamp(getctime(output_path))
    prev_output_path = output_path + '_' + last_changed.strftime('%Y_%m_%d__%H_%M_%S')
    os.rename(output_path, prev_output_path)
    return prev_output_path

def swap_prev_symlink(output_path, prev_output_path):
    prev_output_link = get_prev(output_path)
    if exists(prev_output_link):
        os.remove(prev_output_link)
    if prev_output_path:
        os.symlink(prev_output_path, prev_output_link)

def get_prev(fpath):
    return fpath + '_prev'


class BaseTestCase(unittest.TestCase):
    script = None

    data_dir = 'data'
    results_dir = 'results'
    gold_standard_dir = 'gold_standard'

    remove_work_dir_on_success = False

    def setUp(self):
        if not isdir(self.data_dir):
            os.makedirs(self.data_dir)
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _check_file_throws(self, fpath, ignore_matching_lines=None, wrapper=None, cmp_line_number_only=True,
                           check_diff=True):
        assert isfile(fpath), 'file does not exist: ' + fpath
        assert getsize(fpath) > 0, 'file is empty: ' + fpath

        if check_diff:
            cmp_fpath = None
            if isdir(self.gold_standard_dir):
                cmp_fpath = join(self.gold_standard_dir, relpath(fpath, self.results_dir))
            elif isfile(get_prev(fpath)):
                cmp_fpath = get_prev(fpath)
            assert cmp_fpath, 'cmp_file not found (no gold standard dir or *_prev file)'
            assert isfile(cmp_fpath), 'cmp_file ' + cmp_fpath + ' not found'

            cmdl = 'diff'
            if ignore_matching_lines:
                if isinstance(ignore_matching_lines, six.string_types):
                    ignore_matching_lines = [ignore_matching_lines]
                for r in ignore_matching_lines:
                    cmdl += ' -I ' + subprocess.list2cmdline([r])
            if wrapper:
                wrapper = subprocess.list2cmdline(wrapper)
                if not fpath.endswith('.gz'):
                    fpath = '<(cat ' + fpath + ' | ' + wrapper + ')'
                    cmp_fpath = '<(cat ' + cmp_fpath + ' | ' + wrapper + ')'
                else:
                    fpath = '<(gunzip -c ' + fpath + ' | ' + wrapper + ')'
                    cmp_fpath = '<(gunzip -c ' + cmp_fpath + ' | ' + wrapper + ')'
            elif fpath.endswith('.gz'):
                fpath = '<(gunzip -c ' + fpath + ')'
                cmp_fpath = '<(gunzip -c ' + cmp_fpath + ')'
            cmdl += ' ' + fpath + ' ' + cmp_fpath
            ret_code = call(cmdl, suppress_output=not is_travis())
            assert ret_code == 0, 'diff returned non-zero: ' + fpath

    @staticmethod
    def _check_dir_not_empty(dirpath, description=None):
        assert verify_dir(dirpath, description=description), dirpath
        contents = [join(dirpath, fname) for fname in os.listdir(dirpath)
                    if not fname.startswith('.')]
        assert len(contents) >= 1, dirpath + ': ' + str(contents)
        assert all(verify_file(realpath(fpath), is_critical=True)
                   for fpath in contents
                   if isfile(realpath(fpath))), dirpath + ': ' + str(contents)
