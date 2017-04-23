"""Centralize running of external commands, providing logging and tracking.
"""
import collections
import os
import subprocess
from os.path import isfile

from ngs_utils.file_utils import file_transaction, verify_file
from ngs_utils.logger import info, err


def run(cmd, output_fpath=None, input_fpath=None, checks=None, stdout_to_outputfile=True,
        stdout_tx=True, reuse=False, stderr_fpath=None, env_vars=None):
    """Run the provided command, logging details and checking for errors.
    """
    if output_fpath and reuse:
        if verify_file(output_fpath, silent=True):
            info(output_fpath + ' exists, reusing')
            return output_fpath
        if not output_fpath.endswith('.gz') and verify_file(output_fpath + '.gz', silent=True):
            info(output_fpath + '.gz exists, reusing')
            return output_fpath

    env = os.environ.copy()
    if env_vars:
        for k, v in env_vars.items():
            if v is None:
                if k in env:
                    del env[k]
            else:
                env[k] = v

    if checks is None:
        checks = [file_nonempty_check]

    def _try_run(_cmd, _output_fpath, _input_fpath, _stderr_fpath):
        try:
            info(' '.join(str(x) for x in _cmd) if not isinstance(_cmd, basestring) else _cmd)
            _do_run(_cmd, checks, env, _output_fpath, _input_fpath, _stderr_fpath)
        except:
            raise

    if output_fpath:
        if isfile(output_fpath):
            os.remove(output_fpath)
    if output_fpath:
        if stdout_tx:
            with file_transaction(None, output_fpath) as tx_out_file:
                if stdout_to_outputfile:
                    cmd += ' > ' + tx_out_file
                else:
                    cmd += '\n'
                    cmd = cmd.replace(' ' + output_fpath + ' ', ' ' + tx_out_file + ' ') \
                             .replace(' "' + output_fpath + '" ', ' ' + tx_out_file + '" ') \
                             .replace(' \'' + output_fpath + '\' ', ' ' + tx_out_file + '\' ') \
                             .replace(' ' + output_fpath + '\n', ' ' + tx_out_file) \
                             .replace(' "' + output_fpath + '"\n', ' ' + tx_out_file + '"') \
                             .replace(' \'' + output_fpath + '\'\n', ' ' + tx_out_file + '\'') \
                             .replace('\n', '')
                _try_run(cmd, tx_out_file, input_fpath, stderr_fpath)
        else:
            _try_run(cmd, output_fpath, input_fpath, stderr_fpath)

    else:
        _try_run(cmd, None, input_fpath, stderr_fpath)


def find_bash():
    for test_bash in [find_cmd("bash"), "/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")


def find_cmd(cmd):
    try:
        return subprocess.check_output(["which", cmd]).strip()
    except subprocess.CalledProcessError:
        return None


def _normalize_cmd_args(cmd):
    """Normalize subprocess arguments to handle list commands, string and pipes.
    Piped commands set pipefail and require use of bash to help with debugging
    intermediate errors.
    """
    if isinstance(cmd, basestring):
        # check for standard or anonymous named pipes
        if cmd.find(" | ") > 0 or cmd.find(">(") or cmd.find("<("):
            return "set -o pipefail; " + cmd, True, find_bash()
        else:
            return cmd, True, None
    else:
        return [str(x) for x in cmd], False, None


def _do_run(cmd, checks, env=None, output_fpath=None, input_fpath=None, _stderr_fpath=None):
    """Perform running and check results, raising errors for issues.
    """
    cmd, shell_arg, executable_arg = _normalize_cmd_args(cmd)
    s = subprocess.Popen(cmd, shell=shell_arg, executable=executable_arg,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, close_fds=True, env=env)
    debug_stdout = collections.deque(maxlen=100)
    while 1:
        line = s.stdout.readline()
        if line:
            debug_stdout.append(line)
            info('  ' + line.rstrip())
        exitcode = s.poll()
        if exitcode is not None:
            for line in s.stdout:
                debug_stdout.append(line)
            if exitcode is not None and exitcode != 0:
                error_msg = " ".join(cmd) if not isinstance(cmd, basestring) else cmd
                error_msg += "\n"
                error_msg += "".join(debug_stdout)
                s.communicate()
                s.stdout.close()
                raise subprocess.CalledProcessError(exitcode, error_msg)
            else:
                break
    s.communicate()
    s.stdout.close()
    # Check for problems not identified by shell return codes
    if checks:
        for check in checks:
            if not check(output_fpath, input_fpath):
                raise IOError("External command failed")
        # except subprocess.CalledProcessError as e:
        #     e.returncode
        #     e.cmd


def file_nonempty_check(output_fpath=None, input_fpath=None):
    if output_fpath is None:
        return True
    ok = file_exists_check(output_fpath)
    if not ok:
        err("Did not find non-empty output file {0}".format(output_fpath))
    return ok


def file_exists_check(output_fpath=None, input_fpath=None):
    if output_fpath is None:
        return True
    ok = os.path.exists(output_fpath)
    if not ok:
        err("Did not find output file {0}".format(output_fpath))
    return ok


def file_reasonable_size(output_fpath, input_fpath):
    ok = file_exists_check(output_fpath)
    if not ok:
        return ok
    # named pipes -- we can't calculate size
    if input_fpath.strip().startswith("<("):
        return True
    if input_fpath.endswith((".bam", ".gz")):
        scale = 7.0
    else:
        scale = 10.0
    orig_size = os.path.getsize(input_fpath) / pow(1024.0, 3)
    out_size = os.path.getsize(output_fpath) / pow(1024.0, 3)
    if out_size < (orig_size / scale):
        err("Output file unexpectedly small. %.1fGb for output versus "
            "%.1fGb for the input file. This often indicates a truncated "
            "BAM file or memory errors during the run." % (out_size, orig_size))
        return False
    else:
        return True
