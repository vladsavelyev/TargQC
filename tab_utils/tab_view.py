#!/usr/bin/env python
'''
A data aware tab-delimited file viewer

Works by reading in the first few lines to determine the appropriate widths
for each of the columns.  It will then display the data with the appropriate
spacing to keep columns together.  If a future cell is larger than the
predetermined size, it is truncated.

This can then be fed into something like 'less' for paging
'''

import sys
import os
import math
from tab_utils.support import gzip_opener

def usage():
    print(__doc__)
    print("""Usage: %s {opts} filename.tab

Options:
-l lines    The number of lines to read in to estimate the size of a column.
            [default 100]
-d delim    Use this (opposed to a tab) for the delimiter

-max size   The maximum length of a column (default: unlimited)
-min size   The minimum length of a column (default: 0)

""" % os.path.basename(sys.argv[0]))
    sys.exit(1)

def main(argv):
    fname = '-'
    preview_lines = 100
    delim = '\t'
    max_size = None
    min_size = 0
    last = None
    for arg in argv:
        if arg in ['-h','--help']:
            usage()
        elif last == '-l':
            preview_lines = int(arg)
            last = None
        elif last == '-d':
            delim = arg
            last = None
        elif last == '-min':
            min_size = int(arg)
            last = None
        elif last == '-max':
            max_size = int(arg)
            last = None
        elif arg in ['-l','-d','-max','-min']:
            last = arg
        elif arg == '-':
            fname = '-'
        elif os.path.exists(arg):
            fname = arg

    tab_view(fname, preview_lines, delim, max_size, min_size)

def tab_view(fname, preview_lines, delim, max_size, min_size):
    colsizes = []
    coltypes = []
    preview_buf = []
    prev_count = 0
    inpreview = True

    try:
        f = gzip_opener(fname).open()
        for line in f:
            if inpreview and line[0] == '##':
                preview_buf.append(line)
            else:
                if inpreview:
                    cols = line.rstrip().split(delim)

                    for i, col in enumerate(cols):
                        if len(colsizes) <= i:
                            colsizes.append(len(col))
                            coltypes.append('i')
                        elif len(col) > colsizes[i]:
                            colsizes[i] = len(col)
                        try:
                            v = int(col)
                        except:
                            coltypes[i] = 't'

                    preview_buf.append(line)
                    prev_count += 1
                    if prev_count >= preview_lines:
                        if max_size:
                            colsizes = [ min(max_size, int(math.ceil(x))) for x in colsizes ]
                        else:
                            colsizes = [ max(min_size, int(math.ceil(x))) for x in colsizes ]
                        for preview in preview_buf:
                            _write_cols(preview, colsizes, coltypes)
                        preview_buf = None
                        inpreview = False
                else:
                    _write_cols(line, colsizes, coltypes)
        if f != sys.stdin:
            f.close()
        if preview_buf:
            colsizes = [ int(math.ceil(x)) for x in colsizes ]
            for preview in preview_buf:
                _write_cols(preview, colsizes, coltypes)
    except KeyboardInterrupt:
        print("")
        pass
    except IOError:
        print("")
        pass


def _write_cols(line, colsizes, coltypes):
    cols = line.rstrip().split('\t')

    if len(cols) == 1:
        sys.stdout.write('%s\n' % cols[0])
        return

    while len(cols) < len(colsizes):
        cols.append('')


    for i, col in enumerate(cols):
        if i >= len(colsizes):
            val = '%s' % col  # for headers w/o values
        elif len(col) > colsizes[i]:
            # if too big, show as much as possible, and indicate the
            # truncation with '$'
            val = '%s$' % col[:colsizes[i]-1]
        elif coltypes[i] == 'i':
            # numbers right justified
            val = col.rjust(colsizes[i])
        else:
            # text left justified
            val = col.ljust(colsizes[i])

        sys.stdout.write(val)
        if i < (len(cols) - 1):
            sys.stdout.write('  ')
    sys.stdout.write('\n')

if __name__ == '__main__':
    main(sys.argv[1:])

