#!/usr/bin/env python
'''
Adds a new constant column to each line in a tab-delimited file
'''

import sys
import os

from tab_utils.support import gzip_opener

def tab_tag(fname, colname, colvalue, colidx, delim='\t', noheader=False):
    f = gzip_opener(fname).open()
    header = not noheader
    for line in f:
        if line[0] == '#':
            sys.stdout.write(line)
            continue

        cols = line.rstrip('\n').split(delim)
        if header:
            cols.insert(colidx, colname)
            header = False
        else:
            cols.insert(colidx, colvalue)

        sys.stdout.write('%s\n' % '\t'.join(cols))
 
    f.close()
    
def usage(msg=""):
    if msg:
        print(msg)
    print(__doc__)
    print("""Usage: %s {opts} filename.tab tag

Note: columns start at 1

Options:
    -d delim    Use this (opposed to a tab) for the delimiter
    -name val   The name for this field (in header) (default: tag)
    -pos idx    Insert the new column as column #idx (default: 1)
    -noheader   This file doesn't have a header


""" % os.path.basename(sys.argv[0]))
    sys.exit(1)
    
def main(argv):
    fname = None
    delim = '\t'
    noheader  = False
    name = None
    tag = None
    pos = 0

    last = None

    for arg in argv:
        if arg in ['-h','--help']:
            usage()
        elif last == '-d':
            delim = arg
            last = None
        elif last == '-name':
            name = arg
            last = None
        elif last == '-pos':
            pos = int(arg) - 1
            last = None
        elif arg in ['-d', '-name', '-pos']:
            last = arg
        elif arg == '-noheader':
            noheader = True
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        elif not tag:
            tag = arg
        else:
            usage("Unknown option: %s" % arg)
    
    if not fname or not tag:
        usage("Missing input file or tag!")

    if not name:
        name = tag

    tab_tag(fname, name, tag, pos, delim, noheader)
    
if __name__ == '__main__':
    main(sys.argv[1:])

