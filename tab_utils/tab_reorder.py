#!/usr/bin/env python
'''
Re-orders columns in a tab delimited file

(Reorders based upon column index, not name)
'''

import sys
import os

from tab_utils.support import gzip_opener

def tab_reorder(fname, column_order, delim='\t'):
    f = gzip_opener(fname).open()
    for line in f:
        if line[0] == '#':
            sys.stdout.write(line)
            continue

        cols = line.rstrip('\n').split(delim)

        cols_used = set()
        outcols = []

        for col_idx in column_order:
            if col_idx == '*':
                for i, val in enumerate(cols):
                    if not i in cols_used:
                        outcols.append(val)
                break
            else:
                outcols.append(cols[col_idx])
                if col_idx < 0:
                    cols_used.add(len(cols) + col_idx)
                else:
                    cols_used.add(col_idx)


        sys.stdout.write('%s\n' % '\t'.join(outcols))
    f.close() 
    
def usage(msg=""):
    if msg:
        print(msg)
    print(__doc__)
    print("""Usage: %s {opts} filename.tab col1,col2,col3

Columns should be specified by index:
    1,5,3 (etc)

Alternatively, ranges may be specified using from:to notation:
    1:3

Also, negative indexes may be used to indicate columns from the right side
    -3:-1

Finally, the catch-all '*' can be used to output all columns which haven't been included yet.

Note: columns start at 0

Options:
    -d delim    Use this (opposed to a tab) for the delimiter

""" % os.path.basename(sys.argv[0]))
    sys.exit(1)
    
def main(argv):
    fname = None
    delim = '\t'
    last = None

    column_order = []

    for arg in argv:
        if arg in ['-h','--help']:
            usage()
        elif last == '-d':
            delim = arg
            last = None
        elif arg in ['-d']:
            last = arg
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        else:
            for val in arg.split(','):
                if ':' in val:
                    l,r = [int(x) for x in val.split(':')]
                    for i in range(l,r+1):
                        column_order.append(i)
                elif val != '*':
                    column_order.append(int(val))
                else:
                    column_order.append(val)
    
    if not fname or not column_order:
        usage("Missing input file or columns!")

    tab_reorder(fname, column_order, delim)
    
if __name__ == '__main__':
    main(sys.argv[1:])

