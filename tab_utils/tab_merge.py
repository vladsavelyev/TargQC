#!/usr/bin/env python
'''
Merges tab-delimited files together, combining common columns
'''
import sys
import os
import gzip
from tab_utils.support import filenames_to_uniq


class MergeException(Exception):
    pass
    
def merge_files(fnames,common_cols,uncommon_cols, keycols, noheader=False,collate=True,headercomment=False,keydesc=False,nomissing=False, given_names=None):
    if given_names:
        names = given_names
    else:
        names = filenames_to_uniq([os.path.basename(x) for x in fnames])
    files = []
    file_col_count = 0
    for fname in fnames:
        if fname[-3:] == '.gz':
            files.append(gzip.open(fname))
        else:
            files.append(open(fname))

    commented_lines = {}
    for name,f in zip(names,files):
        for line in f:
            if line[0] == '#':
                if not name in commented_lines:
                    commented_lines[name] = []
                    
                commented_lines[name].append(line[1:])
            else:
                break
        f.seek(0)
    
    if commented_lines:
        commentsout = False
        for name in commented_lines:
            for line in commented_lines[name][:-1]:
                commentsout = True
                sys.stdout.write('# %s %s' % (name,line))
            if not headercomment:
                commentsout = True
                sys.stdout.write('# %s %s' % (name,commented_lines[name][-1]))
                
        if commentsout:
            sys.stdout.write('\n')
    
    headers = True

    # header is the last commented line
    if headercomment:
        header_cols = []
        for name in commented_lines:
            cols = commented_lines[name][-1].rstrip().split('\t')
            for i in common_cols:
                header_cols.append(cols[i])
            if not collate:
                for j in uncommon_cols:
                    for name in names:
                        header_cols.append('%s %s' % (name,cols[j]))
            else:
                for name in names:
                    for j in uncommon_cols:
                        header_cols.append('%s %s' % (name,cols[j]))
            break
        sys.stdout.write('\t'.join(header_cols))
        sys.stdout.write('\n')
        headers = False
    
    #no header... just show column #
    
    elif noheader:
        header_cols = []
        for common in common_cols[1:]:
            header_cols.append('')
            sys.stdout.write('\t')
            
        if len(uncommon_cols) > 1:
            if not collate:
                for j in uncommon_cols:
                    for name in names:
                        header_cols.append('%s col %s' % (name,j))
            else:
                for name in names:
                    for j in uncommon_cols:
                        header_cols.append('%s col %s' % (name,j))
                
        else:
            for name in names:
                header_cols.append('%s' % (name))

        sys.stdout.write('\t'.join(header_cols))
        sys.stdout.write('\n')
        headers = False
        
    lines = []
    for f in files:
        lines.append(None)
        
    while True:
        missing_values = False
        
        for i,f in enumerate(files):
            try:
                while not lines[i]:
                    lines[i] = f.next()
                    
                    if lines[i][0] == '#' or lines[i].strip() == '':
                        lines[i] = None
                    
            except Exception as e:
                lines[i] = None
                pass
        
        good = False
        for line in lines:
            if line:
                good = True
        if not good:
            break
        
        outcols = []
        values = []
        common_keys = []
        num_of_columns = len(common_cols) + len(uncommon_cols)
        
        # look for missing values
        for i,line in enumerate(lines):
            if not line:
                missing_values = True
                values.append(['',] * num_of_columns)
                continue
                
            cols = line.rstrip().split('\t')
            if not cols:
                missing_values = True
                values.append(['',] * num_of_columns)
            
            while len(cols) < file_col_count:
                cols.append('')
                
            values.append(cols)
            
            keys = []
            for j,num in zip(keycols[0],keycols[1]):
                if num and not headers:
                    keys.append(float(cols[j]))
                else:
                    keys.append(cols[j])
            common_keys.append((keys,i))
            
        common_keys.sort()
        if keydesc:
            common_keys.reverse()
            
        for keys,i in common_keys:
            if keys == common_keys[0][0]:
                lines[i] = None
                if not outcols:
                    for j in common_cols:
                        outcols.append(values[i][j])
            else:
                missing_values = True
                values[i] = ['',] * num_of_columns

        if not values:
            continue

        if nomissing and missing_values:
            continue

        # first line is header
        if headers:
            headers = False
            if not file_col_count:
                file_col_count = len(values[0])
            if not collate:
                for j in uncommon_cols:
                    for name in names:
                        outcols.append('%s %s' % (name,values[0][j]))
            else:
                for name in names:
                    for j in uncommon_cols:
                        outcols.append('%s %s' % (name,values[0][j]))
        else:
            if not collate:
                for j in uncommon_cols:
                    for cols in values:
                        if len(cols)>j:
                            outcols.append(cols[j])
                        else:
                            outcols.append('')

            else:
                for cols in values:
                    for j in uncommon_cols:
                        if len(cols)>j:
                            outcols.append(cols[j])
                        else:
                            outcols.append('')
        if outcols:
            sys.stdout.write('%s\n' % '\t'.join(outcols))

    
    for f in files:
        f.close()
        
def usage():
    print(__doc__)
    print("""\
Usage: %s {opts} common_cols merge_cols files

common_cols and merge_cols should be a comma-separated list of column numbers.

Files must be in the same sort order (given as  if not ascending text).
If the case of extra rows in one of the files, blank values will be 
substituted. New column names will be guessed based upon the filenames. 
Commented lines and blank lines are ignored, except for any commented lines 
that are at the begining of the file - these are kept.  This assumes that the 
first non-comment, non-blank row is the header.

Options:
    -headercomment     the header is last commented line ('#')
    -noheader          the files have no header row
    -nomissing         discard rows with missing values
    -collate           order uncommon values by file first, not column
    
    -keycols col,col   if there are missing values, use these columns to 
                       determine which file has missing data.  If the col ends 
                       in 'n', this is taken to be a number.
                       (defaults to common cols, in order)
                       
    -keydesc           keys are sorted descending order
                       (defaults to ascending)

Valid column definitions:
    1,2,3,4
    1-4
    1-4,5

""" % os.path.basename(sys.argv[0]))

def _split_cols(arg):
    ret=[]
    nums=[]
    for x in arg.split(','):
        num = False
        if x[-1] == 'n':
            num = True
            x = x[:-1]
        
        if '-' in x:
            s,e = x.split('-')
            vals = list(range(int(s)-1,int(e)))
            ret.extend(vals)
            for v in vals:
                nums.append(num)
        else:
            ret.append(int(x)-1)
            nums.append(num)
    return ret,nums

def main(argv):
    noheader=False
    collate = False
    common = None
    uncommon = None
    keycols = None
    keydesc = False
    headercomment = False
    nomissing = False
    files=[]
    given_names=[]
    
    last = None
    
    for arg in argv:
        if arg == '-h':
            usage()
            sys.exit(1)
        if last == '-keycols':
            keycols = _split_cols(arg)
            last = None
        elif last == '-names':
            given_names = [x.strip() for x in arg.split(',')]
            last = None
        elif arg in ['-keycols', '-names']:
            last = arg
        elif arg == '-keydesc':
            keydesc = True
        elif arg == '-nomissing':
            nomissing = True
        elif arg == '-headercomment':
            headercomment = True
        elif arg == '-noheader':
            noheader = True
        elif arg == '-collate':
            collate = True
        elif not common:
            common = _split_cols(arg)
        elif not uncommon:
            uncommon = _split_cols(arg)
        elif os.path.exists(arg):
            files.append(arg)
        else:
            usage()
            sys.exit(1)
    if not files or not common or not uncommon:
        usage()
        sys.exit(1)
    
    if not keycols:
        keycols = common
    
    merge_files(files,common[0],uncommon[0],keycols,noheader,collate,headercomment,keydesc,nomissing, given_names)

if __name__ == '__main__':
    main(sys.argv[1:])
