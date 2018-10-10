import sys
import os
import gzip

def open_gzipsafe(f, mode='r'):
    if f.endswith('.gz') or f.endswith('.gz.tx') or f.endswith('.gzip.tx'):
        try:
            h = gzip.open(f, mode=mode + 't')
        except IOError as e:
            # print('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
            return open(f, mode=mode)
        else:
            if 'w' in mode:
                return h
            else:
                try:
                    h.read(1)
                except IOError as e:
                    # print('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
                    h.close()
                    return open(f, mode=mode)
                else:
                    h.close()
                    h = gzip.open(f, mode=mode + 't')
                    return h
    else:
        return open(f, mode=mode)

class gzip_opener:
    """
    A Python 2.6 class to handle 'with' opening of text files that may
    or may not be gzip compressed.
    """
    def __init__(self, fname):
        self.fname = fname

    def open(self):
        return self.__enter__()

    def __enter__(self):
        if self.fname == '-':
            self.f = sys.stdin
        else:
            self.f = open_gzipsafe(os.path.expanduser(self.fname))
        return self.f

    def __exit__(self, type, value, traceback):
        if self.f != sys.stdin:
            self.f.close()
        return False

def filenames_to_uniq(names,new_delim='.'):
    '''
    Given a set of file names, produce a list of names consisting of the
    uniq parts of the names. This works from the end of the name.  Chunks of
    the name are split on '.' and '-'.
    
    For example:
        A.foo.bar.txt
        B.foo.bar.txt
        returns: ['A','B']
    
        AA.BB.foo.txt
        CC.foo.txt
        returns: ['AA.BB','CC']
    
    '''
    name_words = []
    maxlen = 0
    for name in names:
        name_words.append(name.replace('.',' ').replace('-',' ').strip().split())
        name_words[-1].reverse()
        if len(name_words[-1]) > maxlen:
            maxlen = len(name_words[-1])

    common = [False,] * maxlen
    for i in range(maxlen):
        last = None
        same = True
        for nameword in name_words:
            if i >= len(nameword):
                same = False
                break
            if not last:
                last = nameword[i]
            elif nameword[i] != last:
                same = False
                break
        common[i] = same

    newnames = []
    for nameword in name_words:
        nn = []
        for (i, val) in enumerate(common):
            if not val and i < len(nameword):
                nn.append(nameword[i])
        nn.reverse()
        newnames.append(new_delim.join(nn))
        
    return newnames
