#!/usr/bin/env python
"""
02/04/03 legrand@embl-grenoble.fr
"""
import  sys,os,fnmatch
import  xupy

facteur_repr = [['Compl','compl',"%6.1f"],
                ['ComplL','complL',"%6.1f"],
                ['Unique','unique',"%6d"],
                ['Total','total',"%6d"],
                ['Compar','compar',"%6d"],
                ['Rsym','rsym',"%6.1f"],
                ['RsymL','rsymL',"%6.1f"],
                ['Rmeas','rmeas',"%6.1f"],
                ['I/sig','isig',"%6.1f"],
                ['I/sigL', 'isigL',"%6.1f"],
                ['Total3','total3',"%6d"],
                ['Compl3','compl3',"%6.1f"],
                ['Compa3','compar3',"%6d"],
                ['Rsym3','rsym3',"%6.1f"],
                ['Rmea3','rmeas3',"%6.1f"],
                ['Misfit','misfit',"%6d"],
                ['Absent','absent',"%6d"],
                ['AbsIav','AbsIav',"%6.1f"]]


def list_lp(pattern='*'):
    result = []
    names = os.listdir(".")
    for name in names:
            if fnmatch.fnmatch(name, pattern) and os.path.isfile(name):
                    result.append(name)
    result.sort()
    return result
    
if __name__ == '__main__':
    digit_template = ".???"
    if "-old" in sys.argv:
        digit_template = ".*"
        sys.argv.remove("-old")
    
    facteur_names = []
    for a in facteur_repr: facteur_names.append(a[0])
    
    if len(sys.argv) >= 2:
        file_list = sys.argv[1:]
    else:
        file_list = list_lp("CORRECT.LP"+digit_template)
        
    print "   "+len(facteur_names)*"%7s" % tuple(facteur_names)
    for lp in file_list:
        stat = xupy.resum_scaling(lpf=lp)
        print "%3s" % (lp.split(".")[-1]),
        if not stat: print "Error in CORRECT."; continue
        for k in facteur_repr:
            print k[2] % stat[k[1]],
        print
