#!/usr/bin/env python
# __init__.py
# Maintained by P.Legrand
# 14th April 2005
#

#import psyco
#psyco.full()

import os, sys
import  __init__

def _print_autopar(parDict):
    K = parDict.keys()
    K.sort()
    for k in K:
        print "%s = %s" % (k, parDict[k])

def _test1(filename):
    #dc = __init__.Collect(filename)
    try:
        im = __init__.Image(filename)
    except __init__.XIOError:
        print __init__.Image(filename)
        im =  __init__.Image(filename, doInterpret=False)
        if im.type == 'unknown':
            print "\nUnrecognise file type for image %s.\nSorry." % filename
        else:
            print "\nNo plugin to interpret '%s' image type.\nSorry." % im.type
        sys.exit()
    im.info()
    im.data_info()
    
def _export(filename, format):
    try:
        datacoll = __init__.Collect(filename)
        datacoll.intepretImage()
        datacoll.lookup_imageRanges()
        
    except __init__.XIOError:
        print "\nError while trying to acceess %s.\nSorry." % filename
        sys.exit()
    
    newPar = datacoll.export(format)   
    if format == "adp":
        newPar["anomalous"] = "off"
        newPar["spg"] = 0
        
    elif format == "xds":
        newPar["DELPHI"] = newPar["OSCILLATION_RANGE"]*3
        
    _print_autopar(newPar)

if __name__ == "__main__":

    if len(sys.argv) > 1:
        if "-xds" in sys.argv:
            sys.argv.remove("-xds")
            filename = sys.argv[1]
            _export(filename, "xds")
        if "-adp" in sys.argv:
            sys.argv.remove("-adp")
            filename = sys.argv[1]
            _export(filename, "adp")
        else:
            filename = sys.argv[1]
            _test1(filename)
                       
