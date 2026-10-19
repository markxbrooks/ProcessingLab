#!/usr/bin/env python

__version__ = "0.0.1"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "06-12-2004"
__copyright__ = "Copyright (c) 2004  Pierre Legrand"
__license__ = "LGPL"

import os, sys
from xupy import XParam, xdsInp2Param, opReadCl, opWriteCl
#import TaskManager
import popen2
import types

XDSHOME = '/home/legrand/prog/pxsoft/XDS/XDS-linux_ifc_Intel'

class XDS(XParam):

    def __init__(self, obj=None):
        """Constructor for the Param classes from file or string."""
        self.__listeners = []
        self.__cancelled = 0
        self.__execfile = os.path.join(XDSHOME,"xds_par")
        if type(obj) == types.FileType:
            exec obj.read() in self.__dict__
            obj.close()
        if type(obj) == type(""):
            exec obj in self.__dict__
     
    def addListener (self, listener):
        ll = self.__listeners
        if listener not in ll:
            ll.append(listener)

    def addListeners (self, listeners):
        for l in listeners:
            self.addListener(l)

    def removeListener (self, listener):
        ll = self.__listeners
        ll.remove(listener)

    def removeListeners (self, listeners):
        for l in listeners:
            self.removeListener(l)

    def fireListeners (self, value):
        for func in self.__listeners:
            func(value)

    def cancel (self):
        self.__cancelled = 1

    def run(self, runXDS_args):
        rdir, rsave = runXDS_args
        self.__cancelled = 0
        self.fireListeners(0)
        
        In = xdsInp2Param(inp_str=opReadCl("XDS.INP"))
        self.mix(In)
        
        if 0: #rdir:
            if not os.path.exists(rdir):
                try: os.mkdir(rdir)
                except:
                    print "STOP! Can't creat xds working directory:",rdir
                    sys.exit()
            if os.path.isdir(rdir):
                os.chdir(rdir)
        
        opWriteCl("XDS.INP", "%s" % self)
        xdsProcess = popen2.Popen3(self.__execfile)
        while 1:
            if xdsProcess.poll() != -1:
                break
            if self.__cancelled:
                os.kill(xdsProcess.pid, 9)
                break
            nl = xdsProcess.fromchild.readline()
            if nl: print nl,

        if rsave: saveLastVersion(LP_names)
        #return "XDS finished for integration of images %4d to %4d" % \
        #     (self.DATA_RANGE[0],self.DATA_RANGE[1])
        os.chdir("../")
        return "toto"
        # use iteration to calculate -- may take quite a while
        while 1:
            last = 0
            for x in range(1, value + 1):
                if self.__cancelled: break
                result = result * x
                next = x * 100 / value
                if next != last:
                    self.fireListeners(next)
                    last = next
        self.fireListeners(100)
        if self.__cancelled: result = -1
        return result


if __name__ == "__main__":
    #In = xdsInp2Param(inp_str=opReadCl("XDS.INP"))
    newrun = XDS()
    #newrun.mix(In)
    #print newrun
    newrun.run(("xds_run1",None))
