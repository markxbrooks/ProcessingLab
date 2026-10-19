#!/usr/bin/env python

__version__ = "0.1.0"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "20-12-2004"
__copyright__ = "Copyright (c) 2004  Pierre Legrand"
__license__ = "LGPL"

from XDS import XDS
from xupy import xdsInp2Param, opReadCl
from qt import * # Generally advertised as safe
from TaskManager import TaskManager
import re
import sys


class XDSgui(QWidget):
    def __init__(self, *args):
        self.__running = 0
        QWidget.__init__(self, *args)
        self.setCaption("XDS logging application")
        self.layout = QGridLayout(self, 6, 2, 5, 10)
        self.tsdisp = QTextEdit(self)
        self.tsdisp.setMinimumSize(600, 500)
        self.tsdisp.setTextFormat(Qt.LogText)
        self.tsdisp.setFont(QFont("Courier",9))
        self.tscount = QLabel("", self)
        self.tscount.setFont(QFont("Sans", 60))
        self.tsstep = QLabel("", self)
        self.tsstep.setFont(QFont("Sans", 20))
        self.tsstep.setMinimumSize(160, 100)
        self.progress = QProgressBar(self)
        self.log = QPushButton("&Start XDS", self)
        self.cancel = QPushButton("&Cancel XDS", self)
        self.quit = QPushButton("&Quit", self)
        self.layout.addMultiCellWidget(self.progress, 5, 5, 0, 1)
        self.layout.addMultiCellWidget(self.tsdisp, 0, 4, 1, 1)
        self.layout.addWidget(self.tscount, 0, 0)
        self.layout.addWidget(self.tsstep, 1, 0)
        self.layout.addWidget(self.log, 2, 0)
        self.layout.addWidget(self.cancel, 3, 0)
        self.layout.addWidget(self.quit, 4, 0)
        self.connect(self.log, SIGNAL("clicked()"),
                     self.start_XDS)
        self.connect(self.cancel, SIGNAL("clicked()"),
                     self.doCancel)
        self.connect(self.quit, SIGNAL("clicked()"),
                     self.finishing)
        # A timer to periodically call periodicCall :-)
        self.timer = QTimer()
        QObject.connect(self.timer,
              SIGNAL("timeout()"),
              self.periodicCall)
        # Start the timer -- this replaces the initial call 
        # to periodicCall
        self.timer.start(200)
        self.tm = TaskManager(1)
        print "Contructed..."
    
    def periodicCall(self):
        """
        Check every 100 ms if there is something new in the queue.
        """
        if self.__running:
            newoutp = self.newrun.getOutp()
            if newoutp:
                self.tsdisp.append(newoutp)
                self.progress.setProgress(100.*self.newrun.step/self.nsteps)
                self.tscount.setText("%d/%d" % (self.newrun.step,self.nsteps))
                self.tsstep.setText(self.newrun.step_name)
            #if self.newrun.running == 0:
            #    self.__running = 0
            #    self.tsdisp.append("\n")
            
    def start_XDS(self):
        Inp = xdsInp2Param(inp_str=opReadCl("XDS.INP"))
        self.newrun = XDS()
        self.newrun.inpParam.mix(Inp)
        
        JOB = self.newrun.inpParam.JOB
        if JOB == "ALL" or JOB == "XDS": self.nsteps = 9
        elif type(JOB) == str: self.nsteps = 2
        else: self.nsteps = len(JOB) +1
        
        self.tm.runTask(self.newrun.run, ("xds_run1",))
        self.__running = 1
        self.tscount.setText("%d/%d" % (0,self.nsteps))
        #self.tsdisp.append(`self.newrun`)

    def doCancel(self):
        print "Trying to cancel..."
        self.__running = 0
        self.newrun.cancel()
    
    def finishing(self):
        self.close()
        if self.__running: self.doCancel()
        self.tm.terminate()
        print "Bye-bye!"

    def __del__(self):
        if self.__running: self.doCancel()
        self.tm.terminate()
        print "Tchaowww!"
                

if __name__ == "__main__":
   app = QApplication(sys.argv)
   app.connect(app, SIGNAL('lastWindowClosed()'), app,
               SLOT('quit()'))
   xdsgui = XDSgui()
   #QObject.connect(xdsgui.quit, SIGNAL("clicked()"), xdsgui.finishing)
   #QObject.disconnect(xdsgui.log, SIGNAL("clicked()"),xdsgui.log_timestamp)
   app.connect(app, SIGNAL('lastWindowClosed()'), xdsgui.__del__)
   app.setMainWidget(xdsgui)
   xdsgui.show()
   sys.exit(app.exec_loop())

