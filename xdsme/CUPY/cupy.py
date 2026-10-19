#! /usr/bin/env python
# 
# CuPy is distributed under GNU General Public License, 
# available at http://www.gnu.org/copyleft/gpl.html
#
# Credits: Cupy uses YAPTU code by Alex Martelli.
# See http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52305
#
# TODO:
#   - get informations from cns pdb header
#   - when writing files, check if file allready exist...
#   - rethink file handling
#   - rethink file naming
#   - better handling of errors (add try, Exept...)
#   - test ... [cns example files, with another pdb, hkl]
#   - add coments to code
#   - translate all french in english
#       and check spelling !
#   - better handling of argument
#      - add a -i interactive filling...
#      - add -o option for file redirection 
#   - speed improvement, add caching system for template ?

import sys, re, string, types, os

__version__ = "0.4.0"
__author__ = "Pierre Legrand (pierre.legrand@crchul.ulaval.ca)"
__date__ = "24 Aout 2002"
__copyright__ = "Copyright (c) 2001 Pierre Legrand"
__license__ = "GPL"


_usage = 'Usage:  %s cns_input.inp [param.dat ["bmin = 5"]]' % (sys.argv[0])
_help = """
 CuPy is a versatil CNS templating Utility for/with Python

 TCPuPy is a simple hack to ease the scripting and
 templating of CNS using the python language.
 Cupy needs a python distribution (version 1.52 or higher) to work.
 CNS executable file is traced using the CNS_SOLVE environement varaible.

 Usage 1:  tcpupy.py  cns_input.inp
 
    Process cns_input.inp file and output 2 files:
    cns_input_def.par: contaning the defaults parameters as
                       define in cns_input.inp, and
    cns_input.tpl:     a cns_input template file.
      

 Usage 2:  tcpupy.py  cns_input.inp myparam.dat ["high_res = 1.9; bmin = 5"]
 
    Run cns_input.inp using defaults parameters exepte if they
    are define in "par.py".
    Parameters can also be input at the command line. If a third argument
    is given it will be evaluated and the new parameter(s) value(s) will
    resplace those defined before. Serveral parameters can be difined this
    way using the python synthax:
            "high_res = 1.9; bmin = 5"
        or  "high_res, bmin = 1.9, 5"
  ##  Output cns log in "output_file" or if not specified in "cns_input.out"
    

 Usage 3:  From inside python. For example:
 
     from tcpupy import run_cns
     class MyPar:
         a, b, c = 61.76, 40.73, 26.74
         alpha, beta, gamma = 90.0, 90.0, 90.0
         sg = "P2(1)2(1)2(1)"
         reflection_infile_1 = "hkl.cv"
         low_res, high_res = 20.0, 4.0
         structure_infile = "generate.mtf"
         parameter_infile_1 = "CNS_TOPPAR:protein_rep.param"
         parameter_infile_2 = "toppar/ligand.par"
         parameter_infile_3 = "CNS_TOPPAR:water_rep.param"
         coordinate_infile = "generate.pdb"
         coordinate_outfile = "rigid_1.pdb"
         r_group_1 = "segid A"
         r_group_2 = "segid B"
         
     run_cns("rigid.inp", MyPar)
     
     MyPar.coordinate_infile = "rigid_1.pdb"
     MyPar.coordinate_outfile = "rigid_2.pdb"
     MyPar.high_res = 2.0
     
     run_cns("rigid.inp", MyPar, "rigid_2.out")
     
     MyPar.coordinate_infile = "rigid_2.pdb"
     MyPar.coordinate_outfile = "minimize.pdb"
     MyPar.minimize_nstep = 50
     
     run_cns("minimize.inp", MyPar)
     
     etc... """

# set logical variables
true  = yes = on  = "true"
false = no  = off = "false"

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def opWriteCl(filename, _str):
    f = open(filename,"w")
    r = f.write(_str)
    f.close()

def Walk(root, recurse=0, pattern='*', return_folders=0):
    import fnmatch, os, string
    # initialize
    result = []
    # must have at least root folder
    try:
        names = os.listdir(root)
    except os.error:
        print ">>> Warning, can't find %s directory" % (root)
        return result
    # expand pattern
    pattern = pattern or '*'
    pat_list = string.splitfields(pattern, ';')
    # check each file
    for name in names:
        fullname = os.path.normpath(os.path.join(root, name))
        # grab if it matches our pattern and entry type
        for pat in pat_list:
            if fnmatch.fnmatch(name, pat):
                if os.path.isfile(fullname) or (return_folders and os.path.isdir(fullname)):
                    result.append(fullname)
                continue		    
        # recursively scan other folders, appending results
        if recurse:
            if os.path.isdir(fullname) and not os.path.islink(fullname):
                result = result + Walk(fullname, recurse, pattern, return_folders)    
    return result

def fix_pdb(pdbin, changes = []):
    """Verify if a trailling END is present in pdbin, and add on if not.
    As a second argument you can add a list of string couple to make
    replacement in the pdb file. Exemple:
    fix_pdb("toto.pdb", [["OH2 HOH Z","OH2 WAT Z"],["CD  LEU","CD1 LEU"]])
    In case of changes original pdbin file is save with the ".bck"
    extention.
    """
    inf = open(pdbin).read()
    ind = string.find(string.rstrip(inf[-100:]),"END")
    if ind == -1 or changes:
        opWriteCl(pdbin+".bck",inf)
        if changes:
            for change in changes:
                if len(change) == 2:
                    inf = string.replace(inf,change[0],change[1])
        add = ""
        if ind == -1:
            inf = string.rstrip(inf)
            if inf[-1:] != "\n": add = "\nEND\n"
            else: add = "END\n"
        opWriteCl(pdbin,inf + add)
    
# utility stuff to avoid tests in the mainline code
class _nevermatch:
    "Polymorphic with a regex that never matches"
    def match(self, line):
        return None
_never = _nevermatch()     # one reusable instance of it suffices
def identity(string, why):
    "A do-nothing-special-to-the-input, just-return-it function"
    return string
def nohandle(string):
    "A do-nothing handler that just re-raises the exception"
    raise

# and now the real thing
class copier:
    "Smart-copier (YAPTU) class"
    def copyblock(self, i=0, last=None):
        "Main copy method: process lines [i,last) of block"
        def repl(match, self=self):
            "return the eval of a found expression, for replacement"
            # uncomment for debug: print '!!! replacing',match.group(1)
            expr = self.preproc(match.group(1), 'eval')
            try: return str(eval(expr, self.globals, self.locals))
            except: return str(self.handle(expr))
        block = self.locals['_bl']
        if last is None: last = len(block)
        while i<last:
            line = block[i]
            match = self.restat.match(line)
            if match:   # a statement starts "here" (at line block[i])
                # i is the last line to _not_ process
                stat = match.string[match.end(0):].strip()
                j=i+1   # look for 'finish' from here onwards
                nest=1  # count nesting levels of statements
                while j<last:
                    line = block[j]
                    # first look for nested statements or 'finish' lines
                    if self.restend.match(line):    # found a statement-end
                        nest = nest - 1     # update (decrease) nesting
                        if nest==0: break   # j is first line to _not_ process
                    elif self.restat.match(line):   # found a nested statement
                        nest = nest + 1     # update (increase) nesting
                    elif nest==1:   # look for continuation only at this nesting
                        match = self.recont.match(line)
                        if match:                   # found a contin.-statement
                            nestat = match.string[match.end(0):].strip()
                            stat = '%s _cb(%s,%s)\n%s' % (stat,i+1,j,nestat)
                            i=j     # again, i is the last line to _not_ process
                    j=j+1
                stat = self.preproc(stat, 'exec')
                stat = '%s _cb(%s,%s)' % (stat,i+1,j)
                # for debugging, uncomment...: print "-> Executing: {"+stat+"}"
                exec stat in self.globals,self.locals
                i=j+1
            else:       # normal line, just copy with substitution
                self.ouf.write(self.regex.sub(repl,line))
                i=i+1
    def __init__(self, regex=_never, dict={},
            restat=_never, restend=_never, recont=_never, 
            preproc=identity, handle=nohandle, ouf=sys.stdout):
        "Initialize self's attributes"
        self.regex   = regex
        self.globals = dict
        self.locals  = { '_cb':self.copyblock }
        self.restat  = restat
        self.restend = restend
        self.recont  = recont
        self.preproc = preproc
        self.handle  = handle
        self.ouf     = ouf
    def copy(self, block=None, inf=sys.stdin):
        "Entry point: copy-with-processing a file, or a block of lines"
        if block is None: block = inf.readlines()
        self.locals['_bl'] = block
        self.copyblock()

def inp2tpl(inp_f):
    """ Seperate a CNS input file in a template and a paramter file."""
    assign = re.compile("[\s\};][\w|\.]+=.+;")
    number = re.compile("[-|+]*[\d]+\.?\d*")
    n_par, tpl_s, par_s = 1, "", ""
    while 1:
        line = inp_f.readline()
        if not line: break
        if string.count(line, "define("):
            tpl_s = tpl_s + line
            while 1:
                line = inp_f.readline()
                if not line: break
                n_par = ( n_par + string.count(line, "(")
                                - string.count(line, ")"))
                if not n_par: break
                m = assign.search(line)
                if  m:
                    substr = line[m.start():m.end()].split(";")
                    for subsubstr in substr[:-1]:
                        m_eq = string.index(subsubstr, "=")
                        var_name = string.strip(subsubstr[:m_eq])
                        # Translate dot containing names (canot be python variable names)
                        if var_name.count("."):
                           var_name_corr = var_name.replace(".","_d_")
                        else: var_name_corr = var_name
                        var_value = string.strip(subsubstr[m_eq+1:])
                        fmt = " %s=\"@%s@\""
                        var_value_lim = var_value[0]+var_value[-1]
                        if number.match(var_value):
                            fmt = " %s=@%s@"
                            par_s = par_s + "%s = %s\n" % (var_name_corr, var_value)
                        elif var_value_lim == "\"\"":
                            par_s = par_s + "%s = %s\n" % (var_name_corr, var_value)
                        else:
                            par_s = par_s + "%s = \"%s\"\n" % (var_name_corr, var_value)
                        if ( var_value_lim == "()" or 
                            var_value == "true" or 
                            var_value == "false" ): fmt = " %s=@%s@"
                        line = string.replace(
                            line, subsubstr, fmt % (var_name, var_name_corr))
                tpl_s = tpl_s + line
        tpl_s = tpl_s + line
    inp_f.close()
    return tpl_s, par_s

def find_inpfn(inp_n):
    try:
        f = open(inp_n)
	f.close()
	return inp_n
    except IOError:
        root = os.path.join(os.environ["CNS_SOLVE"], "inputs")
        files = Walk(root, 1, '*.inp;*.INP')
        #print '>>  File %s not found in current directory.' % (inp_n)
	#print '>>? Trying to find it in location: %s' % (root)
        for file in files:
               if string.count(file,"/"+inp_n):
	           #print ">>! Matching input file : %s" % (file)
	           return file
	print ">># Error ! Can't find file %s." % (inp_n)
        sys.exit()

def tpl2inp(template_s, param, inp):
    """Fill the template with param and write the resulting cns input file"""
    import re
    rex=re.compile('@([^@]+)@')
    t = open("CuPy.tpl", "w")
    t.write(template_s)
    t.close()
    i = open(inp, "w")
    cop = copier(rex, vars(param), ouf=i)
    cop.copy(inf=open("CuPy.tpl"))
    i.close()

def exec_cns(inp_n, out_n=None, exec_name = None):
    if not out_n : out_f = " "
    else: out_f = " > " + out_n
    if not exec_name: exec_name = "cns"
    #cns_bin = os.path.join(os.environ["CNS_SOLVE"], os.environ["CNS_ARCH"], "bin","cns")
    #os.system(cns_bin+" < "+inp_n + out_f) # use popen instead ?
    os.system(exec_name +" < "+inp_n + out_f) # use popen instead ?

class Param:
    """A simple class to handle the parameters that fills templates"""
    def __init__(self, obj=None):
        """Constructor for the Param classes from file or string."""
        if type(obj) == types.FileType:
            exec obj.read() in self.__dict__
        if type(obj) == types.StringType:
            exec obj in self.__dict__
    def mix(self, custom):
        """Complete and/or replace the current param with custom params."""
        self.__dict__.update(custom.__dict__)
            
def run_cns(inp_fn, my_par, out_f = None, exec_name = "cnsi"):
    new_inp = "CuPy_"+inp_fn
    inp_fn = find_inpfn(inp_fn)
    tpl_s, par_s = inp2tpl(open(inp_fn))
    def_par = Param(par_s)
    def_par.mix(my_par)
    tpl2inp(tpl_s, def_par, new_inp)
    
    exec_cns(new_inp, out_f, exec_name)

if __name__=='__main__':
    
    if len(sys.argv) == 1: print _usage; sys.exit()
    if sys.argv[1] == '-help' or \
       sys.argv[1] == '--help' or \
       sys.argv[1] == '-h': print _help; sys.exit()

    if sys.argv[1] == 'test':
        root = os.path.join(os.environ["CNS_SOLVE"], "inputs")
        files = Walk(root, 1, '*.inp;*.INP')
        for file in files:
               print file
        #inp_f = open("bindividual.inp")
        #tpl, par = inp2tpl(inp_f)
        #P = CnsParam(par)
        #print " P >> ", P.a, P.b, P.c
        #print " P >> ", P.anom_library, P.sg
        sys.exit()
        
    if len(sys.argv) == 2:
        file_name = find_inpfn(sys.argv[1])
        base_name = os.path.splitext(os.path.split(file_name)[1])[0]
        tpl, par = inp2tpl(open(file_name))
        #open(base_name+".tpl","w").write(tpl)
        open(base_name+"_par.py","w").write(par)
        print par
    
    if len(sys.argv) >= 3:
        file_name = sys.argv[1]
        base_name = os.path.splitext(os.path.split(file_name)[1])[0]
        param_name = sys.argv[2]
        my_par = Param(open(param_name))
        if len(sys.argv) >= 4: my_par.mix(Param(sys.argv[3]))
        run_cns(file_name, my_par)
