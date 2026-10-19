#!/usr/bin/env python

__version__ = "0.4.6"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "16-11-2004"
__copyright__ = "Copyright (c) 2004 Pierre Legrand"
__license__ = "LGPL"


options = ["CCP4","CNS","FALL","SHELX","SOLVE"]
usage   = """
>>>   Usage : %s filein.hkl -[a|n] [free_hkl_to_inherit] MODE \n
        MODE can be on off these:\n
          %s\n
        option:
          -n : normal output (defaulf)
          -a : anomal output\n
>>>   Cell parameters and space group number are taken
      from the input file header.
"""

fmt_inp = """
<== Input file:           %(file_name_in)s
<<< Space group number:      %(spgn_in)s
<<< Cell parameters:         %(cell_str)s
<<< Friedel's law:           %(anomal_in)s
<<< Merge:                   %(merge_in)s
<<< Dirname:                 %(dirname)s
<<< Wavelength               %(wavelength).4f
"""

fmt_outp = """==> Output file:          %(file_name_out)s
==> Output mode:             %(mode)s
>>> Friedel's law:           %(anomal_out)s
>>> Merge:                   %(merge_out)s
>>> Resolution limit:        %(res_low).2f - %(res_high).2f Angstroem

>>> From Wavelength          %(wavelength).4f
     ??? Guessed Atome Edge       %(ha_name)s
     ???         F'               %(fp).2f
     ???         F"               %(fpp).2f
"""

xdsconv_script = """UNIT_CELL_CONSTANTS=  %(cell_str)s
SPACE_GROUP_NUMBER=   %(spgn_in)s
INPUT_FILE=%(file_name_in)s %(file_type)s 100 1.0
OUTPUT_FILE=%(dir_mode)s%(file_name_out)s %(mode)s FRIEDEL'S_LAW=%(anomal_out)s MERGE=%(merge_out)s
"""

xprep_script = """%(file_name_out)s
4
%(cell_str)s\n\n
H\nB\nS\nS\n\n\n
D\nS\nN\nE
A\nA\nA
\n\n\n\n
%(ident)s\nY\n
Pt\n4\n\n
S\nN\n\nE
P\nR\n15 3
P\n%(ident)s.pattX.ps\nX\n\n\n
P\n%(ident)s.pattY.ps\nY\n\n\n
P\n%(ident)s.pattZ.ps\nZ\n\n\n
E\nQ
"""

f2mtz_script = """
TITLE data from XDS 
FILE %(file_name_out)s
SYMMETRY    %(spgn_in)s
CELL    %(cell_str)s
LABOUT  H K L FP%(lbl)s SIGFP%(lbl)s %(cinp_ano)s FreeRflag
CTYPOUT H H H  F   Q    D      Q     %(cinp_ano2)s  X
END
"""

solve_script = """#!/bin/sh
set noclobber

# run solve

solve << eof-solve > solve.out
symfile %(spg_name)s.sym
cell  %(cell_str)s
resolution %(res_low).3f %(res_high).3f
logfile ./solve.logfile

readformatted
unmerged
read_intensities

checksolve
mad_atom %(ha_name)s

lambda 1 
label sad wavelength = %(wavelength)s
rawmadfile %(file_name_out)s
fixscattfactors
fprprv_mad %(fpp)s
nsolsite_deriv 4
acceptance 0.35

!nres 200

SAD
eof-solve

# run resolve

resolve << eof-resolve > resolve.out
!solvent_content 0.30
!seq_file protein.seq
eof-resolve

fft HKLIN resolve.mtz MAPOUT resolve.ccp4map << eof-fft > resolve_fft.out 
TITLE    from resolve
LABIN   F1=FP PHI=PHIM W=FOMM
END
eof-fft
"""

crossec_script = """
ATOM %(ha_name)s
NWAV 1 %(wavelength)s
END
"""

cns_par = """a, b, c, alpha, beta, gamma = %(cns_cell)s
sg = %(cns_sg)s
low_res, high_res = %(res_low).3f %(res_high).3f
reflection_infile_1 = %(file_name_out)s
obs_f, obs_sigf, = "fobs", "sigma"
test_set, test_flag = "test", 1
"""

cad_script = "LABIN FILE 1 ALL\nEND\n"

HAd = {0.72227:("U ",92,"L3"),0.72766:("Y ",39,"K "),0.76973:("Sr",38,"K "),
0.81554:("Rb",37,"K "),0.86552:("Kr",36,"K "),0.92040:("Br",35,"K "),
0.92340:("Bi",83,"L3"),0.95073:("Pb",82,"L3"),0.97974:("Se",34,"K "),
1.00910:("Hg",80,"L3"),1.04000:("Au",79,"L3"),1.04500:("As",33,"K "),
1.07230:("Pt",78,"L3"),1.10580:("Ir",77,"L3"),1.11658:("Ge",32,"K "),
1.14080:("Os",76,"L3"),1.17730:("Re",75,"L3"),1.19580:("Ga",31,"K "),
1.21550:("W ",74,"L3"),1.25530:("Ta",73,"L3"),1.28340:("Zn",30,"K "),
1.29720:("Hf",72,"L3"),1.34050:("Lu",71,"L3"),1.38059:("Cu",29,"K "),
1.38620:("Yb",70,"L3"),1.43340:("Tm",69,"L3"),1.48350:("Er",68,"L3"),
1.48807:("Ni",28,"K "),1.53680:("Ho",67,"L3"),1.59160:("Dy",66,"L3"),
1.60815:("Co",27,"K "),1.64970:("Tb",65,"L3"),1.71170:("Gd",64,"L3"),
1.74346:("Fe",26,"K "),1.77610:("Eu",63,"L3"),1.84570:("Sm",62,"L3")}

# A dictionary containing all sg_number:sg_name from the spacegroup.lib
cns_sg_lib = {1:"P1",2:"P-1",3:"P2",4:"P2(1)",5:"C2",6:"PM",7:"PC",8:"CM",
9:"CC",10:"P2/M",11:"P2(1)/M",12:"C2/M",13:"P2/C",14:"P2(1)/C",15:"C2/C",
16:"P222",17:"P222(1)",18:"P2(1)2(1)2",19:"P2(1)2(1)2(1)",20:"C222(1)",
21:"C222",22:"F222",23:"I222",24:"I2(1)2(1)2(1)",25:"PMM2",26:"PMC2(1)",
27:"PCC2",28:"PMA2",29:"PCA2(1)",30:"PNC2",31:"PMN2(1)",32:"PBA2",
33:"PNA2(1)",34:"PNN2",35:"CMM2",36:"CMC2(1)",37:"CCC2",38:"AMM2",
39:"ABM2",40:"AMA2",41:"ABA2",42:"FMM2",43:"FDD2",44:"IMM2",45:"IBA2",
46:"IMA2",47:"PMMM",48:"PNNN",49:"PCCM",50:"PBAN",51:"PMMA",52:"PNNA",
53:"PMNA",54:"PCCA",55:"PBAM",56:"PCCN",57:"PBCM",58:"PNNM",59:"PMMN",
60:"PBCN",61:"PBCA",62:"PNMA",63:"CMCM",64:"CMCA",65:"CMMM",13065:"AMMM",
66:"CCCM",67:"CMMA",68:"CCCA",69:"FMMM",70:"FDDD",71:"IMMM",72:"IBAM",
73:"IBCA",74:"IMMA",75:"P4",76:"P4(1)",77:"P4(2)",78:"P4(3)",79:"I4",
80:"I4(1)",81:"P-4",82:"I-4",83:"P4/M",84:"P4(2)/M",85:"P4/N",86:"P4(2)/N",
87:"I4/M",88:"I4(1)/A",89:"P422",90:"P42(1)2",91:"P4(1)22",92:"P4(1)2(1)2",
93:"P4(2)22",94:"P4(2)2(1)2",95:"P4(3)22",96:"P4(3)2(1)2",97:"I422",
98:"I4(1)22",99:"P4MM",100:"P4BM",101:"P4(2)CM",102:"P4(2)NM",103:"P4CC",
104:"P4NC",105:"P4(2)MC",106:"P4(2)BC",107:"I4MM",108:"I4CM",109:"I4(1)MD",
110:"I4(1)CD",111:"P-42M",112:"P-42C",113:"P-42(1)M",114:"P-42(1)C",
115:"P-4M2",116:"P-4C2",117:"P-4B2",118:"P-4N2",119:"I-4M2",120:"I-4C2",
121:"I-42M",122:"I-42D",123:"P4/MMM",124:"P4/MCC",125:"P4/NBM",126:"P4/NNC",
127:"P4/MBM",128:"P4/MNC",129:"P4/NMM",130:"P4/NCC",131:"P4(2)/MMC",
132:"P4(2)/MCM",133:"P4(2)/NBC",134:"P4(2)/NNM",135:"P4(2)/MBC",
136:"P4(2)/MNM",137:"P4(2)/NMC",138:"P4(2)/NCM",139:"I4/MMM",140:"I4/MCM",
141:"I4(1)/AMD",142:"I4(1)/ACD",143:"P3",144:"P3(1)",145:"P3(2)",146:"R3",
20146:"R3R",147:"P-3",148:"R-3",20148:"R-3R",149:"P312",150:"P321",
151:"P3(1)12",152:"P3(1)21",153:"P3(2)12",154:"P3(2)21",155:"R32",
20155:"R32R",156:"P3M1",157:"P31M",158:"P3C1",159:"P31C",160:"R3M",
20160:"R3MR",161:"R3C",20161:"R3CR",162:"P-31M",163:"P-31C",164:"P-3M1",
165:"P-3C1",166:"R-3M",20166:"R-3MR",167:"R-3C",20167:"R-3CR",168:"P6",
169:"P6(1)",170:"P6(5)",171:"P6(2)",172:"P6(4)",173:"P6(3)",174:"P-6",
175:"P6/M",176:"P6(3)/M",177:"P622",178:"P6(1)22",179:"P6(5)22",
180:"P6(2)22",181:"P6(4)22",182:"P6(3)22",183:"P6MM",184:"P6CC",
185:"P6(3)CM",186:"P6(3)MC",187:"P-6M2",188:"P-6C2",189:"P-62M",
190:"P-62C",191:"P6/MMM",192:"P6/MCC",193:"P6(3)/MCM",194:"P6(3)/MMC",
195:"P23",196:"F23",197:"I23",198:"P2(1)3",199:"I2(1)3",200:"PM-3",
201:"PN-3",202:"FM-3",203:"FD-3",204:"IM-3",205:"PA-3",206:"IA-3",
207:"P432",208:"P4(2)32",209:"F432",210:"F4(1)32",211:"I432",
212:"P4(3)32",213:"P4(1)32",214:"I4(1)32",215:"P-43M",216:"F-43M",
217:"I-43M",218:"P-43N",219:"F-43C",220:"I-43D",221:"PM-3M",
222:"PN-3N",223:"PM-3N",224:"PN-3M",225:"FM-3M",226:"FM-3C",227:"FD-3M",
228:"FD-3C",229:"IM-3M",230:"IA-3D"}

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r
    
def opWriteCl(filename, _str):
    f = open(filename,"w")
    r = f.write(_str)
    f.close()

def guessHA(wavelength):
    "Trying to find the closest HA edge in HAd from wavelength"
    diffs = map(lambda x: abs(x - wavelength),HAd.keys())
    return HAd[HAd.keys()[diffs.index(min(diffs))]]

def get_crossec(P):
    opWriteCl("crossec.inp", crossec_script % vars(P))
    try:
        os.system("crossec<crossec.inp>crossec.out")
    except: pass
    txt = opReadCl("crossec.out")
    idx = txt.index(" $$\n<B>")
    if idx != -1:
        try:
            fp, fpp = map(float, txt[idx-23:idx-2].split())
            return fp, fpp
        except: pass
    else: return 0., 0.

    
class Dumy: pass

class DoMode:
    def __init__(self, P):
        self.mode = P.mode
        self.dir_mode = string.lower(self.mode)+"/"
        P.dir_mode = self.dir_mode
        P.ha_name = ""
        try:
            from xdsHKLinfos import get_info
            infos = get_info(P.file_name_in)
            # print infos
            P.res_high =  infos["res_high"]
            P.res_low = infos["res_low"]
        except: 
            P.res_high = 2.2
            P.res_low = 30
        if P.wavelength:
            self.wavelength = P.wavelength
            self.HA = guessHA(self.wavelength)
            P.ha_name = self.HA[0]
            P.fp, P.fpp = get_crossec(P)
        else:
           self.wavelength = 1.5418
           P.ha_name = Unknown
           P.fp, P.fpp = 0., 0.
        if self.mode == "CNS":
            self.name_ext = ".cv"
            P.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
        if self.mode == "SHELX":
            self.name_ext = ".hklf4"
        if self.mode == "SOLVE":
            self.name_ext = ".hkli"
            #P.anomal_out = "FALSE"
        if self.mode == "CCP4":
            self.name_ext = ".mtz"
            P.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
            if P.anomal_out == "FALSE":
                P.cinp_ano = "DANO%(lbl)s SIGDANO%(lbl)s ISYM%(lbl)s" \
                              % vars(P)
                P.cinp_ano2 = "Y"
            else: P.cinp_ano = P.cinp_ano2 = ""
            P.merge_out = "TRUE"
        
        P.ident = ".".join(P.file_name_in.split(".")[:-1])
        #P.file_name_out = P.ident+"_"+P.ha_name+self.name_ext
        P.file_name_out = P.ident + self.name_ext
        if self.mode == "CCP4":
            self.last_name = P.file_name_out
            P.file_name_out = "F2MTZ.HKL.TMP"
        opWriteCl("XDSCONV.INP", xdsconv_script % vars(P) + P.free_out)
        
    def run(self):
        ls = os.listdir(".")
        if not self.dir_mode[:-1] in ls:
            os.makedirs(self.dir_mode)
        if self.mode not in ("SOLVE",):
            os.system("xdsconv")
            #os.chdir(P.dirmode)
     
    def post_run(self, P):
        if self.mode == "SHELX":
            # run xprepx
            opWriteCl("%s/xprepx.inp" % P.dir_mode, xprep_script % vars(P))
            os.system("cd %s;xprepx<xprepx.inp>xprep.log" % P.dir_mode)
        elif self.mode == "CNS":
            #file_name_out = P.dirmode+dirname+".cv"
            P.cns_sg = cns_sg_lib[int(XC.spgn_in)]
            P.cns_cell = 6*"%.2f, " % tuple(P.cell)
            opWriteCl("%s/cns_xtal.par" % P.dir_mode, cns_par % vars(P))
        elif self.mode == "CCP4":
            opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_script % vars(P))
            opWriteCl("%s/cad.inp" % P.dir_mode, cad_script)
            os.system("cd %s;f2mtz hklout TMP.MTZ<f2mtz.inp >f2mtz.log" % P.dir_mode)
            os.system("cd %s;cad hklin1 TMP.MTZ hklout %s <cad.inp >cad.log"
                          % (P.dir_mode, self.last_name))
            #os.system("rm -f F2MTZ.INP ccp4/TMP.MTZ ccp4/F2MTZ.HKL.TMP")
        elif self.mode == "SOLVE":
            os.system("grep -v \! %s > solve/%s" % \
                                      (P.file_name_in, P.file_name_out))
            P.spg_name = (xupy.SPGlib[int(P.spgn_in)][0]).lower()
            opWriteCl("%s/solve.sh" % P.dir_mode, solve_script % vars(P))
            print 

import xupy
import string, sys, os

if __name__=='__main__':
    #### Default values
    XC = Dumy()
    
    XC.file_type = "XDS_ASCII"
    XC.anomal_out = ""
    XC.merge_out = ""
    XC.cell = ""
    XC.lbl = ""
    XC.free_out = ""

    if len(sys.argv) == 1:
        print usage % (sys.argv[0], "|".join(options))
        sys.exit()
    if sys.argv.count("-a"):
        sys.argv.remove("-a")
        XC.anomal_out = "FALSE"
    if sys.argv.count("-n"):
        sys.argv.remove("-n")
        XC.anomal_out = "TRUE"
    if sys.argv.count("-u"):
        sys.argv.remove("-u")
        XC.merge_out = "FALSE"
    if sys.argv.count("-m"):
        sys.argv.remove("-m")
        XC.merge_out = "TRUE"
    if sys.argv.count("-f"):
        sys.argv.remove("-f")
        XC.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
    for arg in sys.argv:
        if arg.count("-l="):
            XC.lbl = "_"+arg[3:]
            sys.argv.remove(arg)
    if len(sys.argv) > 3:
        # test input file type of inherited reflections
        add = "INHERIT_TEST_REFLECTIONS_FROM_FILE=%s CNS\n" % (sys.argv[3])
        xdsconv_script += add
    
    XC.dirname = os.path.split(os.getcwd())[1]
    XC.file_name_in = sys.argv[1]
    sys.argv.remove(XC.file_name_in)
    
    if len(sys.argv) >= 2:
        mode = string.upper(sys.argv[1])
        XC.mode = mode
    else: mode =  ""
    if not options.count(mode):
        print "\n>>> Error !!! Unknown mode: ",mode," ???"
        print usage % (sys.argv[0], "|".join(options))
        sys.exit()
    
    H = xupy.read_xdsascii_head(XC.file_name_in)
    XC.spgn_in = H["sym"]
    XC.cell = map(float,H["cell"].split())
    XC.cell_str = 6*"%.2f  " % tuple(XC.cell)
    XC.anomal_in = H["friedels_law"]
    XC.merge_in = H["merge"]
    
    if not H["wavelength"]:
    # Try to catch the wavelength from 
        if H["inputfile_name"] \
           and os.path.exists(H["inputfile_name"]) \
           and H["inputfile_type"] == "XDS_ASCII":
               _h = xupy.read_xdsascii_head(H["inputfile_name"])
               if _h["wavelength"] : XC.wavelength = _h["wavelength"]
    else: XC.wavelength = H["wavelength"]
    if H["friedels_law"] == "TRUE" and XC.anomal_out == "FALSE":
        print "\n>>> WARNING!  The input file does not contain Friedel's mate."
    if H["merge"] == "TRUE" and XC.merge_out == "FALSE":
        print "\n>>> WARNING!  The input file does not unmerged reflections."
       
    if not XC.anomal_out: XC.anomal_out = H["friedels_law"]
    if not XC.merge_out: XC.merge_out = H["merge"]
    
    # Depending on the chosen mode, predefine suffix_name, merge_out, anomal_out...

    print fmt_inp % vars(XC)
            
    #-----------------------------
    
    E = DoMode(XC)
    print fmt_outp % vars(XC)
    E.run()
    E.post_run(XC)
    
    sys.exit()
    
def get_wavelength():
    try: xparm = open("../XPARM.XDS").readlines()
    except: return 1.5418
    return string.atof(xparm[1].split()[0])

def get_amore_symop(symop_file,sg_num):
    patt = "# %3s " % sg_num
    symopF = open(symop_file,"r").read()
    start1 = string.index(symopF,patt)
    start2 = string.index(symopF,"x,y,z",start1)
    stop  = string.index(symopF,"end",start2)
    print symopF[start1:start2]
    return symopF[start2:stop+4]

def read_cns(file_name_in):
    import struct
    nskip = 3
    fmt_line = '6x5s5s5s6x12s7x12s1x'
    raw = open(file_name_in,"r")
    nbr_refl = string.atoi(string.split(raw.readline())[1])
    for a in range(nskip) : raw.readline(),
    raw_read   = raw.read(struct.calcsize(fmt_line) * nbr_refl)
    raw_unpack = struct.unpack((fmt_line* nbr_refl), raw_read)
    data = reshape(array(map(string.atof,raw_unpack)),(nbr_refl,5))
    raw.close()
    return data
    
def fall2shelx_ano(file_name_in,file_name_out):

    """ Reads a fall anomalous file and output a
        shelx anomalous file (H,K,L,DANO,SIG(DANO))
    """
    #from Numeric import *
    import  datafile, Numeric
    d = datafile.datafile()
    d.read(file_name_in)
    print Numeric.shape(d.data)
    d.data = Numeric.concatenate((d.data[:,:3],d.data[:,5:7]),1)
    fitre = Numeric.greater_equal(d.data[:,-1],0)
    d.data = Numeric.compress(fitre,d.data,0)
    fmt = "%4d%4d%4d%8.2f%8.2f\n" * d.data.shape[0]
    print Numeric.shape(d.data)
    open(file_name_out,"w").write(fmt % tuple(Numeric.ravel(d.data)))

def fall2ccp4_ano(file_name_in,file_name_out):
    """ Reads a fall anomalous file and output a
        pre-ccp4 anomalous file (H,K,L,F,SIG(F),DANO,SIG(DANO),ISYM)
    """
    #from Numeric import *
    import  datafile, Numeric
    d = datafile.datafile()
    d.read(file_name_in)
    isym = (Numeric.less(d.data[:,-1],0)*Numeric.less(d.data[:,-2],0)*2 + 
            Numeric.less(d.data[:,-1],0)*Numeric.greater(d.data[:,-2],0)+
            Numeric.equal(d.data[:,-1],0))
    print Numeric.shape(d.data)
    d.data = Numeric.concatenate((d.data[:,:3],d.data[:,-4:]),1)
    isym = Numeric.transpose(Numeric.array([isym]))
    d.data = Numeric.concatenate((d.data,isym),1)
    #fitre = Numeric.greater_equal(d.data[:,-1],0)
    #d.data = Numeric.compress(fitre,d.data,0)
    fmt = "%4d%4d%4d%13.5e%13.5e%13.5e%13.5e%2d\n" * d.data.shape[0]
    #print Numeric.shape(d.data)
    open(file_name_out,"w").write(fmt % tuple(Numeric.ravel(d.data)))

def cns2mtz(file_name_in,file_name_out,cell,sym):
    """ Reads a cns file and output a
        mnz file (H,K,L,F,SIG(F))
        with both F+ and F-.
    """
    inp = """#!/bin/csh -f
f2mtz hklin %s hklout %s << EOF
    SKIP 4
    TITLE mtz from xds, xscale and xdsconv
    CELL %s
    SYMM %s
    LABOUT H   K  L   FP SIGFP  
    CTYPE  H   H  H   F  Q 
    FORMAT '(6X,3F5.0,6X,F12.4,7X,F12.4)'
EOF\n""" % (file_name_in,file_name_out,cell,sym)
    open("xds2mtz.com","wb").write(inp)
    os.system("csh xds2mtz.com")

def cns2mtz_ano(file_name_in,file_name_out,cell,sym):
    """ Reads a cns file and output a
        mnz file (H,K,L,F,SIG(F))
        with both F+ and F-.
    """
    inp = """#!/bin/csh -f
f2mtz hklin %s hklout %s << EOF
    TITLE mtz from xds, xscale and xdsconv
    CELL %s
    SYMM %s
    LABOUT H   K  L   FP SIGFP DANO SANO ISYM 
    CTYPE  H   H  H   F  Q     D    Q    Y
    FORMAT '(3F4.0,E13.5,E13.5,E13.5,E13.5,F2.0)'
EOF\n""" % (file_name_in,file_name_out,cell,sym)
    open("xds2mtz.com","wb").write(inp)
    os.system("csh xds2mtz.com")

def mtz2shelx(file_name_in,file_name_out):
    inp = """#!/bin/csh -f
mtz2various hklin %s hklout %s << EOF
    OUTPUT SHELX
    LABIN  FP=DANO SIGFP=SANO
EOF\n""" % (file_name_in,file_name_out)
    open("mtz2shelx.com","wb").write(inp)
    os.system("csh mtz2shelx.com")

def cns2SnB_ano(file_name_in,file_name_out):
    """ Reads a cns anomalous file and output a
        SnB anomalous file (H,K,L,F,SIG(F))
        with both F+ and F-.
        Can divide the sigma, in order to avoid
        to many rejection by SnB...
    """
    devid_factor = 1.0
    data = read_cns(file_name_in)
    data[:,-1] = data[:,-1] / devid_factor
    fmt = "%5.0f%5.0f%5.0f%14.3f%14.3f\n" * shape(data)[0]
    open(file_name_out,"w").write(fmt % tuple(ravel(data)))


if mode == "AMORE":
    post_run      = "AMORE"
    mode          = "MIR"
    anomal_out  = "TRUE"
    file_name_out = dirmode+"hkl.tmp"
    
elif mode == "EPMR":
    post_run      = "EPMR"
    mode          = "MIR"
    anomal_out  = "TRUE"
    file_name_out = dirmode+"hkl"
        
elif mode == "CNS":
    file_name_out = dirmode+dirname+".cv"
    script = script+"GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n\n"
    cns_par = 'a, b, c, alpha, beta, gamma = '+string.join(file_cc.split(),", ")
    cns_par += '\nsg = "%s"\n' % cns_sg_lib[string.atoi(spgn_in)]
    cns_par += 'low_res, high_res = 1000.0, 0.9\n'
    cns_par += 'reflection_infile_1 = "%s' % dirname+'.cv"\n'
    cns_par += 'obs_f, obs_sigf, = "fobs", "sigma"\n'
    cns_par += 'test_set, test_flag = "test", 1\n'
    opWriteCl("cns/cns_xtal.par", cns_par)
    
elif mode == "SNB":
    mode          = "CNS"
    file_name_out = dirmode+"hkl.tmp"
    post_run      = "SNB_ANO"
    
elif mode == "FALL":
    file_name_out = dirmode+dirname+".fall"
    
elif mode == "SHELX" and anomal_out == "TRUE":
    file_name_out = dirmode+"shelxl.hkl"
    
elif mode == "SHELX" and anomal_out == "FALSE":
    post_run      = "SHELX_ANO"
    mode          = "FALL"
    file_name_out = dirmode+"hkl.tmp"


#========================================================================#
#======        POST PROCESSING    ====    POST PROCESSING        ========#
#========================================================================#

#=====================#
if post_run == "AMORE":
#=====================#
    os.system("""echo "FORMAT (3I5,2E12.4)" > format.dat""")
    os.system("cat format.dat hkl.tmp > hkl.d ; rm -f format.dat")
    amore_data_fmt = " ** %s **\n%s\n%s 0\n95.  0.\n15. 3.5\n1   1\n"
    symop = get_amore_symop(os.environ["AMORE"]+"/symop.d",spgn_in)
    opWriteCl("data.d",(amore_data_fmt % (dirname,cell,symop)))

#====================#
if post_run == "EPMR":
#====================#
    opWriteCl("cell", "%s  %s\n" % (cell,spgn_in))

#====================#
if post_run == "CCP4":
#====================#
    print "\nOUTPUT CCP4: F,SIG(F)n"
    ccp4_data = """
set cell = '%s'\nset spacegroup =  %s\nset sym =  $spacegroup\nset resmin_m = %.2f\n
set resmin = %.2f\nset resmax = %.2f\nset resol  = "$resmin_m $resmax"
set labin = "FP=FP SIGFP=SIGFP FREE=FreeR_flag"\n"""
    opWriteCl("source.dat", ccp4_data % (cell,spgn_in,500,500,1.0))
    cns2mtz("hkl.tmp",dirname+".mtz",cell,spgn_in)

#========================#
if post_run == "CCP4_ANO":
#========================#
    print "\nOUTPUT CCP4_ANO : F,SIG(F),DANO,SIG(DANO),ISYM\n"
    fall2ccp4_ano("hkl.tmp","hkl.tmp")
    ccp4_data = """
set cell = '%s'\nset spacegroup =  %s\nset sym =  $spacegroup\nset resmin_m = %.2f\n
set resmin = %.2f\nset resmax = %.2f\nset resol  = "$resmin_m $resmax"
set labin = "FP=FP SIGFP=SIGFP DP=DANO SIGDP=SANO ISYM=ISYM"\n"""
    open("source.dat","wb").write(ccp4_data % (cell,spgn_in,500,500,1.0))
    cns2mtz_ano("hkl.tmp",dirname+"_ano.mtz",cell,spgn_in)

#=========================#
if post_run == "SHELX_ANO":
#=========================#
    shelxinp = ""
    wavelength = get_wavelength()
    print wavelength
    fall2ccp4_ano("hkl.tmp","hkl.tmp")
    cns2mtz_ano("hkl.tmp",dirname+"_ano.mtz",cell,spgn_in)
    mtz2shelx(dirname+"_ano.mtz",dirname+"_ano_mtz.hkl")
    print "\nOUTPUT SHELX_ANO : F,SIG(F),DANO,SIG(DANO)\n"
    raw = open(dirname+"_ano_mtz.hkl")
    line = raw.readline()
    while line[0] != " ":
        if line[:4] == "CELL": 
            line = "CELL %.5f %s\n" % (wavelength,cell)
        print line,
        shelxinp += line
        line = raw.readline()
    shelxinp += "SFAC SE\nUNIT 10\n"
    shelxinps = shelxinp + "OMIT 0\nPATT 20 30 3\nHKLF 3\n"
    shelxinpd = shelxinp + "PATS\nFIND 8\nMIND -3.5\nNTRY 100\nHKLF 3\n"
    opWriteCl("shelxs.ins", shelxinps)
    opWriteCl("shelxd.ins", shelxinpd)
    fall2shelx_ano("hkl.tmp","shelxd.hkl")
    if "shelxs.hkl" in os.listdir("."):
        os.remove("shelxs.hkl")
    os.symlink("shelxd.hkl","shelxs.hkl")
    
#=======================#
if post_run == "SNB_ANO":
#=======================#
    print "\nOUTPUT SNB_ANO : F,SIG(F)\n"
    cns2SnB_ano("hkl.tmp",dirname+"_ano.snb")

os.system("rm -f hkl.tmp")
