#!/usr/bin/env python
"""
Last version
01/10/02 legrand@lccp.ibs.fr
 TODO:
    - calculate patterson normal and anomal
    - Add this for cns2mtz with same Free reflections.
    SKIP 5
    TITLE mtz from xds, xscale and xdsconv
    CELL 154.55 154.55 154.55 90 90 90
    SYMM 208
    LABOUT H   K  L   FP SIGFP   FreeR_flag
    CTYPE  H   H  H   F  Q       X
    FORMAT '(6X,3F5.0,6X,F12.4,7X,F12.4,7X,F1.0)'
"""
options = ["AMORE","CCP4","CNS","EPMR","FALL","SHELX","SNB"]
usage   = """
>>>   Usage : xdsconv.py filein.hkl [X] option [free_hkl_to_inherit]\n
        option = %s
          if X = N : normal output (defaulf)
          if X = A : anomal output\n
>>>   Cell parameters and space group number are taken
      from the header of input file.
"""
import xupy
import  string,sys,os,struct

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
    from Numeric import *
    import  datafile
    d = datafile.datafile()
    d.read(file_name_in)
    print shape(d.data)
    d.data = concatenate((d.data[:,:3],d.data[:,5:7]),1)
    fitre = greater_equal(d.data[:,-1],0)
    d.data = compress(fitre,d.data,0)
    fmt = "%4d%4d%4d%8.2f%8.2f\n" * d.data.shape[0]
    print shape(d.data)
    open(file_name_out,"w").write(fmt % tuple(ravel(d.data)))

def fall2ccp4_ano(file_name_in,file_name_out):
    """ Reads a fall anomalous file and output a
        pre-ccp4 anomalous file (H,K,L,F,SIG(F),DANO,SIG(DANO),ISYM)
    """
    from Numeric import *
    import  datafile
    d = datafile.datafile()
    d.read(file_name_in)
    isym = (less(d.data[:,-1],0)*less(d.data[:,-2],0)*2 + 
            less(d.data[:,-1],0)*greater(d.data[:,-2],0)+
            equal(d.data[:,-1],0))
    print shape(d.data)
    d.data = concatenate((d.data[:,:3],d.data[:,-4:]),1)
    isym = transpose(array([isym]))
    d.data = concatenate((d.data,isym),1)
    #fitre = greater_equal(d.data[:,-1],0)
    #d.data = compress(fitre,d.data,0)
    fmt = "%4d%4d%4d%13.5e%13.5e%13.5e%13.5e%2d\n" * d.data.shape[0]
    print shape(d.data)
    open(file_name_out,"w").write(fmt % tuple(ravel(d.data)))

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

#### Default values
file_type    = "XDS_ASCII"
outp_anomal  = "TRUE"
cell         = ""
post_run     = ""

if sys.argv.count("A") or sys.argv.count("a"):
    sys.argv.remove("A")
    outp_anomal = "FALSE"
if sys.argv.count("N") or sys.argv.count("n"):
    sys.argv.remove("N")
    outp_anomal = "TRUE"
    
dirname = os.path.split(os.getcwd())[1]
file_name_in = sys.argv[1]
sys.argv.remove(file_name_in)
if len(sys.argv) >= 2:
    mode = string.upper(sys.argv[1])
    dirmode = string.lower(mode)+"/"
else: mode =  ""
if not options.count(mode):
    print "\n>>> Erreur !!! Mode ",mode," ???"
    print usage % string.join(options,"|")
    sys.exit()


file_cc, file_spaceGroup, file_anomal, merge = xupy.read_xdsascii_head(file_name_in)
a, b, c, alpha, beta, gamma = file_cc.split()
for p in file_cc.split() : cell = cell + (("%.2f  ") % string.atof(p))
if file_anomal == "TRUE" and outp_anomal == "FALSE":
     print "\n!!! Error, The input file does not contain Friedel's mate"
     sys.exit()
else: outp_anomal = file_anomal

print """
>>> Format choisi        =  %s 
>>> Num. groupe d'espace =  %s
>>> Maille               =  %s
>>> Repertoir            =  %s
""" % (mode,file_spaceGroup,cell,dirname)

inp = open("XDSCONV.INP","wb")
script = """
UNIT_CELL_CONSTANTS=  %(cell)s
SPACE_GROUP_NUMBER=   %(file_spaceGroup)s
INPUT_FILE=  %(file_name_in)s  %(file_type)s 100 1.0
OUTPUT_FILE= %(file_name_out)s  %(mode)s  FRIEDEL'S_LAW= %(outp_anomal)s
"""
if len(sys.argv) >= 3:
    print sys.argv[2]
    add = "INHERIT_TEST_REFLECTIONS_FROM_FILE=%s CNS\n" % (sys.argv[2])
    script = script + add

ls = os.listdir(".")
if not dirmode[:-1] in ls:
    os.makedirs(dirmode)

if mode == "AMORE":
    post_run      = "AMORE"
    mode          = "MIR"
    outp_anomal  = "TRUE"
    file_name_out = dirmode+"hkl.tmp"
    
elif mode == "EPMR":
    post_run      = "EPMR"
    mode          = "MIR"
    outp_anomal  = "TRUE"
    file_name_out = dirmode+"hkl"
    
elif mode == "CCP4" and outp_anomal == "TRUE":
    post_run      = "CCP4"
    mode          = "CNS"
    file_name_out = dirmode+"hkl.tmp"

elif mode == "CCP4" and outp_anomal == "FALSE":
    post_run      = "CCP4_ANO"
    mode          = "FALL"
    file_name_out = dirmode+"hkl.tmp"
    
elif mode == "CNS":
    file_name_out = dirmode+dirname+".cv"
    script = script+"GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n\n"
    cns_par = 'a, b, c, alpha, beta, gamma = '+string.join(file_cc.split(),", ")
    cns_par += '\nsg = "%s"\n' % cns_sg_lib[string.atoi(file_spaceGroup)]
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
    
elif mode == "SHELX" and outp_anomal == "TRUE":
    file_name_out = dirmode+"shelxl.hkl"
    
elif mode == "SHELX" and outp_anomal == "FALSE":
    post_run      = "SHELX_ANO"
    mode          = "FALL"
    file_name_out = dirmode+"hkl.tmp"

inp.write(script % vars())
inp.close()

os.system("xdsconv")
os.chdir(dirmode)

#========================================================================#
#======        POST PROCESSING    ====    POST PROCESSING        ========#
#========================================================================#

#=====================#
if post_run == "AMORE":
#=====================#
    os.system("""echo "FORMAT (3I5,2E12.4)" > format.dat""")
    os.system("cat format.dat hkl.tmp > hkl.d ; rm -f format.dat")
    amore_data_fmt = " ** %s **\n%s\n%s 0\n95.  0.\n15. 3.5\n1   1\n"
    symop = get_amore_symop(os.environ["AMORE"]+"/symop.d",file_spaceGroup)
    opWriteCl("data.d",(amore_data_fmt % (dirname,cell,symop)))

#====================#
if post_run == "EPMR":
#====================#
    opWriteCl("cell", "%s  %s\n" % (cell,file_spaceGroup))

#====================#
if post_run == "CCP4":
#====================#
    print "\nOUTPUT CCP4: F,SIG(F)n"
    ccp4_data = """
set cell = '%s'\nset spacegroup =  %s\nset sym =  $spacegroup\nset resmin_m = %.2f\n
set resmin = %.2f\nset resmax = %.2f\nset resol  = "$resmin_m $resmax"
set labin = "FP=FP SIGFP=SIGFP FREE=FreeR_flag"\n"""
    opWriteCl("source.dat", ccp4_data % (cell,file_spaceGroup,500,500,1.0))
    cns2mtz("hkl.tmp",dirname+".mtz",cell,file_spaceGroup)

#========================#
if post_run == "CCP4_ANO":
#========================#
    print "\nOUTPUT CCP4_ANO : F,SIG(F),DANO,SIG(DANO),ISYM\n"
    fall2ccp4_ano("hkl.tmp","hkl.tmp")
    ccp4_data = """
set cell = '%s'\nset spacegroup =  %s\nset sym =  $spacegroup\nset resmin_m = %.2f\n
set resmin = %.2f\nset resmax = %.2f\nset resol  = "$resmin_m $resmax"
set labin = "FP=FP SIGFP=SIGFP DP=DANO SIGDP=SANO ISYM=ISYM"\n"""
    open("source.dat","wb").write(ccp4_data % (cell,file_spaceGroup,500,500,1.0))
    cns2mtz_ano("hkl.tmp",dirname+"_ano.mtz",cell,file_spaceGroup)

#=========================#
if post_run == "SHELX_ANO":
#=========================#
    shelxinp = ""
    wavelength = get_wavelength()
    print wavelength
    fall2ccp4_ano("hkl.tmp","hkl.tmp")
    cns2mtz_ano("hkl.tmp",dirname+"_ano.mtz",cell,file_spaceGroup)
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
