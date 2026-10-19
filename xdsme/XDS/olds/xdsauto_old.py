#!/usr/bin/env python
# TODO.
# Verifications:
#   - verify the input parameter from auto.par (creat a dic of alowed keys)
#       give warning for unalowed keys.
#   - try the xds_par or xds, using which (/...xds or /...xds_par)
#   - Make a test procedur.
# Step III and IV: If more than 1 sector (see lookup_frame())
#               III  -> Profile FIT with each sector with (SPOT RANGE...)
#               III  -> Integrate separetly each sector (INTEGRATE)
#                IV  -> Merge and reindex with XSCALE,,,
# look at statistics of occuring space-group to give hint when systematic (see pdb)
#  abscence are not 

"""
Try to automate XDS data processing
"""

import string, sys, shutil
import os, re
from xupy import *
import BeamLines

__version__ = "0.5.2"
__author__ = "Pierre Legrand (legrand@emble-grenoble.fr)"
__date__ = "26-05-2003"
__copyright__ = "Copyright (c) 2003 Pierre Legrand"
__license__ = "GPL"

_job_angle = 4., 6., 10., 15, 15
_nopt = 2         # Default number of profil optimization steps.
_idxref_fit_limit = 40
_nspot_slice = 2  # number of frames/sequence in COLSPOT and Init. Prof determination.
_wcut_rmea3 = 1.2
_nspot_min = 200  # Minimium number of spots for autoindexing.


def read_adpParams(inp_str):
    """ Read ADP input parameters """
    blk_xds = re.compile('^[ \t]*$')  # Blank line
    par_s = ""
    for line in inp_str.splitlines():
        # removes comments and spaning white space
        sline = re.split('(#)',string.strip(line))[0]
        # skiping blanc line
        if not (blk_xds.match(sline)):
            keys = re.split(r"([a-zA-Z0-9_.\- ]+=)",sline)
            # removes leading empty keys
            if not keys[0]: keys = keys[1:]
            if (len(keys)%2): print "WARNING ... odd number of keys!"
            keys[1] = string.strip(keys[1])
            fmt_v = ' %s\n'
            try: string.atof(keys[1])
            except ValueError:  fmt_v = ' "%s"\n'
            par_s = par_s + ("%s" + fmt_v) % (keys[0],keys[1])
    return par_s

def get_num(seq,num):
    L = []
    for i in num: L.append(seq[i])
    return L

def unpack_latticefit(_str):
    ss = _str.split()
    return (float(ss[2]),int(ss[0]),ss[1],map(float,ss[3:9]),
            tuple(map(int,ss[9:])))
            
def get_maxResolution(infile="INTEGRATE.LP"):
    p_lp = opReadCl(infile)[:1115]
    to_find = "RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)"
    resol = p_lp[p_lp.index(to_find)+len(to_find):].split()[1]
    return float(resol)
        
def resum_idxref(idxref="IDXREF.LP"):
    list_latticesFit = []
    i_lp = opReadCl(idxref)
    st1 = i_lp.index("LATTICE-  BRAVAIS-")
    idxref_latticesFit = i_lp[st1+172:st1+5012]
    return map(unpack_latticefit, idxref_latticesFit.splitlines())

def choose_sequence(seq,min_f=None,max_f=None):
    min_c, max_c = 1e10, 1
    min_c, max_c = seq[0][0], seq[-1][-1]
    if min_f: min_c = max(min_f,min_c)
    if max_f: max_c = min(max_f,max_c)
    return [min_c, max_c]

def fix_cell(bravais, cell):
    lattice = bravais[0]
    a, b, c, al, be, ga = cell
    if lattice == "a":
        return a, b, c, al, be, ga 
    if lattice == "m":
        return a, b, c, 90., be, 90.
    if lattice == "o":
        return a, b, c, 90., 90., 90.
    if lattice == "t":
        return (a+b)/2., (a+b)/2., c, 90., 90., 90.
    if lattice == "h":
        a = (a+b)/2.
        return a, a, c, 90., 90., 120.
    if lattice == "c":
        a = (a+b+c)/3.
        return a, a, a, 90., 90., 90.

def select_lattices(limit = _idxref_fit_limit):    
            selected = []
            selection = Param()
            for _latt in resum_idxref():
                 if _latt[0] <= limit and _latt[0] != 0.0:
                     selected.append(_latt)
            return selected

def print_o(s,files):
    from string import join
    if type(s) == type(()) or type(s) == type([]):
        s = string.join(map(str,s))
    if len(files) <= 2:
        for file in files: file.write(s)
    elif len(files) == 3:
        for file in files[:2]: file.write(s)
        files[2].write(s)
        #files[2].write(s.replace("\n","<br>\n"))

def print_o2(t,v,f):
    print_o("  %-49s%24s\n" % (t,v),f)

def print_o3(t,v1,v2,f):
    print_o("\t   %-22s %14s    %-s\n" % (t,v1,v2),f)
    
def print_stage(s,files):
    print_o(s,files[:2])
    print_o("<h2>"+s+"</h2>\n",(files[-1],))
    
def stop(s=""):
    flog.close()
    fhtml.write("</pre></BODY>")
    fhtml.close()
    os.chdir("../")
    saveLastVersion(["xdsauto.log","xdsauto.html"])
    print "  %s  STOP\n" % s
    sys.exit()
  
if __name__=='__main__':
    step = 0
    if "-test" in sys.argv:
        print "test"
        get_maxResolution()
        sys.exit()
    if len(sys.argv) >= 1:
        flog = open("xdsauto.log","w")
        fhtml = open("xdsauto.html","w")
        fhtml.write("<HTML> <HEAD>\n<TITLE>Automated XDS Data Processing")
        fhtml.write('</TITLE>\n</HEAD>\n<BODY BGCOLOR="#ffffff">\n<pre>')
        
        #dirname = os.path.split(os.getcwd())[1]
        f1, f2, f3 = (sys.stdout), (sys.stdout,flog), (sys.stdout,flog,fhtml)
        print_o("\n",f3)
        
        for arg in sys.argv[1:]:
            if arg[:2] == "-O" and len(arg) == 3 and arg[-1].isdigit():
                opt = int(arg[-1])
                if 5 >= opt >= 1:
                    _nopt = opt
                    print_o2("INFO: Number of profil optimization steps",_nopt,f3)
                    del sys.argv[sys.argv.index(arg)]
                else:
                    print "\n  ERROR: Illegal argument for the -O option\n"
                    sys.exit()
        for arg in sys.argv[1:]:
            if len(arg) == 2 and arg[-1].isdigit():
                step = int(arg[-1])
                print_o2("INFO: Jumping directly to step",step,f3)
                
    if step >= 0:
        fmt_shape = "beam divergence / reflecting range"
        fmt_err = "\n  ERROR!  %s\n\n"
        fmt_stg = "\n\n\n\t Stage %25s\n\n"
        fmt_link_xds = '  <A href="xds/%s">XDS log file</A> for this stage.<br>'
        fmt_cell = "%6.1f%6.1f%6.1f"
        parf = "auto.par"
        # Defaults values:
        BRAVAIS = ""
        auto = XParam()
        auto.img_dir = "../img/"
        auto.beamline = "ESRF_BM30A"
        auto.spg = 0
        auto.anomalous = "no"
        auto.resolution_range = "30 0.0"
        
        print_o("\n",f3)
        print_o2("Reading Data collection parameters in file", parf,f3)
        print_o("\n\n",f3)
        try: adpParams_raw = opReadCl(parf)
        except: 
            print_o(fmt_err % ("Can't open file %s\n" % parf),f3)
            stop()        
        try: exec read_adpParams(adpParams_raw) in auto.__dict__
        except:
            print_o(fmt_err%("Problems while reading the file %s\n"%parf),f3)
            stop()
        try: shutil.copyfile("XDS.INP","XDS.INP.old")
        except:  pass

        if not "xds" in os.listdir("."): os.mkdir("xds")
        os.chdir("xds")
        
        if auto.img_dir[0] != "/": auto.img_dir = "../" + auto.img_dir
        if auto.img_dir[-1] != "/": auto.img_dir += "/"
        auto.img_template = auto.img_dir+auto.prefix+"???."+auto.suffix

        keys_show = ["beamline","img_template","frame_first",
                     "frame_last","delta_phi","anomalous","wavelength",
                     "distance","beam_x","beam_y","resolution_range","spg"]
        for par in keys_show:
            if auto[par]: print_o("\t%-20s\t%s\n" % (par, auto[par]),f3)
                
        xpar = xdsInp2Param(inp_str=xdsinp_base)
        BL = getattr(BeamLines,auto.beamline)
        xpar.mix(BL)

           
        seq = lookup_frame_sequences(auto.img_template)
        if not seq:
            print fmt_err%("No frame found with template %s"%auto.img_template)
            stop()
        print_o("\n\n",f3)
        print_o2("INFO: First and last frame number found",\
                       "%d - %d" % (seq[0][0],seq[-1][-1]),f3)
        if len(seq) >= 2:
            print_o2("WARNING! Discontinuous frame sequence",seq,f3)
        collect_range = choose_sequence(seq,min_f=auto.frame_first, \
                                           max_f=auto.frame_last)
        xpar.DATA_RANGE = collect_range
        
        nmin_f = 2*_nspot_slice+1
        if collect_range[1]-collect_range[0]+1 >= nmin_f:
            # use two range ex: i-i+2, f-2,f
            # with f at maximum 90 degre distance of i
            max_frame = min(collect_range[1],
                         collect_range[0]+int(89./auto.delta_phi+_nspot_slice))
            xpar.SPOT_RANGE=(collect_range[0],collect_range[0]+_nspot_slice-1),\
                            (max_frame-_nspot_slice+1,max_frame)
        else:
            xpar.SPOT_RANGE = (collect_range[0], 
                           min(collect_range[0]+nmin_f-1,collect_range[1])),
        
        xpar.STARTING_FRAME = auto.frame_first
        xpar.BACKGROUND_RANGE = auto.frame_first, auto.frame_first+4
        xpar.OSCILLATION_RANGE = auto.delta_phi
        xpar.DELPHI= 4*auto.delta_phi
        xpar.STARTING_ANGLE = auto.phi_init
        xpar.X_RAY_WAVELENGTH = auto.wavelength
        xpar.DETECTOR_DISTANCE = auto.distance
        xpar.NAME_TEMPLATE_OF_DATA_FRAMES = auto.img_template
        xpar.ORGX = auto.beam_y/xpar.QX
        xpar.ORGY = auto.beam_x/xpar.QY
        r = map(float,auto.resolution_range.split())
        xpar.INCLUDE_RESOLUTION_RANGE = max(r), min(r)
        if min(r) != 0.0:
             print_o2("INFO: Truncating the highest resolution to ",\
                                    "%.2f" % min(r),f3)
        xpar.FRIEDELS_LAW = "TRUE" # set to true first for spacegroupe comparaison
        if auto.suffix == "pck":
            xpar.NAME_TEMPLATE_OF_DATA_FRAMES += " CCP4"
        
    if step <= 1:
        Stage = " I -  Peak Search"
        #==================================#
        n = 0
        xdsout = "xds_i_%d.out" % n
        print_stage(fmt_stg % Stage,f3)
        fhtml.write(fmt_link_xds % xdsout)
        
        s = ""
        for r in xpar.SPOT_RANGE: s += "%d to %d and " % tuple(r)
        print_o2("Using images","%s" % s[:-5],f3)
        
        xpar.JOB = "XYCORR", "INIT", "COLSPOT"
        run_xds(xpar,inp_f=None,out_f=xdsout)
        
        initlp = opReadCl("INIT.LP")
        colspotlp = opReadCl("COLSPOT.LP")
        ind1 = initlp.find("MEAN GAIN VALUE")
        ind2 = initlp.find("AVERAGE BACKGROUND COUNTS IN A PIXEL")
        ind3 = colspotlp.find("NUMBER OF DIFFRACTION SPOTS ACCEPTED")
        nspot = int(colspotlp[ind3:ind3+58].split()[-1])
        m_gain = float(initlp[ind1:ind1+51].split()[-1])
        m_count = float(initlp[ind2:ind2+58].split()[-1])
        
        print_o2("Mean Gain value", "%8.2f" % m_gain ,f3)
        print_o2("Average background counts","%8.1f" % m_count ,f3)
        print_o2("Number of diffraction spots found","%12d" % nspot ,f3)
        
        while nspot <= _nspot_min:
            n += 1
            xpar.STRONG_PIXEL -= 2
            xpar.MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT -= 1
            if xpar.STRONG_PIXEL < 3. or \
               xpar.MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT < 5: pass
            print_o("  WARNING: Number of spots found is too low (%d <= %d)\n"\
                                              % (nspot,_nspot_min),f3)
            print_o("  Recollecting spots with minimum I/sigma = %.1f\n\n"\
                                               % xpar.STRONG_PIXEL,f3)
            xpar.JOB = "COLSPOT"
            xdsout = "xds_i_%d.out" % n
            fhtml.write(fmt_link_xds % xdsout)
            run_xds(xpar,inp_f=None,out_f=xdsout)
            colspotlp = opReadCl("COLSPOT.LP")
            nspot = int(colspotlp[ind3:ind3+58].split()[-1])
            
        #print_o2("Number of diffraction spots found","%12d" % nspot ,f3)

    else: opWriteCl("XDS.INP","%s" % xpar)
    
    if step <= 2:
        Stage = " II -  Auto-indexing"
        #==================================#
        xdsout = "xds_ii.out"
        print_stage(fmt_stg % Stage,f3)
        fhtml.write(fmt_link_xds % xdsout)
        
        xpar.JOB = "IDXREF"
        run_xds(xpar,out_f=xdsout)
        
        ilp = opReadCl("IDXREF.LP")
        ilp_err = ilp.find("!!! ERROR !!!")
        if ilp_err >= 0:
            print_o("\n  WARNING %s\n\n" % ilp[ilp_err:].splitlines()[0],f3)
        ilp_ind = ilp.find("INDEXING OF OBSERVED SPOTS IN SPACE GROUP #   1")+54
        indexed = ilp[ilp_ind:ilp_ind+38].split()
        num_indexed, num_total = int(indexed[0]),  int(indexed[3])
        ilp_xyz = ilp.index("SPOT    POSITION (PIXELS)", ilp_ind)+25
        ilp_xyz = ilp[ilp_xyz:ilp_xyz+67].split()
        stdev_xy, stdev_z = float(ilp_xyz[0]), float(ilp_xyz[-1])
        
        reduce_cell = tuple(map(float,\
                        opReadCl("XPARM.XDS").splitlines()[7].split()[1:]))
        vol_redcell = cell_volum(reduce_cell)
        
        print_o2("Selected reduced cell axes", fmt_cell % reduce_cell[:3],f3)
        print_o2("Selected reduced cell angles", fmt_cell % reduce_cell[3:],f3)
        print_o2("Fraction of indexed spots","(%d/%d) %7.1f%%" % \
                    (num_indexed, num_total,100.*num_indexed/num_total),f3)
        print_o2("Standard deviation of spot    position (pixels) ",\
                                                       "%8.2f" % stdev_xy,f3)
        print_o2("Standard deviation of spindel position (degrees)",\
                                                       "%8.2f" % stdev_z,f3)
        print_o("  \n\n  Selected list of possible Bravais lattices\n",f3)
        print_o(("\n  Bravais Lattice Symm   fit"+15*" "+"cell"+17*" "+
                "ORDER\n\t"+26*" "+"a     b     c  alp  bet  gam\n\n"),f3)
                
        for LAT in select_lattices(limit=100):
             rc_au_vol, au_vol = "", ""
             a,b,c,al,be,ga = fix_cell(LAT[2],LAT[3])
             vol = cell_volum((a,b,c,al,be,ga))
             brav = BravaisDico[LAT[2][0]]
             mult = SPGlib[Bravais_to_Laue[LAT[2]][0][3][0]][2]
             if LAT[0] <= _idxref_fit_limit:
                  au_vol = mult/(vol/cell_volum(reduce_cell))
                  rc_au_vol = int(round(au_vol))*"*"
                  au_vol = "%.0f" % au_vol
             fmt_l = "  %1s %-13s (%2s) %5.1f "+3*"%6.0f"+3*"%5.0f"+"%3s  %s\n"
             print_o(fmt_l % (LAT[2][1],brav,LAT[2],LAT[0],\
                              a,b,c,al,be,ga,au_vol,rc_au_vol),f3)
        print_o("\n  ORDER is defined as the ratio: reduce cell volum/asym. unit volum.\n", f3)
        
    if step <= 3:
        Stage = " III - Integrate and Optimize Profiles"
        #==============================================#
        print_stage(fmt_stg % Stage,f3)
        
        xpar.JOB = "DEFPIX", "INTEGRATE"
        xpar.DATA_RANGE = collect_range
        reduce_cell = tuple(map(float,\
              opReadCl("XPARM.XDS").splitlines()[7].split()[1:]))
        print_o2("Selected reduced cell axes",fmt_cell % reduce_cell[:3],f3)
        print_o2("Selected reduced cell angles",fmt_cell % reduce_cell[3:],f3)
        print_o2("Runing integration with Space Group",SPGlib[1][1],f3)
        print_o2("Anomalous set","OFF",f3)
        print_o2("Number of optimization steps","%s" % _nopt,f3)
        
        for i in range(_nopt):
            seq = lookup_frame_sequences(auto.img_template)
            collect_range = choose_sequence([seq[0]],min_f=auto.frame_first,
                                    max_f=auto.frame_last)
            nframe = int(_job_angle[i]/auto.delta_phi)-1       
            collect_range[1] = min(collect_range[1],collect_range[0]+nframe)
            xpar.DATA_RANGE = collect_range
            print_o("\n",f3)
            print_o2("Step %d/%d  Integrate image range" % (i+1,_nopt),\
                                 "%d to %d" % tuple(xpar.DATA_RANGE),f3)
            xdsout = "xds_iii_%d.out"%i
            fhtml.write(fmt_link_xds % xdsout)
            run_xds(xpar, out_f=xdsout)
            xpar.JOB = "INTEGRATE"
            if i == 0:
                ipar = getProfilRefPar(infile="INTEGRATE.LP",init="yes")
                print_o2("Initial fitted "+fmt_shape, "%.3f / %.3f" % \
                            (ipar.BEAM_DIVERGENCE,ipar.REFLECTING_RANGE),f3)
            xpar.mix(getProfilRefPar())
            print_o2("New refined    "+fmt_shape, "%.3f / %.3f" % \
                            (xpar.BEAM_DIVERGENCE,xpar.REFLECTING_RANGE),f3)
        xpar.REFINE_INTEGRATE = "ORIENTATION", "BEAM",  "AXIS" #"DISTANCE",

    if step <= 4:
        Stage = " IV -  Testing Laue Symmetries"
        #=====================================================#
        print_stage(fmt_stg % Stage,f3)

        selected = select_lattices()
        print_o(("\n Lat  fit"+13*" "+"cell"+15*" "+
                "Symm  Compar Compl  Rsym Rmea3 Misf ORDER\n"+15*" "+
                "a    b    c  al  be  ga"+15*" "+3*"   (%)"+"\n\n"),f3)

        stage4, listRmeas3, BestOrder = [], [], []
        reduce_cell = map(float,opReadCl("XPARM.XDS").splitlines()[7].split()[1:])
        vol_redcell = cell_volum(reduce_cell)
        
        for _latt in selected:
             newcell = fix_cell(_latt[2],_latt[3])
             print_o(("%3s %5.1f" % (_latt[2],_latt[0]),\
                      " %5.0f%5.0f%5.0f%4.0f%4.0f%4.0f " % newcell),f3)
             vol_cell = cell_volum(newcell)
             
             n = 0
             for laue in Bravais_to_Laue[_latt[2]]:
                 mult = SPGlib[laue[1]][2]
                 au_fraction = int(round(mult/(vol_cell/vol_redcell)))
                 if n != 0: print_o(39*" ",f3)
                 print_o("%6s" % (laue[2]),f3)
                 
                 xpar.JOB = "CORRECT"
                 xpar.SPACE_GROUP_NUMBER = laue[1]
                 xpar.UNIT_CELL_CONSTANTS = newcell
                 xpar.REIDX = _latt[4]
                 run_xds(xpar, out_f="xds_iv.out")
                 
                 s = resum_scaling()
                 if not s:
                     print_o(fmt_err % "in CORRECT step",f3)
                     continue
                 print_o(" %7d" % (s.compar),f3)
                 for k in ['compl','rsym','rmeas3']:
                     print_o( "%6.1f" % (s[k]),f3)
                 print_o("%5d   %d\n" % (s.misfit, au_fraction),f3)
                 shutil.copyfile("CORRECT.LP","CORRECT.LP"+"_%s"%laue[2])
                 stage4.append([newcell,_latt[2],laue,s])
                 BestOrder.append(s.compl)
                 n = 1
                 if s.compar3 > 0:
                     listRmeas3.append([s.rmeas3,s.compar3])

        print_o("\n  ORDER is defined as the ratio: reduce cell volum/asym. unit volum.\n", f3)
        wMeanRmea3, wRMSRmea3 = wMean(listRmeas3), wStandardDeviation(listRmeas3)
        print_o("\n            Mean and RMS for Rmeas3 =  %7.2f  %7.2f\n" % \
                           (mean(listRmeas3), standardDeviation(listRmeas3)),f3)
        print_o("   Weighted Mean and RMS for Rmeas3 =  %7.2f  %7.2f\n" % \
                           (wMeanRmea3, wRMSRmea3),f3)
        
        i = 0
        print_o("\n  Chosen Lattice by order of preference:\n\n",f3)
        if len(listRmeas3) >= 3:
            BestOrder.sort()
            for a in range(len(stage4)):
                for L in stage4:
                    if L[-1].compl == BestOrder[-1] \
                    and 0 < L[-1].rmeas3 <= wMeanRmea3+wRMSRmea3*_wcut_rmea3 \
                    and i <= 3:
                       if not auto.spg: auto.spg = L[2][1]
                       cell = "%7.1f%7.1f%7.1f%6.1f%6.1f%6.1f" % L[0]
                       if i == 0: print_o("<b>",(fhtml,))
                       print_o("  %d) %s %5s %6s %4d" % \
                                (i+1,cell,L[2][0],L[2][2],L[2][1]),f3)
                       print_o("%5.1f %5.1f\n" % (L[3].compl,L[3].rmeas3),f3)
                       if i == 0: print_o("</b>",(fhtml,))
                       i += 1
                       del BestOrder[-1]
            
        
    if step <= 5:
        Stage = " V -  Integration"
        #=========================#
        print_stage(fmt_stg % Stage,f3)

        if xpar.has_key("REIDX"): del xpar.REIDX
        if auto.spg:
            BRAVAIS = ""
            spgn = 0
            if type(auto.spg) == type(""):
                auto.spg = auto.spg.upper()
                for _spg in SPGlib.keys():
                    if auto.spg == SPGlib[_spg][0]:
                        spgn = _spg
            elif type(auto.spg) == type(1):
                spgn = auto.spg
                if spgn not in SPGlib.keys():
                    print_o(fmt_err % \
                    "Space group not possible for biological macromolecules\n",f3)
                    stop()
            if not spgn:
                print_o(fmt_err % ("Unconsistant spage group: %s\n"%auto.spg),f3)
                stop()
            SPG = SPGlib[spgn]
        if not BRAVAIS:
            for bravais in Bravais_to_Laue.keys():
                for laue in Bravais_to_Laue[bravais]:
                    #print laue
                    #print spgn, laue[3]
                    if spgn in laue[3]:
                        BRAVAIS = bravais
                        #print " hit!",bravais
        
        if BRAVAIS == "mI": BRAVAIS = "mC"
        
        selected1 = select_lattices(limit=_idxref_fit_limit)
        selected2 = []
        bestfit = 1000
        #print "selected1", selected1
        #print "Bravais = ",BRAVAIS
        for _latt in selected1:
            if _latt[2] == BRAVAIS:
                if _latt[0] <= bestfit: bestfit = _latt[0]
                selected2.append(_latt)
        nsel = len(selected2)
        if nsel == 0:
            print_o("  No lattice could be selected for Space Group %s\n" % \
                                                SPG[1],f3)
            stop()
        for _latt in selected2:
            if _latt[0] == bestfit: selected = _latt
        
        newcell = fix_cell(selected[2],selected[3])
        xpar.UNIT_CELL_CONSTANTS = newcell
        xpar.SPACE_GROUP_NUMBER = spgn
        if auto.anomalous == "yes":
            xpar.FRIEDELS_LAW = "FALSE"
            ano = "ON"
        else:
            xpar.FRIEDELS_LAW = "TRUE"
            ano = "OFF"
        seq = lookup_frame_sequences(auto.img_template)
        xpar.DATA_RANGE = choose_sequence(seq,min_f=auto.frame_first,\
                                              max_f=auto.frame_last)
        
        print_o2("Selected cell axes",fmt_cell % newcell[:3],f3)
        print_o2("Selected cell angles",fmt_cell % newcell[3:],f3)
        print_o2("Selected Bravais lattice",BRAVAIS,f3)
        print_o2("Runing integration with Space Group", SPG[1],f3)
        print_o2("","(number %d)" % spgn,f3)
        print_o2("Using frames number","%d to %d" % tuple(xpar.DATA_RANGE),f3)
        print_o2("Anomalous set",ano,f3)
        
        xpar.mix(getProfilRefPar())
        print_o2("Using refined "+fmt_shape, "%.3f / %.3f" %\
                             (xpar.BEAM_DIVERGENCE, xpar.REFLECTING_RANGE),f3)
        
        xpar.JOB = "IDXREF"
        xdsout = "xds_v_1.out"
        fhtml.write(fmt_link_xds % xdsout)
        run_xds(xpar,inp_f=None,out_f=xdsout)
        
        xdsout = "xds_v_2.out"
        xpar.JOB = "DEFPIX", "INTEGRATE", "CORRECT"
        fhtml.write(fmt_link_xds % xdsout)
        run_xds(xpar,inp_f=None,out_f=xdsout)
        
        s = resum_scaling()
        if not s:
            print_o(fmt_err % "In CORRECT step.",f3)
            stop()

        print_o("\n\n\n\t\t\tScaling Statistics\n\n",f3)
        print_o3("Resolution","%.2f - %.2f" % (max(r),s.reso),\
                              "(%.2f - %.2f)" % (s.resoL,s.reso),f3)
        print_o3("Completeness","%.1f%%" % s.compl, "(%.1f%%)" % s.complL,f3)
        print_o3("Rmeas","%.1f%%" % s.rmeas, "(%.1f%%)" % s.rmeasL,f3)
        print_o3("Rsym","%.1f%%" % s.rsym, "(%.1f%%)" % s.rsymL,f3)
        print_o3("Compared","%d   " % s.compar, "(%d)" % s.comparL,f3)
        print_o3("I/sigma(I)","%.1f " % s.isig, "(%.1f)" % s.isigL,f3)
        print_o3("Multiplicity","%.1f " % (s.total/s.unique),"",f3)
        print_o3("Unique","%d   " % s.unique,"",f3)
        print_o3("Rejected misfits","%d   " % s.misfit,"",f3)
        print_o3("Systematic absents","%d   " % s.absent,"",f3)
        if s.absent:
            print_o3("   with <Iabs>/<I>","%.1f%%" % s.AbsIav,"",f3)
        if ano == "ON":
            print_o("\n",f3)
            print_o3("Anomalous contrib.","%.1f " % s.anom, "",f3)
        print_o("\n",f3)
        print_o3("Wilson scaling (B/Corr)","%.1f    " % s.wilson_b,\
                                         " %.2f" % s.wilson_corr,f3)
            
        print_o("\n\n",f3)
        xpar.mix(getProfilRefPar())
        print_o2("Final refined  "+fmt_shape, "%.3f / %.3f" % \
                             (xpar.BEAM_DIVERGENCE, xpar.REFLECTING_RANGE),f3)
        stop()
