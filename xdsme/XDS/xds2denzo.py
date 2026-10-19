#!/usr/bin/env python
"""
    10/06/03 legrand@embl-grenoble.fr
    28/06/03 last modifications.
    
    TODO:
        - detector conversion
        - finding the good film rotation for each detector...
"""

import Numeric
from math import sin, tan, atan, pi
import sys, os.path
from Scientific.Geometry.VectorModule import *
from Scientific.Geometry.TensorModule import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
            Going from XDS to Denzo:

  A program to convert the orientation matix, extract information
  from XDS output files and write denzo pseudo .x files:

    Usage:   %s   CORRECT.LP [-1|-1:10|-1,90,...]
    
  Image number selection: -n            (select image n)
                          -n1:n2        (select images n1 to n2)
                          -n1,n2,n3     (select image n1, n2 and n3)
        
""" % _progname

                    
mosaicity_conversion_factor = 4.410
# In XDS the mosaicity is discribed as:
# "the standard deviation (e.s.d) of a Gaussian modeling the rocking curve.
# If the limit of integration is set (by default) to 2% of reflection profile
# maximum, from the guassian function this correspond to a mosaic spread of
# 4.41 times the e.s.d.

class XP: pass

FMT = XP()
FMT.dotx_mat = 3*"%15.8f"+3*"%10.6f"+"\n"
FMT.dotx_par = 4*"%12.5f"+4*"%10.4f"+6*"%2d"+"\n"
FMT.refl = ("%4d%4d%4d" +
            " 0%8.1f   100.0   1.00%6.1f %5.3f" +
            "%7.1f%7.1f 1.000    99.9\n")

FMT.dotx_end = """ 999%(boxe_size_x)4d%(boxe_size_y)4d
%(spot_profile)s
spot elliptical %(spot_radius)5.3f %(spot_radius)5.3f    0.0
background elliptical %(background_radius)5.3f %(background_radius)5.3f    0.0
goniostat single axis
goniostat orientation    0.000    0.000
motor axis  0.000000  1.000000  0.000000
profile fitting radius   22.00
resolution limits  %(resolution_limits)s
wavelength  %(wavelength).5f
monochromator   %(monochromator).3f
oscillation start  %(phi_start).2f end  %(phi_end).2f
mosaicity   %(mosaicity).3f
spindle axis   0   0   1 vertical axis   1   0   0
unit cell    %(cell)s
crystal rotx    0.000 roty    0.000  rotz  0.000
sector        1
raw data file '%(image_template)s'
format %(detector)s
space group %(spg)s
cassette rotx   0.0 roty    0.0 rotz    %(cassette_rotz).2f 2 theta    0.00
distance %(distance_mm)7.2f
x beam %(beam_x)8.3f  y beam %(beam_y)8.3f
film rotation   %(film_rotation).2f
crossfire y  0.000 x  0.000 xy  0.000
"""

FMT.strategy_inp = """XFILE %(xfile_name)s
COMPLETENESS 100 99 98 97 96 95
NBINS 10
>NFRAME %(nframe)d
PSFILE %(psfile)s
FORMAT %(detector)s
SPACE GROUP %(spg)s
>RESOLUTION %(resolution_limits)s
%(anomalous)s
"""

SPGlib = {
# each SGnumber entry maps [SGsymbol1,SGsymbol2,SGorder]
1:['P1','P1',1],3:['P2','P2',2],4:['P21','P2(1)',2],5:['C2','C2',4],
5:['C2','C2',4],16:['P222','P222',4],17:['P2221','P222(1)',4],
18:['P21212','P2(1)2(1)2',4],19:['P212121','P2(1)2(1)2(1)',4],
21:['C222','C222',8],20:['C2221','C222(1)',8],22:['F222','F222',16],
23:['I222','I222',8],24:['I212121','I2(1)2(1)2(1)',8],
75:['P4','P4',4],76:['P41','P4(1)',4],77:['P42','P4(2)',4],
78:['P43','P4(3)',4],89:['P422','P422',8],90:['P4212','P42(1)2',8],
91:['P4122','P4(1)22',8],92:['P41212','P4(1)2(1)2',8],
93:['P4222','P4(2)22',8],94:['P42212','P4(2)2(1)2',8],
95:['P4322','P4(3)22',8],96:['P43212','P4(3)2(1)2',8],79:['I4','I4',8],
80:['I41','I4(1)',8],97:['I422','I422',16],98:['I4122','I4(1)22',16],
143:['P3','P3',3],144:['P31','P3(1)',3],145:['P32','P3(2)',3],
149:['P312','P312',6],150:['P321','P321',6],151:['P3112','P3(1)12',6],
152:['P3121','P3(1)21',6],153:['P3212','P3(2)12',6],
154:['P3221','P3(2)21',6],168:['P6','P6',6],169:['P61','P6(1)',6],
170:['P65','P6(5)',6],171:['P62','P6(2)',6],172:['P64','P6(4)',6],
173:['P63','P6(3)',6],177:['P622','P622',12],178:['P6122','P6(1)22',12],
179:['P6522','P6(5)22',12],180:['P6222','P6(2)22',12],
181:['P6422','P6(4)22',12],182:['P6322','P6(3)22',12],
146:['H3','R3',9],155:['H32','R32',18],195:['P23','P23',12],
198:['P213','P2(1)3',12],207:['P432','P432',24],208:['P4232','P4(2)32',24],
212:['P4332','P4(3)32',24],213:['P4132','P4(1)32',24],196:['F23','F23',48],
209:['F432','F432',96],210:['F4132','F4(1)32',96],197:['I23','I23',24],
199:['I213','I2(1)3',24],211:['I432','I432',48],214:['I4132','I4(1)32',48]}

detector2denzo = {"CCDCHESS":"ccd unsupported-mar","ADSC":"ccd adsc unsupported-q4","RAXIS":"RAXIS",
                    "MAR":"SMALLMAR","MAR345":"mar345 345mm 100 marpck","MAC":"DIP2000",
                    "CCDBRANDEIS":"ADSC"}

def detector2strategy(det_type, detector_diameter):
    if det_type == "CCDCHESS" and 166. > detector_diameter > 159.:
        return "MARCCD165"
    elif det_type == "CCDCHESS" and 135. > detector_diameter > 130.:
        return "MARCCD133"
    elif det_type == "ADSC" and 189. > detector_diameter > 187.:
        return "ADSCQ4R"
    elif det_type == "ADSC" and 211. > detector_diameter > 209.:
        return "ADSCQ210"
    elif det_type == "ADSC" and 316. > detector_diameter > 314.:
        return "ADSCQ315"
    elif det_type == "CCDBRANDEIS" and 201. > detector_diameter > 199.:
        return "BRAND2X2"
    elif det_type == "CCDBRANDEIS" and  96. > detector_diameter >  94.:
        return "CCDB"
    elif det_type == "MAC":
        return "DIP2"
    #"RAXIS":"RAXIS",
    #"MAR":"SMALLMAR","MAR345":"mar345 345mm 100 marpck"
                    
film_rotation_DB = {"MARCCD133":-90., "MARCCD165":0.,
                    "ADSCQ4R":90., "ADSCQ210":90.,
                    "ADSCQ315":90., "DIP2":-90.,}
                    
cassette_rotz_DB = {"MARCCD133":0., "MARCCD165":90.,
                    "ADSCQ4R":0., "ADSCQ210":0.,
                    "ADSCQ315":0., "DIP2":0.,}
                   

def spot_profile(boxe_size, background_radius, spot_radius):
    x0, y0 = boxe_size[0]/2.-0.5,boxe_size[1]/2.-0.5
    def cercle_2D(x, y):
        return Numeric.sqrt((x-x0)**2 + (y-y0)**2) <= radius
    radius = background_radius
    background_boxe = (Numeric.fromfunction(cercle_2D, boxe_size) == 0)
    radius = spot_radius
    spot_boxe = 2*Numeric.fromfunction(cercle_2D, boxe_size)
    boxe = spot_boxe+background_boxe
    boxe_str = ""
    line_fmt = " "+ boxe.shape[1]*"%1d"+ "\n"
    for l in boxe:
         boxe_str += line_fmt % tuple(l)
    return  boxe_str

def write_dotx(dnz, format, dotx_fname, hklxy_str):
    dotx = open(dotx_fname,"w")
    dotx.write("   Generated with XDS2DENZO\n")
    for l in range(3):
        dotx.write(format.dotx_mat % 
                 tuple(Numeric.ravel((dnz.Amat[l],dnz.Umat[l]))))
                 
    dotx.write(format.dotx_par % (dnz.phi_start, dnz.phi_end,
                   dnz.distance_pix,
                   dnz.wavelength, 0.0, 0.0, 0.0, 
                   dnz.mosaicity,
                   -1,1,-1,-1,-1,-1))
    dotx.write(hklxy_str)    
    dotx.write(format.dotx_end % vars(dnz))
    dotx.close()
    
def d_min(det_diameter, det_distance, wavelength):
    return 1/(2*sin(atan((det_diameter/2.)/det_distance)/2.)/wavelength)

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r
    
def getXpar(_str,match,limit=80,func=float):
    start = _str.index(match)+len(match)
    tmp = _str[start:start+limit].split("\n")[0].split()
    return map(func,tmp)

def read_correct_info(infname="CORRECT.LP"):
    xdsPar = {}
    corr = opReadCl(infname)
    corrp = corr[corr.index("PARAMETERS USING ALL IMAGES")+100:
           corr.index("DETERMINATION OF RECIPROCAL")-100]
    corri = corr[:1500]
    xdsPar["rot"] = getXpar(corrp,"ROTATION AXIS")
    xdsPar["beam"] = getXpar(corrp,"COORDINATES (REC. ANGSTROEM)")
    xdsPar["distance"] = getXpar(corrp,"DETECTOR DISTANCE (mm)")[0]
    xdsPar["origin"] = getXpar(corrp,"(PIXELS) OF DIRECT BEAM")
    xdsPar["a"] = getXpar(corrp,"CELL A-AXIS")
    xdsPar["b"] = getXpar(corrp,"CELL B-AXIS")
    xdsPar["c"] = getXpar(corrp,"CELL C-AXIS")
    xdsPar["cell"] = getXpar(corrp,"UNIT CELL PARAMETERS")
    xdsPar["mosaicity"] = getXpar(corrp,"CRYSTAL MOSAICITY (DEGREES)")[0]
    iqx, iqy  = corri.index("QX=")+3, corri.index("QY=")+3
    xdsPar["pixel_size"] = float(corri[iqx:iqx+10]),float(corri[iqy:iqy+10])
    itemp = corri.index("TEMPLATE_OF_DATA_FRAMES=")+24
    xdsPar["template"] = corri[itemp:itemp+80].split("\n")[0].split()[0].replace("?","#")
    xdsPar["symmetry"] = int(getXpar(corrp,"SPACE GROUP NUMBER")[0])
    xdsPar["detector"] = getXpar(corri,"DETECTOR=",40,str)[0]
    xdsPar["phi_init"] = getXpar(corri,"STARTING_ANGLE=",15)[0]
    xdsPar["num_init"] = getXpar(corri,"STARTING_FRAME=",15,int)[0]
    xdsPar["delta_phi"] = getXpar(corri,"OSCILLATION_RANGE=",15)[0]
    xdsPar["divergence_esd"] = getXpar(corri,"BEAM_DIVERGENCE_E.S.D.=",15)[0]
    xdsPar["polarization"] = getXpar(corri,"FRACTION_OF_POLARIZATION=",15)[0]
    xdsPar["resolution_range"] = getXpar(corri,"INCLUDE_RESOLUTION_RANGE=",20)
    xdsPar["npixel_x"] = getXpar(corri,"NX=",7,int)[0]
    xdsPar["npixel_y"] = getXpar(corri,"NY=",7,int)[0]
    xdsPar["friedel"] = getXpar(corri,"FRIEDEL'S_LAW=",7,str)[0]
    return xdsPar

def read_xds_hkl(hkl_str, ftype = "INTEGRATE"):
    hkl_str = hkl_str[hkl_str.index("!END_OF_HEADER")+15:
                      hkl_str.index("!END_OF_DATA")]
    hkl_num = Numeric.array(map(float,hkl_str.split()))
    _len = len(hkl_num)
    if ftype == "INTEGRATE":hkl_num.shape = (_len/19,19)
    if ftype == "XDS_ASCII":hkl_num.shape = (_len/11,11)
    return hkl_num

if __name__=='__main__':

    if len(sys.argv) == 1:
        print _usage
        sys.exit()
    image_num_list = 1,
    for arg in sys.argv:
        if arg[0] == "-":
            if arg.count(":") == 1:
                try:
                    II, IF = map(int,(arg[1:].split(":")))
                    image_num_list = range(II, IF+1)
                    sys.argv.remove(arg)
                except: pass
            elif arg.count(",") >= 1:
                try:
                    image_num_list = map(int,(arg[1:].split(",")))
                    sys.argv.remove(arg)
                except: pass
            else:
                try:
                    image_num_list = (float(arg[1:]),)
                    sys.argv.remove(arg)
                except: pass
    print ">>> Selected Image(s) = ", image_num_list
    
    xds_dir, infile_type = os.path.split(sys.argv[1])
    if infile_type[:10] == "CORRECT.LP":
        xdsPar = read_correct_info(sys.argv[1])
    else:
        print "\nERROR! Unrecognize commande: ",sys.argv[1:]
        print _usage
        sys.exit()

    if len(sys.argv) >= 2:
        
        #hkl_name = "XDS_ASCII.HKL"
        hkl_name = "INTEGRATE.HKL"
        if hkl_name == "INTEGRATE.HKL":
            iZD, irXYD, ncol = 7, 5, 18
        if hkl_name == "XDS_ASCII.HKL":
            iZD, irXYD, ncol = 7, 4, 10
        if xds_dir: hkl_name = xds_dir + "/" + hkl_name
        if os.path.exists(hkl_name) and os.path.isfile(hkl_name):
            print ">>> Reading reflection file ", hkl_name
            hkl = opReadCl(hkl_name)
            hkl_xds = read_xds_hkl(hkl,"INTEGRATE")
            D = xdsPar["distance"]/xdsPar["pixel_size"][0]
            XY = hkl_xds[:,irXYD:irXYD+2][:,::-1]
            X, Y = XY[:,0], XY[:,1]
            X0, Y0 = xdsPar["origin"][1], xdsPar["origin"][0]
            # OBL is the cosine of incidence angle at the detector
            # In this case the Detector is supposed perpendicular to the beam.
            OBL = D/Numeric.sqrt(D**2+(X-X0)**2+(Y-Y0)**2)
            OBL.shape = (OBL.shape[0],1)
            hklxy = Numeric.concatenate((hkl_xds[:,:5], OBL, XY),1)
            # it seems that denzo start counting pixels at [0,0] whereas XDS
            # start at [1,1] ???
            # hklxy[:,5:6] = hklxy[:,5:6] -1
        else:
            print ">>> Warning: Can't find reflection file: ", hkl_name
            hkl = ""
            hkl_xy = [[]]
            hklxy_str = ""


    # Transforming XDS orientation matrices into
    # denzo A and U matrices
    
    a = Vector(xdsPar["a"])
    b = Vector(xdsPar["b"])
    c = Vector(xdsPar["c"])
    beam = Vector(xdsPar["beam"])
    rot = Vector(xdsPar["rot"])

    volum = a.cross(b)*c
    ar = b.cross(c)/volum
    br = c.cross(a)/volum
    cr = a.cross(b)/volum
    wavelength = 1/beam.length()

    # cammat: Creat orthonormal axes with
    # y // to the rotation (spindel) axis
    # z // to the beam
    # x perpendicular to z and to the beam
    if 1:
        cammat_Y = rot.normal()
        cammat_X = cammat_Y.cross(beam).normal()
        cammat_Z = cammat_X.cross(cammat_Y)
        #print "Signes of cammat:",cammat_X.cross(cammat_Y)*cammat_Z
   
    cammat = Tensor([cammat_X,cammat_Y,cammat_Z])

    Amat = Numeric.transpose([cammat*ar,cammat*br,cammat*cr])*wavelength
    #for l in Amat:  print 3*"%15.8f" % tuple(l)
    
    tU_z = c.normal()
    tU_x = ar.normal()
    tU_y = tU_z.cross(tU_x)

    cammat = Tensor([cammat_Z,cammat_X,cammat_Y])
    Umat =  Numeric.transpose([cammat*tU_x,cammat*tU_y,cammat*tU_z])
    #for l in Umat:  print 3*"%10.6f" % tuple(l)


    # Preparing Denzo parameters.
    dnz = XP()
    dnz.Umat = Umat
    dnz.Amat = Amat
    dnz.wavelength = wavelength
    dnz.delta_phi = xdsPar["delta_phi"]
    dnz.pixel_size_x = xdsPar["pixel_size"][0]
    dnz.pixel_size_y = xdsPar["pixel_size"][1]
    dnz.distance_pix = xdsPar["distance"]/dnz.pixel_size_x
    dnz.distance_mm = xdsPar["distance"]
    dnz.phi_init = xdsPar["phi_init"]
    dnz.delta_phi = xdsPar["delta_phi"]
    dnz.mosaicity = mosaicity_conversion_factor*xdsPar["mosaicity"]
    dnz.beam_x = xdsPar["origin"][1]*xdsPar["pixel_size"][1]
    dnz.beam_y = xdsPar["origin"][0]*xdsPar["pixel_size"][0]
    dnz.image_template = xdsPar["template"]
    dnz.spg = SPGlib[xdsPar["symmetry"]][0].lower()
    dnz.monochromator = xdsPar["polarization"]
    dnz.spot_radius = dnz.distance_mm*tan(2.5*xdsPar["divergence_esd"]*pi/180)
    dnz.diameter = min(xdsPar["npixel_x"]*dnz.pixel_size_x,
                        xdsPar["npixel_y"]*dnz.pixel_size_y)
    dnz.cell = 6*"%.3f  " % tuple(xdsPar["cell"])
    dnz.detector = detector2denzo[xdsPar["detector"]]
    resolution_limits = tuple(xdsPar["resolution_range"])
    if min(resolution_limits) == 0.:
        resolution_limits = (max(resolution_limits), 
               d_min(dnz.diameter, dnz.distance_mm, dnz.wavelength))
    dnz.resolution_limits = "%.2f  %.2f" % resolution_limits 
    _detector = detector2strategy(xdsPar["detector"],dnz.diameter)
    dnz.film_rotation = film_rotation_DB[_detector]
    dnz.cassette_rotz = cassette_rotz_DB[_detector]
    

    # Determining denzo spot profile
    spot_radius = int(dnz.spot_radius/(dnz.pixel_size_x)+0.5)
    bkg_radius  = spot_radius+1
    dnz.boxe_size_x = (spot_radius+3)*2
    dnz.boxe_size_y = (spot_radius+3)*2
    dnz.background_radius = bkg_radius*dnz.pixel_size_x
    dnz.spot_profile = spot_profile((dnz.boxe_size_x,dnz.boxe_size_y),
                                     bkg_radius, spot_radius)[:-1]
    

    # Selecting reflections and Writing dotx files
    prefix = dnz.image_template.split("/")[-1].split(".")[0].replace("#","")
    dotx_fname_fmt = prefix + "%.3d.x.xds"
    for image_num in image_num_list:
         dotx_fname = dotx_fname_fmt % image_num
         dnz.phi_end = dnz.phi_init + image_num*dnz.delta_phi
         dnz.phi_start = dnz.phi_end - dnz.delta_phi
         frame_filter = (image_num > hkl_xds[:,iZD])*\
                           (hkl_xds[:,iZD] >= (image_num-1))
         hklxy_select = Numeric.compress(frame_filter,hklxy,0)      
         nref = hklxy_select.shape[0]
         out_fmt1 = ">>> Number of reflection selected for frame %4d : %d"
         print  out_fmt1 % (image_num,nref)
         hklxy_str = nref*FMT.refl  % tuple(Numeric.ravel(hklxy_select))
         write_dotx(dnz, FMT, dotx_fname, hklxy_str)
    

    # Preparing STRATEGY input file.
    STR = XP()
    str_inpfname = "auto_strategy.inp"
    STR.xfile_name = dotx_fname_fmt % image_num_list[0]
    STR.psfile = str_inpfname.split(".")[0]+".ps"
    STR.spg = dnz.spg.upper()
    STR.nframe = 90
    STR.resolution_limits = dnz.resolution_limits
    STR.detector = _detector
    if xdsPar["friedel"] == "FALSE":
        STR.anomalous = "ANOMALOUS"
    else: 
        STR.anomalous = ""
    strat = open("auto_strategy.inp","w")
    strat.write(FMT.strategy_inp % vars(STR))
    strat.close()
