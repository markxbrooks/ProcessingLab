#!/usr/bin/env python
"""
    15/01/03 First version
    xds2mos.py is inspired form the wk2al.f
    program by Richard Kahn (IBS-Grenoble)
    
    Convert: Orientation matrix, Beam center, detector type
    Try to guess: the Omega and 2theta angles.
    
    TODO:
        - Guess the 2theta and add it to the input.

    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
    
    Uses the ScientificPython module by Konrad Hinsen
    http://starship.python.net/crew/hinsen/scientific.html
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "19-10-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "GPL"
__version__ = "0.4.5"


import sys
import os
import math

#try:
#    from cgtypes import vec3 #as Vector
#    from cgtypes import mat3 #as Tensor
#except ImportError:
#    try:
#        from pycgtypes import vec3 #as Vector
#        from pycgtypes import mat3 #as Tensor
#    except ImportError:
#        print "Can't find cgtypes or pycgtype modules... Please check your"
#        print "$PYTHONPATH variable if it is installed somewhere or instal" 
#        print "one of these module for vector and matrix operations."
#        sys.exit()

from ThreeAxisRotation import *
from CrystalOrientation import *

_debug = False
_progname = os.path.split(sys.argv[0])[1]
_usage = """
            Going from XDS to Mosflm:

  A program to convert the orientation matix, extract information
  from XDS output files and write a mosflm input file:

    Usage:   %s  [-m][-Angles] XPARM.XDS|GXPARM.XDS|IDXREF.LP|CORRECT.LP
    
    -m :       Start "ipmosflm < xds2mos.inp" to verify predictions
    -Angles:   Writes out the crystal orientation in xds2mos.umat 
               as setings angles (default is as a U matrix)
        
""" % _progname

detector2scanner = {
   "ADSC":              "ADSC",
   "CCDCHESS":          "MARCCD",
   "RAXIS":             "RAXIS",
   "MAR":               "SMALLMAR",
   "MAR345":            "MAR",
   "MAC":               "DIP2020",
   "SMARTCCD":          "SMART",
   "CCDD2AM":           "ESRF",
   "CCDBRANDEIS":       "B4"}
                    

# Depending on the XDS output read, different mosflm input can be generated.
mosflmInpTemplate1 = """
TILE      Input created from  %(title)s

DEBUG CONTROL

MATRIX            xds2mos.umat

! The DETECTOR type should be detected automaticaly by Mosflm,
! In case of problem, try to comment/uncomment the next line.   
DETECTOR         %(detector)s

! If the image orientation doesn't match the predictions try to
! comment/uncomment the DETECTOR line and play with the OMEGA values
! (0, 90, 180, 270)
! DETECTOR OMEGA  %(omega)12.2f

! TWOTHETA          %(twotheta)12.2f  # This is experimental...

WAVELENGTH        %(wavelength)12.5f
DISTANCE          %(distance)12.3f
BEAM              %(beam_x)12.3f %(beam_y)12.3f
PIXEL             %(pixel_x)12.5f %(pixel_y)12.5f
SYMMETRY          %(symmetry)s
TEMPLATE          %(template)s
EXTENTION         %(extention)s

IMAGE             %(image_numb)d PHI %(phi_i)9.3f TO %(phi_f)9.3f
GO
"""
mosflmInpTemplate2 = """
MOSAICITY        %(mosaicity).3f
"""

mosaicity_conversion_factor = 4.410
# In XDS the mosaicity is discribed as:
# "the standard deviation (e.s.d) of a Gaussian modeling the rocking curve.
# If the limit of integration is set (by default) to 2% of reflection profile
# maximum, from the guassian function this correspond to a mosaic spread of
# 4.41 times the e.s.d.

spg_num2symb = {1:'P1',3:'P2',4:'P21',5:'C2',16:'P222',17:'P2221',18:'P21212',
19:'P212121',20:'C2221',21:'C222',22:'F222',23:'I222',24:'I212121',75:'P4',
76:'P41',77:'P42',78:'P43',79:'I4',80:'I41',89:'P422',90:'P4212',91:'P4122',
92:'P41212',93:'P4222',94:'P42212',95:'P4322',96:'P43212',97:'I422',98:'I4122',
143:'P3',144:'P31',145:'P32',146:'H3',149:'P312',150:'P321',151:'P3112',
152:'P3121',153:'P3212',154:'P3221',155:'R32',168:'P6',169:'P61',170:'P65',
171:'P62',172:'P64',173:'P63',177:'P622',178:'P6122',179:'P6522',180:'P6222',
181:'P6422',182:'P6322',195:'P23',196:'F23',197:'I23',198:'P213',199:'I213',
207:'P432',208:'P4232',209:'F432',210:'F4132',211:'I432',212:'P4332',
213:'P4132',214:'I4132'}


def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r
    
def str_mat(mat, name="", format="%12.6f"):
    "Return a string representation of a 3x3 matrix"
    Pformat = name + 3*(3*format+"\n")
    return Pformat % tuple(mat.toList(1))

def getOmega(xdsPar):
    """Calculate an Omega value (in radian) wich defines how the fast (X) and
    slow (Y) axis of detector files are  orientated toward the camera frame.
    The calculation of this omega value is supposed to reflect the Mosflm
    definition...
    
    But it seems that I get different values from the mosflm defaults... This
    may be due to: A) My missanderstanding of the mosflm documentation, B)
    Some tricks in the image reading routines.
    
    Nonetheless, this calculated value works for translating
    correctly the beam coordinates from XDS to mosflm [at least in the tested
    cases of MARCCD, MAR345 and ADSC detector images].
    
    Reference:
    http://www.ccp4.ac.uk/dist/x-windows/Mosflm/doc/mosflm_user_guide.html#a3
    """
    
    # Xd = CAMERA_y = beam
    # Yd = CAMERA_z = rot
    Xd =  vec3(xdsPar["beam"]).normalize()
    Yd =  vec3(xdsPar["rot"]).normalize()
    CAMERA_x = Xd.cross(Yd)
    CAMERA = mat3(CAMERA_x, Xd, Yd).transpose()
    
    # This is the definition of the fast:X and slow:Y axis for the detector files.
    XDSdetector_X = vec3(xdsPar["detector_X"])
    XDSdetector_Y = vec3(xdsPar["detector_Y"])
            
    # Now this axes are translated in the mosflm Camera frame
    Xs = XDSdetector_X*CAMERA
    Ys = XDSdetector_Y*CAMERA
    
    # Both angles should be identical.
    omegaX = Xd.angle(Xs)
    omegaY = Yd.angle(Ys)
    
    if _debug:
        print "X xds: fast =",XDSdetector_X
        print "Y xds: slow =",XDSdetector_Y
        print "Xs: fast =",XDSdetector_X,"->", Xs
        print "Ys: slow =",XDSdetector_Y,"->", Ys
        print "Xd: ", Xd
        print "Yd: ", Yd
    
        print 
        print "Detector: %s" % detector2scanner[xdsPar["detector_type"]]
        print "OmegaX:   %8.2f" % (omegaX*r2d)
        print "OmegaY:   %8.2f" % (omegaY*r2d)
        
    return omegaX

def getTwoTheta(xdsPar):
    """Tries to calculate the 2theta angle (in radian) of the detector.
    I am not completely sure of this calculation. How 2theta is precisely
    geometricaly defined in mosflm ?
    I need to look in the mosflm code where it is taken into account.
    """
    
    beam = vec3(xdsPar["beam"])
    rot  = vec3(xdsPar["rot"]).normalize()
    camY = rot.cross(beam)
    
    XDSdetector_X = vec3(xdsPar["detector_X"]).normalize()
    XDSdetector_Y = vec3(xdsPar["detector_Y"]).normalize()
    #XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)
    
    #print beam.angle(XDSdetector_Z)*r2d
    if abs(rot * XDSdetector_X) - 1 <= 0.05:
        detecorVector = -XDSdetector_Y
        #print 1
    elif abs(rot * XDSdetector_Y) - 1 <= 0.05:
        detecorVector = XDSdetector_X
        #print 2
    else:
        raise Exception, "Can't calculate TwoTheta angle"    
    
    return camY.angle(detecorVector)
    

def getBeamOrigin(xdsPar):
    """Calculate the direct beam coordinates on the detector."""
    
    distance = xdsPar["distance"]
    beam = vec3(xdsPar["beam"])
    
    XDSdetector_X = vec3(xdsPar["detector_X"])
    XDSdetector_Y = vec3(xdsPar["detector_Y"])
    XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)
    
    # Calculate the direct beam coordinates on the detector
    beamOx = xdsPar["originXDS"][0]*xdsPar["pixel_size"][0]
    beamOy = xdsPar["originXDS"][1]*xdsPar["pixel_size"][1]
    beam0z = beam*XDSdetector_Z
    
    beamX = beamOx + beam*XDSdetector_X*distance/beam0z
    beamY = beamOy + beam*XDSdetector_Y*distance/beam0z
    
    beamXp = beamX/xdsPar["pixel_size"][0]
    beamYp = beamY/xdsPar["pixel_size"][1]
    
    if _debug:
        if "origin" in xdsPar:
            print "\nBEAM center read from XDS in pixel:",
            print  "%9.2f %9.2f" % tuple(xdsPar["origin"])
            # When given by XDS, verifies that my calculation is correct
            assert (xdsPar["origin"][0] - beamXp) < 0.05
            assert (xdsPar["origin"][1] - beamYp) < 0.05
        print "BEAM center calculated in pixel:\t%9.2f %9.2f" % (beamXp,beamYp)
        print "BEAM center calculated in mm:\t\t%9.2f %9.2f\n" % (beamX,beamY)
    
    return beamX, beamY
    

def OM_xds2mos(xdsPar):
    """ Convert the XDS direct space Orientation Matrix to a mosflm OM
    
    Mosflm CAMERA coordinate frame has orthonormal axes with:
      z // rotation axis
      y perpendicular to z and to the beam
      x perpendicular to y and z (along the beam)
    
    For more details see the mosflm documentation:
    http://www.ccp4.ac.uk/dist/x-windows/Mosflm/doc/mosflm_user_guide.html#a3
    """
    
    a = vec3(xdsPar["a"])
    b = vec3(xdsPar["b"])
    c = vec3(xdsPar["c"])
    beam = vec3(xdsPar["beam"])
    rot = vec3(xdsPar["rot"])
    wavelength = 1/beam.length()
    
    volum = a.cross(b)*c
    ar = b.cross(c)/volum
    br = c.cross(a)/volum
    cr = a.cross(b)/volum
    
    CAMERA_z = rot.normalize()
    CAMERA_y = CAMERA_z.cross(beam).normalize()
    CAMERA_x = CAMERA_y.cross(CAMERA_z)
    CAMERA = mat3(CAMERA_x,CAMERA_y,CAMERA_z).transpose()
    
    UBxds = mat3(ar,br,cr)
    UBmos = CAMERA*UBxds*wavelength
        
    return UBmos

def getXpar(_str,match,limit=70,func=float):
    start = _str.index(match)+len(match)
    tmp = _str[start:start+limit].splitlines()[0].split()
    return map(func,tmp)

def write_mosflm_matrices(filename, UB, U, phi123, cell):
    "Writes a mosflm formated crystal oritentation matrice file"
    
    matf = open(filename,"w")
    fmt_cell = 6 * "%12.4f" + "\n"
    fmt_angles = 3 * "%12.3f" + "\n"
    matf.write(str_mat(UB, name="", format="%12.8f"))
    matf.write(fmt_angles % (0,0,0))
    matf.write(str_mat(U, name="", format="%12.8f"))
    matf.write(fmt_cell % tuple(cell))
    matf.write(fmt_angles % tuple(phi123))
    matf.close()

def read_xdsinp_info(infname="XDS.INP"):
    xdsPar = {}
    xdsinp = opReadCl(infname)
    xdsPar["detector_type"] = getXpar(xdsinp,"DETECTOR=",40,str)[0]
    it = xdsinp.index("TEMPLATE_OF_DATA_FRAMES=")+24
    xdsPar["template"] = xdsinp[it:it+60].split()[0].replace("?","#")
    return xdsPar

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
    xdsPar["originXDS"] = getXpar(corrp,"ORIGIN (PIXELS) AT ")
    xdsPar["a"] = getXpar(corrp,"CELL A-AXIS")
    xdsPar["b"] = getXpar(corrp,"CELL B-AXIS")
    xdsPar["c"] = getXpar(corrp,"CELL C-AXIS")
    xdsPar["cell"] = getXpar(corrp,"UNIT CELL PARAMETERS")
    xdsPar["mosaicity"] = getXpar(corrp,"CRYSTAL MOSAICITY (DEGREES)")[0]
    iqx, iqy  = corri.index("QX=")+3, corri.index("QY=")+3
    xdsPar["pixel_size"] = float(corri[iqx:iqx+9]),float(corri[iqy:iqy+9])
    inx, iny  = corri.index("NX=")+3, corri.index("NY=")+3
    xdsPar["pixel_numb"] = int(corri[inx:inx+8]),int(corri[iny:iny+8])
    it = corri.index("TEMPLATE_OF_DATA_FRAMES=")+24
    xdsPar["template"] = corri[it:it+60].split()[0].replace("?","#")
    xdsPar["symmetry"] = int(getXpar(corrp,"SPACE GROUP NUMBER")[0])
    xdsPar["detector_type"] = getXpar(corri,"DETECTOR=",40,str)[0]
    xdsPar["detector_X"] = getXpar(corri,"DETECTOR_X-AXIS=")
    xdsPar["detector_Y"] = getXpar(corri,"DETECTOR_Y-AXIS=")
    xdsPar["phi_init"] = getXpar(corri,"STARTING_ANGLE=",15)[0]
    xdsPar["num_init"] = getXpar(corri,"STARTING_FRAME=",15,int)[0]
    xdsPar["delta_phi"] = getXpar(corri,"OSCILLATION_RANGE=",15)[0]
    return xdsPar

def read_idxref_info(infname="IDXREF.LP"):
    xdsPar = {}
    idxr = opReadCl(infname)
    idxrp = idxr[idxr.index("FRACTION PARAMETERS USED AT START OF INTEGRAT"):]
    idxri = idxr[:1200]
    xdsPar["rot"] = getXpar(idxrp,"ROTATION AXIS")
    xdsPar["beam"] = getXpar(idxrp,"COORDINATES (REC. ANGSTROEM)")
    xdsPar["distance"] = getXpar(idxrp,"DETECTOR DISTANCE (mm)")[0]
    xdsPar["origin"] = getXpar(idxrp,"(PIXELS) OF DIRECT BEAM")
    xdsPar["originXDS"] = getXpar(idxrp,"ORIGIN (PIXELS) AT ")
    xdsPar["a"] = getXpar(idxrp,"CELL A-AXIS")
    xdsPar["b"] = getXpar(idxrp,"CELL B-AXIS")
    xdsPar["c"] = getXpar(idxrp,"CELL C-AXIS")
    xdsPar["cell"] = getXpar(idxrp,"UNIT CELL PARAMETERS")
    iqx, iqy  = idxri.index("QX=")+3, idxri.index("QY=")+3
    xdsPar["pixel_size"] = float(idxri[iqx:iqx+9]),float(idxri[iqy:iqy+9])
    inx, iny  = idxri.index("NX=")+3, idxri.index("NY=")+3
    xdsPar["pixel_numb"] = int(idxri[inx:inx+8]),int(idxri[iny:iny+8])
    it = idxri.index("TEMPLATE_OF_DATA_FRAMES=")+24
    xdsPar["template"] = idxri[it:it+60].split()[0].replace("?","#")
    xdsPar["symmetry"] = int(getXpar(idxrp,"SPACE GROUP NUMBER")[0])
    xdsPar["detector_X"] = getXpar(idxri,"DETECTOR_X-AXIS=")
    xdsPar["detector_Y"] = getXpar(idxri,"DETECTOR_Y-AXIS=")
    xdsPar["phi_init"] = getXpar(idxri,"STARTING_ANGLE=",15)[0]
    xdsPar["num_init"] = getXpar(idxri,"STARTING_FRAME=",15,int)[0]
    xdsPar["delta_phi"] = getXpar(idxri,"OSCILLATION_RANGE=",15)[0]
    return xdsPar

def read_xparm_info(infname="XPARM.XDS"):
    xdsPar = {}
    xparm = opReadCl(infname).splitlines()
    if len(xparm[0].split()) == 3:
        xdsPar["rot"] = map(float,xparm[0].split())
    if len(xparm[0].split()) == 6:
        xdsPar["rot"] = map(float,xparm[0].split()[3:])
    xdsPar["beam"] = map(float,xparm[1].split()[1:])
    xdsPar["distance"] = float(xparm[3].split()[0])
    xdsPar["originXDS"] = map(float,xparm[3].split()[1:])
    xdsPar["a"] = map(float,xparm[8].split())
    xdsPar["b"] = map(float,xparm[9].split())
    xdsPar["c"] = map(float,xparm[10].split())
    xdsPar["cell"] = map(float,xparm[7].split()[1:])
    xdsPar["pixel_size"] = map(float,xparm[2].split()[2:])
    xdsPar["pixel_numb"] = map(float,xparm[2].split()[:2])
    xdsPar["symmetry"] = int(xparm[7].split()[0])
    xdsPar["num_init"], xdsPar["phi_init"], xdsPar["delta_phi"] = \
                map(float,xparm[0].split()[:3])
    xdsPar["detector_X"] = map(float,xparm[4].split())
    xdsPar["detector_Y"] = map(float,xparm[5].split())
    return xdsPar

def PARS_xds2best(xdsPar, _read_idxref=0,_read_correct=0):
    "Convert XDS output parameters to Mosflm input parameters."
    
    mosPar = {}
    mosPar["title"] = "xds2mos   version: %s" % (__version__)
    mosPar["distance"] = abs(xdsPar["distance"])
    mosPar["wavelength"] = 1/vec3(xdsPar["beam"]).length()
    mosPar["symmetry"] = spg_num2symb[xdsPar["symmetry"]]
    mosPar["omega"] = xdsPar["omega"]*r2d
    mosPar["twotheta"] = xdsPar["twotheta"]*r2d
        
    xc = xdsPar["origin"][0]
    yc = xdsPar["origin"][1]
    
    cosOmega = math.cos(xdsPar["omega"])
    sinOmega = math.sin(xdsPar["omega"])
    
    mosPar["beam_x"] = xc*cosOmega + yc*sinOmega
    mosPar["beam_y"] = xc*sinOmega + yc*cosOmega
    
    mosPar["detector"] = detector2scanner[xdsPar["detector_type"]]
    mosPar["pixel_x"] =  xdsPar["pixel_size"][1]
    mosPar["pixel_y"] =  xdsPar["pixel_size"][0]
    mosPar["template"] = xdsPar["template"]
    mosPar["extention"] = xdsPar["template"].split(".")[-1]
    mosPar["image_numb"] = xdsPar["num_init"]
    mosPar["phi_i"] = xdsPar["phi_init"]
    mosPar["phi_f"] = xdsPar["phi_init"] + xdsPar["delta_phi"]
    
    if _read_correct:
        mosPar["mosaicity"] = mosaicity_conversion_factor*xdsPar["mosaicity"]
        
    return mosPar

if __name__=='__main__':
    _start_mosflm = 0
    _write_out_angles = 0
    _read_xparm = 0
    _read_idxref = 0
    _read_correct = 0
    
    if len(sys.argv) == 1:
        print _usage
        sys.exit()
        
    if len(sys.argv) >= 2:
        if "-m" in sys.argv:
            _start_mosflm = 1
            sys.argv.remove("-m")
        if "-Angles" in sys.argv:
            _write_out_angles = 1
            sys.argv.remove("-Angles")
    
    infilelocation, infilename = os.path.split(sys.argv[1])
    
    if infilename[:10] == "GXPARM.XDS" or infilename[:9] == "XPARM.XDS":
        xdsPar = read_xparm_info(sys.argv[1])
        xdsinpFile = os.path.join(infilelocation, "XDS.INP")
        xdsPar.update(read_xdsinp_info(xdsinpFile))
        _read_xparm = 1
        
    elif infilename[:9] == "IDXREF.LP":
        xdsPar = read_idxref_info(sys.argv[1])
        xdsinpFile = os.path.join(infilelocation, "XDS.INP")
        xdsPar.update(read_xdsinp_info(xdsinpFile))
        _read_idxref = 1
    
    elif infilename[:10] == "CORRECT.LP":
        xdsPar = read_correct_info(sys.argv[1])
        _read_correct = 1
    
    else:
        print "\nERROR! Unrecognize commande: ",sys.argv[1:]
        print _usage
        sys.exit()

    matf_name = "xds2mos.umat"
    mosi_name = "xds2mos.inp"
    
    print "\n   xds2mos version: %s\n" % (__version__)
    print "   Extracting data from:\t\t\t%s" % sys.argv[1]
    if not _read_correct:
        print "   Detector type and image name taken from:\tXDS.INP\n"

    wave = 1/vec3(xdsPar["beam"]).length()
    UB = OM_xds2mos(xdsPar)
    B = BusingLevy(reciprocal(xdsPar["cell"]))
    U = UB * B.inverse() / wave

    xdsPar["origin"] = getBeamOrigin(xdsPar)
    xdsPar["omega"] = getOmega(xdsPar)
    xdsPar["twotheta"] = getTwoTheta(xdsPar)
    
    print "\n   Calculated Omega:     %9.2f degree" % (xdsPar["omega"]*r2d)
    print "   Calculated 2theta:    %9.2f degree\n" % (xdsPar["twotheta"]*r2d)
    
    if _write_out_angles:
        phi123 = ThreeAxisRotation(U.toList(1), inversAxesOrder=1).getAngles()
        phi123 = phi123[0]*r2d
        U = mat3([ex, ey, ez])
    else:    
        phi123 = 0, 0, 0
        
    write_mosflm_matrices(matf_name, UB, U, phi123, xdsPar["cell"])
    mosPar = PARS_xds2best(xdsPar,_read_idxref,_read_correct)
    
    mosi = open(mosi_name,"w")
    if _read_correct:
        mosi.write(mosflmInpTemplate2 % mosPar)
    mosi.write(mosflmInpTemplate1 % mosPar)
    mosi.close()
    
    Ud = UB.decompose()[0]
    Bd = Ud.inverse() * UB / wave
    
    print "   Decomposition Epsilon Ud-U = %8.1e" % diffMAT(Ud, U)
    print "   Decomposition Epsilon Bd-B = %8.1e" % diffMAT(Bd, B)
    print 
    print str_mat(B , '   B decomposition:\n\n', "%14.7f")
    
    print str_mat(UB, "   Mosflm UB:\n\n","%14.7f")
    print "   New Mosflm matix file:  %s" % matf_name
    print "   New Mosflm input file:  %s\n" % mosi_name
    print "   Use -m option to start 'ipmosflm < %s'\n" % mosi_name
    
    if _start_mosflm:
        os.system("ipmosflm <  %s" % mosi_name)
