#!/usr/bin/env python

__version__ = "0.2.3"
__author__ = "Pierre Legrand (legrand@emble-grenoble.fr)"
__date__ = "01-10-2003"
__copyright__ = "Copyright (c) 2003 Pierre Legrand"
__license__ = "GPL"

"""
TODO:	- Verify the coherence between ProDC infos and image headers.
        - Creer des repertoir + fichier de commande pour XDS, mosflm...
"""

_esrf_bm14_beam_x, _esrf_bm14_beam_y = 114.260, 111.550
_esrf_id14eh3_beam_x, _esrf_id14eh3_beam_y = 79.2, 82.3


import sys, os, struct

_progname = os.path.split(sys.argv[0])[1]
usage   = """
>>>   Usage : %s [-a] ProDC_file.log [first_image_location]\n

 -a:  If you want to keep the anomalous information (anomalous = yes). 
 -n:  If you don't want to keep the anomalous information (anomalous = no). 

""" % _progname

# default values
_angle_det = 0.000
_resolution = "50 0.0"
_anomalous = "no"
_suffix = "img"
_beamline = "ESRF_ID29"
_beam_x = 106.0
_beam_y = 103.2
_img_dir_link = "img_dir"
_spg = 0

adp_template = """beamline = %(beamline)s
frame_first = %(frame_first)d
frame_last = %(frame_last)d
phi_init = %(phi_init).3f
delta_phi = %(delta_phi).3f
prefix = %(prefix)s_%(run_number)d_
suffix = %(suffix)s
beam_x = %(beam_x).3f
beam_y = %(beam_y).3f
anomalous = %(anomalous)s
wavelength = %(wavelength).6f
distance = %(distance).2f
angle_det = %(angle_det).2f
img_dir = %(img_dir)s
resolution_range = %(resolution)s
spg = %(spg)d
"""

info_show = """beamline = %(beamline)s
images = %(img_dir)s/%(prefix)s_%(run_number)d_\?\?\?%(suffix)s
frame_first = %(frame_first)d
frame_last = %(frame_last)d
phi_init = %(phi_init).3f
delta_phi = %(delta_phi).3f
anomalous = %(anomalous)s
wavelength = %(wavelength).6f
distance = %(distance).2f
"""

xds_template = """
 JOB= ALL
 
 DATA_RANGE= %(frame_first)d %(frame_last)d
 SPOT_RANGE= %(spot_range_1f)d %(spot_range_1l)d
 SPOT_RANGE= %(spot_range_2f)d %(spot_range_2l)d
 BACKGROUND_RANGE= %(frame_first)d %(bkg_range_l)d
 
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 8
 STRONG_PIXEL= 8.0
 
 OSCILLATION_RANGE= %(delta_phi).3f
 STARTING_ANGLE= %(phi_init).3f
 STARTING_FRAME= %(frame_first)d
 X-RAY_WAVELENGTH=  %(wavelength).5f
 NAME_TEMPLATE_OF_DATA_FRAMES= %(frame_template)s
 
 DETECTOR_DISTANCE= %(distance).2f
 DETECTOR= %(detector_xds)s   MINIMUM_VALID_PIXEL_VALUE= 1   OVERLOAD= 65000
 ORGX= %(orgx).2f    ORGY= %(orgy).2f
 RMAX= 1.1
 NX=  %(ny)d   NY=  %(nx)d   QX= %(pixel_size_y).4f   QY= %(pixel_size_x).4f
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000
 
 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0
 ROTATION_AXIS= 1.0 0.0 0.0
 INCIDENT_BEAM_DIRECTION= 0.0 0.0 1.0
 FRACTION_OF_POLARIZATION= 0.98
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0
 AIR= 0.001
 
 SPACE_GROUP_NUMBER= 0
 UNIT_CELL_CONSTANTS= 0 0 0 0 0 0
 INCLUDE_RESOLUTION_RANGE= 50.0 0.0
 
 REFINE(INTEGRATE)= DISTANCE BEAM ORIENTATION CELL   
 DELPHI= %(delphi).3f
 MAXIMUM_NUMBER_OF_PROCESSORS= 16
 """


class collect:
    def __init__(self):
        self.angle_det = _angle_det
        self.resolution = _resolution
        self.anomalous = _anomalous
        self.beamline = _beamline
        self.beam_x = _beam_x
        self.beam_y = _beam_y
        self.spg = _spg
    def update(self, New):
        self.__dict__.update(vars(New))

class dummy: pass

detectorSN2beamline = {
    428:"ESRF_ID14EH1",
    413:"ESRF_ID14EH2",
    420:"ESRF_ID14EH4",
    444:"ESRF_ID29",
      3:"ESRF_BM14",
      8:"ESRF_ID13"}

mad_beamlines = ["ESRF_BM14","ESRF_BM30A","ESRF_ID14EH4","ESRF_ID29"]

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def opWriteCl(filename, _str):
    f = open(filename,"w")
    r = f.write(_str)
    f.close()

def get_filename(c):
    return "%s%s_%d_%03d.%s" % (c.img_dir,c.prefix,
                                c.run_number, c.frame_first,c.suffix)

prodc_params = {"Run Number::":"run_number",
                "Filename prefix:":"prefix",
                "Starting Frame Number:":"frame_first",
                "Number of Frames:":"nframe",
                "Number of Passes:":"passes",
                "Rotation Oscillation (degree):":"delta_phi",
                "Oscillation overlap (degree):":"overlape",
                "Starting Angle (degree):":"phi_init",
                "Detector Distance (cm):":"distance",
                "Detector Distance (m):":"distance",
                "Wavelength:":"wavelength",
                "Data Directory:":"img_dir",
                "Filename suffix:":"suffix",
                "Run Comments::":"comments"}

def read_ProDC_log(file_name_in):
    logf = opReadCl(file_name_in)
    _logs = [_l for _l in logf.split("Collection finished") if len(_l) > 100]
    return _logs 
    
def guess_imageType(image_name):
    "Return a the image detector type and compress type"
    image = open(image_name,"r")
    head = image.read(9000)
    image.close()
    #print head
    imageType = []
    imageCompress = "NO"
    if  head[:15] == "{\nHEADER_BYTES=" and head.count(";\nDETECTOR_SN=") and \
        head.count(";\nPIXEL_SIZE="): imageType.append("ADSC")
    if  head.count("\nSCANNER") and head.count("mar research") and \
        head.count("\nPROGRAM"): imageType.append("MarIP")
    if  head.count("Readout from 4 Corner") and \
        head[1028:1031] == "MMX" and \
        head.count("CCD Temperature"): imageType.append("MarCCD")
    if head.count("CCP4 packed image"): imageCompress = "PCK"
    #print imageType
    if len(imageType) == 1: return imageType[0], imageCompress
    elif len(imageType) == 0: return "Unknown", imageCompress
    elif len(imageType) > 1:
        print "ERROR: Can't choose the detector type between:", imageType
        sys.exit()

def read_image_header(image_name, image_type):
    if image_type == "ADSC":
        "Return a dictonary"
        image = open(image_name,"r")
        header = image.read(512).split("}")[0][2:-1]
        head_par = {}
        for line in header.splitlines():
            key, val = line[:-1].split("=")
            try: val = float(val)
            except: pass
            head_par[key] = val
        return head_par
    if image_type == "Raxis":
        "Works well with RaxisIIc images"
        
        header_fmt  = '>10s10s20s12s'+6*'f'+ \
              '12sf80s84s12s20s4sf20sf20s4sfff12s80slf56s4s4s' + \
              12*'f'+'200sllffllllfff10s10slll180s' + \
              '20s'+16*'f'+5*'l'+9*'f'+6*'l'+'92s768s'
              
        keys = ['Device name', 'Versions',
        'Crystal name', 'Crystal system', 'Cell Parameter a',
        'Cell Parameter b', 'Cell Parameter c', 'Cell Parameter alpha ',
        'Cell Parameter beta', 'Cell Parameter gamma', 'Space group',
        'Mosaic spread', 'Memo', 'Date',
        'Operator name', 'X-Ray target', 'Wave length',
        'Monochromator', 'Monochromator 2theta', 'Colliomator',
        'Filter', 'Camera length (mm)', 'Tube Voltage  (KV)  ',
        'Tube Current  (mA)', 'X-Ray Focus', 'X-Ray Optics',
        'Detector Shape Code', 'Weissenberg Coupling', 'Mounted axis',
        'X-ray beam axis', 'Phi 0', 'Phi start',
        'Phi end', 'Oscillation number  ', 'Exposure time (min)',
        'Direct beam X(pixel)', 'Z(pixel)', 'Omega',
        'Chi', 'Theta', 'Mu',
        'Number of pixel (X)', 'Number of pixel (Z)', 'Pixel size X (mm)',
        'Pixel size Z (mm)', 'Record length (byte)', 'Number of records',
        'Start line (Z)', 'IP number', 'Ratio',
        'Fading function No.1', 'Fading function No.2', 'Host Type',
        'IP Type', 'Image Direction X', 'Image Direction Z',
        'Image Direction XZ']
        image = open(image_name,"r")
        read_size  = struct.calcsize(header_fmt)                # calcul le nombre d'octet a lire
        read_raw   = image.read(read_size)      # raw read in one pass
        read_unp   = list(struct.unpack(header_fmt, read_raw))  # uncoding from header_fmt
        head_par = {}
        del read_unp[13] ; del read_unp[28] ; del read_unp[41:43]
        for l in range(len(keys)): head_par[keys[l]] = read_unp[l]
        image.close()
        return head_par 
        
    if image_type == "MarCCD":
        "Not finished"
        
        header_fmt = '<1024s28sII44sIIIIIII8sIIIIII500sIII16si12si28si12sii32sii128sI2160s'
        keys = ['tiff_stuff', 'gunk1', 'header_byte_order','data_byte_order',
        'gunk1a','pixel_count_x','pixel_count_y','bytes_per_pixel',
        'bytes_per_line','bits_per_pixel','data_type','saturated_pixel',
        'gunk2','origin_location','detector_orientation','view_direction',
        'overflow_location','overflow_8_bits','overflow_16_bits','gunk3',
        'distance','beam_x','beam_y','gunk4','two_theta','gunk5','phi_start',
        'gunk6','phi_end','gunk7','axis_code','phi_width','gunk8',
        'pixel_size_x','pixel_size_y','gunk10','lambda','gunk11']
        
        image = open(image_name,"r")
        read_size  = struct.calcsize(header_fmt)
        read_raw   = image.read(read_size)      # raw read in one pass
        read_unp   = list(struct.unpack(header_fmt, read_raw))  # uncoding from header_fmt
        if read_unp[2] != 1234:
            header_fmt = ">" + header_fmt[1:]
            read_unp   = list(struct.unpack(header_fmt, read_raw))
        head_par = {}
        for l in range(len(keys)): head_par[keys[l]] = read_unp[l]
        head_par["beam_x"] = 0.001 * head_par["beam_x"]
        head_par["beam_y"] = 0.001 * head_par["beam_y"]
        head_par["two_theta"] = 0.001 * head_par["two_theta"]
        head_par["lambda"] = 0.00001 * head_par["lambda"]
        head_par["distance"] = 0.001 * head_par["distance"]
        head_par["phi_start"] = 0.001 * head_par["phi_start"]
        head_par["phi_width"] = 0.001 * head_par["phi_start"]
        head_par["phi_end"] = 0.001 * head_par["phi_end"]
        head_par["pixel_size_x"] = 0.000001 * head_par["pixel_size_x"]
        head_par["pixel_size_y"] = 0.000001 * head_par["pixel_size_y"]
        
        #for k in keys:
        #    if k[:3] != "gun": print k,head_par[k]
        image.close()
        return head_par

def decompose_file_location(file_location, ndigit=3):
    FL = dummy()
    FL.img_dir, _imageName = os.path.split(file_location)
    _prefix = _imageName.split(".")[0][:(-1*ndigit-1)].split("_")
    FL.suffix = _imageName.split(".")[-1]
    FL.run_number = int(_prefix[-1])
    FL.prefix = "_".join(_prefix[:-1])
    #print vars(FL)
    return FL

if __name__=='__main__':

    file_location = ""
    arg_anomalous = ""
    if sys.argv.count("-a"):
        sys.argv.remove("-a")
        arg_anomalous = "yes"

    if sys.argv.count("-n"):
        sys.argv.remove("-n")
        arg_anomalous = "no"

    if len(sys.argv) == 1:
        print usage
        sys.exit()
    
    if len(sys.argv) == 3:
        file_location = sys.argv[2]
        imageGiven = 1
    else: imageGiven = 0

    try:
        logs = read_ProDC_log(sys.argv[1])
    except IOError:
        print "\n  ERROR can't open ProDC logfile: %s\n" %  sys.argv[1]
        sys.exit()
    collects = []
    for log in logs:

        c = collect()
        if not imageGiven: file_location == ""
        
        for p in prodc_params.keys():
            ip = log.find(p)+len(p)
            sp = log[ip:].splitlines()[0].strip()
            if log.count(p):
                try:
                    spf = float(sp)
                    sp = spf
                except: pass
                #print "%35s  %s" % (p, sp)
                setattr(c,prodc_params[p],sp)
                 
        #print vars(c)
        if len(vars(c)) >= 15:
            if c.img_dir[-1] != "/": c.img_dir = c.img_dir[-1] + "/"
            print "\n\t>>>> NEW COLLECT <<<<\n"
            c.frame_last = c.frame_first + c.nframe
            collects.append(c)
            try: suffix_init = c.suffix
            except: suffix_init = c.suffix = "img"

            if imageGiven:
                #if not file_location.count("/"):
                #    file_location = "./"+file_location
                newc = decompose_file_location(file_location)
                c.update(newc)
            else: file_location = get_filename(c)
            
            possible_par = ({}, {"img_dir":"./"}, {"suffix":"pck"},
                        {"img_dir":os.path.dirname(sys.argv[1])+"/"},
                        {"img_dir":"./", "suffix":"pck"})
            
            good_par, checked = None, []
            for par in possible_par:
                _c = collect()
                _c.update(c)
                _c.__dict__.update(par)
                file_location = get_filename(_c)
                checked.append(file_location)
                if os.path.isfile(file_location):
                      print "  Reading header informations from file %s\n" % \
                            file_location
                      good_par = par.copy()
            
            if not good_par:
                print "  Image file not found in these locations:\n\t%s\n" % \
                                             ("\n\t".join(checked))
                print "  RATS! Try to add the file location " \
                                     +"on the commande line.\n"
                sys.exit()
            else:
                c.__dict__.update(good_par)
                file_location = get_filename(c)
            
            imageType, imageCompress = guess_imageType(file_location)
                                            
            try:
                head = read_image_header(file_location,imageType)
                print "  Guess of the detector type:\t%s\n" % imageType
            except:
                print "  ERROR: Can't read header of file %s." % file_location
            if imageType == "Unknown":
                print "  WARNING: Using approximate beam center coordinates.\n"
                c.beam_x = 106.0
                c.beam_y = 103.2
                c.beamline = "Unknown"
                c.pixel_size_x = c.pixel_size_y = 0.1
            if imageType == "ADSC":                                      
                c.beamline = detectorSN2beamline[head["DETECTOR_SN"]]    
                c.beam_x = head["BEAM_CENTER_X"]
                c.beam_y = head["BEAM_CENTER_Y"]
                c.pixel_size_x = c.pixel_size_y = head["PIXEL_SIZE"]
                c.nx, c.ny = head["SIZE1"], head["SIZE2"]
            if imageType == "MarCCD":
                print "  WARNING: Using approximate beam center coordinates."
                c.pixel_size_x = head["pixel_size_x"]
                c.pixel_size_y = head["pixel_size_y"]
                c.nx = head["pixel_count_x"]
                c.ny = head["pixel_count_y"]
                if 0.064 < head["pixel_size_x"] < 0.065 :                      
                    c.beamline = "ESRF_BM14"                             
                    c.beam_x = _esrf_bm14_beam_x
                    c.beam_y = _esrf_bm14_beam_y
                if 0.7319 < head["pixel_size_x"] < 0.7320:
                    c.beamline = "ESRF_ID14EH3"
                    c.beam_x = _esrf_id14eh3_beam_x
                    c.beam_y = _esrf_id14eh3_beam_y
            
            if c.img_dir[0] != "/": c.img_dir = "../" + c.img_dir
            if c.img_dir == ".././": c.img_dir = "../"
            
            if c.beamline in mad_beamlines:
                 c.anomalous = "yes"
                 print "\n  INFO: MAD Beamline recognized. Turning anomalous on."
            if arg_anomalous == "no": c.anomalous = "no"
            elif arg_anomalous == "yes": c.anomalous = "yes"
            
            char =  97
            dir_prefixs = "adp_", "xds_" #, "mosflm_"
            dir_suffix = "%s_%d%s" % (c.prefix,c.run_number,chr(char))
            while len([p for p in dir_prefixs if os.path.exists(p+dir_suffix)]):
                char += 1
                dir_suffix = "%s_%d%s" % (c.prefix,c.run_number,chr(char))
                
            print "  New directory created: %s\n" % dir_suffix
            for p in dir_prefixs: os.makedirs(p+dir_suffix)
            
            
            # Writing adp input file in adp_"suffix"
            new_adp_dir = "adp_"+dir_suffix
            #os.symlink(c.img_dir,new_adp_dir+"/"+_img_dir_link)
            #c.img_dir_link = _img_dir_link
            print adp_template % vars(c)
            opWriteCl(new_adp_dir+"/auto.par",adp_template % vars(c))
            
            c.spot_range_1f, c.spot_range_1l = c.frame_first, c.frame_first+5
            c.spot_range_2f, c.spot_range_2l = c.frame_last-5, c.frame_last
            c.bkg_range_l = c.frame_first+4
            c.delphi = c.delta_phi*4
            c.detector_xds = "ADSC"
            c.orgx, c.orgy = c.beam_y/c.pixel_size_y , c.beam_x/c.pixel_size_x
            
            # Writing xds input file in xds_"suffix"
            new_xds_dir = "xds_"+dir_suffix
            #os.symlink(c.img_dir,new_xds_dir+"/"+_img_dir_link)
            #c.img_dir = _img_dir_link
            c.frame_template = "%s/%s_%d_%3s.%s" % (c.img_dir,c.prefix,
                                        c.run_number, "???", c.suffix)
            opWriteCl(new_xds_dir+"/XDS.INP", xds_template % vars(c))
            
