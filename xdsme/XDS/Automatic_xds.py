#!/usr/bin/env python
############################################################
# Automatic_xds.py
# Maintained by P. Legrand
# Like the Automatic.py but using XDS instead of mosflm
# 
# This should be run as a "main" routine, viz:
# 
# ./Automatic_xds.py [-directory /data/graeme/images] \
#                [-template foo_1_###.bar] \
#                [-image /data/graeme/images/foo_1_001.bar]
# 
############################################################

import os, sys

if not os.environ.has_key('DNAHOME'):
    raise RuntimeError, 'DNAHOME not defined'

dna = os.environ['DNAHOME']
sys.path.append(dna + '/xsd/python')
#sys.path.append(dna + '/scheduler/DiffractionImage/lib')
#sys.path.append(dna + '/scheduler/Scheduler/XDS')

import XSD
import xupy, xdsSetupDB
import XIO


def do_index_request(inputR):
    _workingDir = "autoindex_xds"
    _idxref_fit_limit = 100
    if inputR.name != 'index_request':
        raise XDSException, 'bad input data'

    images = inputR.getImage()
    images.sort()
    
    fileinfo = inputR.getFileinfo()
    
    DCinfo = xupy.DataCollectInfo(fileinfo.template) 
    DCinfo.directory = fileinfo.directory
    DCinfo.lookup_image_ranges()
    one, last = DCinfo.get_range()

    filename = DCinfo.directory + "/" + DCinfo.format % (one)
    
    # XDS don't read the header: We have to read these informations with Screen.
    # !!!! Getting image header informations !!!!
    img = XIO.Image(filename)
    imageXDSkeys = img.export('xds')
    
    two = DCinfo.getClosestImage(int(one + round(90./img.header['PhiWidth'])))
    
    xpar = xupy.xdsInp2Param(inp_str=xupy.xdsinp_base)
    
    # import Detector + Gogniometer configuration from xdsSetupDB
    xpar.mix(getattr(xdsSetupDB, img.type))
    
    # import experimentatl context information stored in the image header 
    # (wave, dist, pixelsize, osc_range, ...).
    xpar.mix(imageXDSkeys)
    
    if not xpar.X_RAY_WAVELENGTH: xpar.X_RAY_WAVELENGTH = 1.54179
    
    xpar.DATA_RANGE = one, last
    xpar.SPOT_RANGE = (one, one ), (two, two)
    # Look at the Background on the 5 first images if they exist 
    xpar.BACKGROUND_RANGE = one, DCinfo.getClosestImage(min(one+4,last))
            
    if inputR.target:
        sym = (str(inputR.target.getSymmetry())).upper()
    	cell = inputR.target.getCell()
    	if sym and cell:
    	    if sym != 'unknown':
                for spgn in xupy.SPGlib.keys():
                    if sym == xupy.SPGlib[spgn][0]:
                        xpar.SPACE_GROUP_NUMBER = spgn
			break
	if cell:
    	    if cell.a:
		xpar.UNIT_CELL_CONSTANTS = cell.a, cell.b, cell.c, \
		         cell.alpha, cell.beta, cell.gamma

    DCinfo.directory = DCinfo.directory + "/../"
    xpar.NAME_TEMPLATE_OF_DATA_FRAMES = DCinfo.getXDSTemplate()
    if img.compression.upper() == "PCK":
        xpar.NAME_TEMPLATE_OF_DATA_FRAMES += " CCP4" 

    xpar.JOB = "XYCORR", "INIT", "COLSPOT", "IDXREF"
    if not _workingDir in os.listdir("."): os.mkdir(_workingDir)
    os.chdir(_workingDir)
    
    #### RUN XDS #####
    xupy.run_xds(xpar, inp_f=None, out_f="xds_ii.out", save=0)
    #### RUN XDS #####
    
    xupy.saveLastVersion(("XDS.INP",),"_II")
    try:
        xupy.saveLastVersion(("XPARM.XDS","IDXREF.LP"),"_II")
    except:
        raise XDSException, 'bad failed to index'
           
    index_response = XSD.Index_response()
    status = XSD.Status()
    # Starts examination of log files IDXREF.LP
    
    ilp = xupy.opReadCl("IDXREF.LP_II")
    ilp_err = ilp.find("!!! ERROR !!!")
    if ilp_err >= 0:
        status.setCode('error')
        status.setMessage(ilp[ilp_err:].splitlines()[0])
        print "\n  WARNING %s\n\n" % ilp[ilp_err:].splitlines()[0]
    else:
        status.setCode('ok')    
    index_response.setStatus(status)
    
    # Analyse spots uses.
    ilp_ind = ilp.find("INDEXING OF OBSERVED SPOTS IN SPACE GROUP #   1")+54
    indexed = ilp[ilp_ind:ilp_ind+38].split()
    num_indexed, num_total = int(indexed[0]),  int(indexed[3])

    spot_search_response = XSD.Spot_search_response()
    spot_search_response.setFound(num_total)
    spot_search_response.setUsed(num_indexed)
    spot_search_response.setRejected(num_total-num_indexed)
    index_response.setSpot_search_response(spot_search_response)
    
    ilp_xyz = ilp.index("SPOT    POSITION (PIXELS)", ilp_ind)+25
    ilp_xyz = ilp[ilp_xyz:ilp_xyz+67].split()
    stdev_xy, stdev_z = float(ilp_xyz[0]), float(ilp_xyz[-1])
        
    XPARM_II = tuple(map(float, xupy.opReadCl("XPARM.XDS_II").split()))
    
    solution = XSD.Solution()
    shift = XSD.Beam()
    refinement = XSD.Refinement()
    beam_shift = XSD.Beam_shift()
    reflections = XSD.Reflections()
    
    new = XSD.Beam()
    new.setX(XPARM_II[16]*xpar.QY)
    new.setY(XPARM_II[15]*xpar.QX)
    
    shift.setX(abs(xpar.ORGY - XPARM_II[16])*xpar.QY)
    shift.setY(abs(xpar.ORGX - XPARM_II[15])*xpar.QX)

    beam_shift.setNew_beam(new)
    beam_shift.setShift(shift)
    refinement.setBeam_shift(beam_shift)

    reflections.setUsed(num_indexed)
    reflections.setUsed_in_indexing(num_indexed)
    refinement.setReflections(reflections)

    solution.setRefinement(refinement)
    solution.setSpot_deviation(stdev_xy*xpar.QX)
    index_response.setSolution(solution)
    
    initial_orientation_matrix = XSD.Matrix()
    reducedCell =  XPARM_II[27:33]
    reducedCell_volum = xupy.Lattice(reducedCell,"aP").volum()
    fmt_cell = 6*"%7.2f"

    initial = XSD.Initial()
    initial_cell = XSD.Cell()
    attrCell = 'A', 'B', 'C', 'Alpha', 'Beta', 'Gamma'
    for i in range(len(attrCell)):
        method = getattr(initial_cell, 'set%s' % attrCell[i])
        method(reducedCell[i])
    initial.setCell(initial_cell)

    u_matrix = XSD.Matrix()
    initial_orientation_matrix = XSD.Matrix()
    for attr in ['11', '12', '13', '21', '22', '23', '31', '32', '33']:
        method = getattr(initial_orientation_matrix, 'setE%s' % attr)
        method(0.0)
        method = getattr(u_matrix, 'setE%s' % attr)
        method(0.0)
    initial.setOrientation_matrix(initial_orientation_matrix)
    solution.setInitial(initial)
    
    # The XDS A matrix describe DIRECT space vectors coordinate!!!
    a_matrix = XSD.Matrix()
    a_matrix.setE11(XPARM_II[33])
    a_matrix.setE21(XPARM_II[34])
    a_matrix.setE31(XPARM_II[34])
    a_matrix.setE12(XPARM_II[35])
    a_matrix.setE22(XPARM_II[36])
    a_matrix.setE32(XPARM_II[37])
    a_matrix.setE13(XPARM_II[38])
    a_matrix.setE23(XPARM_II[39])
    a_matrix.setE33(XPARM_II[40])
    
    refined_cell = XSD.Cell()
    for i in range(len(attrCell)):
        method = getattr(refined_cell, 'set%s' % attrCell[i])
        method(reducedCell[i])
    
    orientation = XSD.Orientation()
    orientation.setCell(refined_cell)
    orientation.setU_matrix(u_matrix)
    orientation.setA_matrix(a_matrix)
    solution.setOrientation(orientation)
    
    solution.setNumber(0)
    solution.setSpacegroup_number(xpar.SPACE_GROUP_NUMBER)
    solution.setSymmetry(inputR.target.getSymmetry())
    index_response.setSolution(solution)
    index_response.setMosaicity_value(0.)
    
    print "\nSelected reduced cell: ", fmt_cell % reducedCell
    print "Fraction of indexed spots","(%d/%d): %7.1f%%" %  \
                    (num_indexed, num_total,100.*num_indexed/num_total)
    print "Standard deviation of spot    position (in pixels): %8.2f" % stdev_xy,
    print " (in mm): %8.3f" % (stdev_xy*xpar.QX)
    print "Standard deviation of spindel position (degrees):   %8.2f" % stdev_z
    print "  \n\n  Selected list of possible Bravais lattices\n"
    print "\n  Bravais Lattice Symm   fit"+15*" "+"cell"+17*" "+ \
                "ORDER\n\t"+26*" "+"a     b     c  alp  bet  gam\n\n"

    for LAT in xupy.select_lattices(limit=100, idxref="IDXREF.LP_II"):
        rc_au_vol, au_vol = "", ""
        LAT.idealize()
        brav = xupy.BravaisDico[LAT.Bravais_type[0]]
        if LAT.fit <= _idxref_fit_limit:
             au_vol = LAT.multiplicity/(LAT.volum()/reducedCell_volum)
             rc_au_vol = int(round(au_vol))*"*"
             au_vol = "%.0f" % au_vol
        fmt_l = "  %1s %-13s (%2s) %5.1f "+"%s"+"%3s  %s"
        print fmt_l % (LAT.Bravais_type[1],brav,LAT.Bravais_type,LAT.fit,
                         LAT.prt(fmt=3*"%6.0f"+3*"%5.0f"),au_vol,rc_au_vol)
    print "\n  ORDER is defined as the ratio: reduce cell volum/asym. unit volum.\n"
    
    
    return index_response

if __name__ == '__main__':

    template = None
    directory = None
    args = sys.argv
    a, b, c, alpha, beta, gamma = 0, 0, 0, 0, 0, 0
    spg = "unknown"

    for i in range(len(args)):
        if args[i] == '-directory':
            directory = args[i + 1]
        if args[i] == '-template':
            template == args[i + 1]
        if args[i] == '-image':
            image = args[i + 1]
        if args[i] == '-cell':
            cell = args[i + 1]
            a, b, c, alpha, beta, gamma = map(float, cell.split())
        if args[i] == '-spg':
            spg = args[i + 1]


    if template and directory:
        imgI = xupy.DataCollectInfo(template)
        imgI.directory = directory
    elif image:
        imgI = xupy.DataCollectInfo(image)
        template = imgI.getMosflmTemplate()
        directory = imgI.directory
    else:
        print "Not enougth information to start"
        sys.exit()

    images = imgI.lookup_image_numbers()
    if len(images) == 0:
        raise RuntimeError, 'No images match the template %s in directory %s' % \
              (template, directory)

    imgI.lookup_image_ranges()
    one, last = imgI.get_range()

    filename = imgI.directory + "/" + imgI.format % (one)
    header = XIO.Image(filename).header

    two = imgI.getClosestImage(int(one + round(90./header['PhiWidth'])))

    index_xml = "<?xml version='1.0'?><!DOCTYPE index_request><index_request>\
    <fileinfo><template>%s</template><directory>%s</directory></fileinfo>\
    <target><symmetry>%s</symmetry><cell><a>%f</a><b>%f</b><c>%f</c>\
    <alpha>%f</alpha><beta>%f</beta><gamma>%f</gamma></cell></target>\
    <image>%d</image><image>%d</image></index_request>" % \
    (template, directory, spg, a, b, c, alpha, beta, gamma, one, two)

    index_request = XSD.Index_request()
    index_request.unmarshal(index_xml)
    index_response = do_index_request(index_request)
    print index_response.marshal()
