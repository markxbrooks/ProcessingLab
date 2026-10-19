
__version__ = "0.2.2"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "28-04-2008"
__copyright__ = "Copyright (c) 2008 Pierre Legrand"
__license__ = "LGPL"


def getEdgeResolution(PixelX, Width, Distance, Wavelength):
            "Calculate EdgeResolution: Graeme's Method"
            from math import sin, atan
            if Distance > 0.0:
                r = 0.5 * float(PixelX) * int(Width)
                return float(Wavelength)/sin(atan(r/float(Distance)))
            else:
                return 0.

def endian(code):
    if code == 'big_endian': return '>'
    else: return '<'
    
class Interpreter:
    
    HTD = {
    # The adsc Header Translator Dictionary.
    # Potential problems:
    # - There are multiple SIZE1, SIZE2 instances.
    # = The orientation of SIZE1 and SIZE2 is unknown
    #     Not a problem as long as SIZE1 = SIZE2..

    'ExposureTime':(['TIME'],float),
    'BeamX':(['BEAM_CENTER_X'], float),
    'BeamY':(['BEAM_CENTER_Y'], float),
    'Distance':(['DISTANCE'], float),
    'Wavelength':(['WAVELENGTH'], float),
    'PixelX':(['PIXEL_SIZE'], float),
    'PixelY':(['PIXEL_SIZE'], float),
    'Width':(['SIZE1'], int),
    'Height':(['SIZE2'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['OSC_START'], float),
    'PhiEnd':(['OSC_START','OSC_RANGE'], lambda x,y: float(x)+float(y)),
    'PhiWidth':(['OSC_RANGE'], float),
    'EdgeResolution':(['PIXEL_SIZE','SIZE1','DISTANCE','WAVELENGTH'],
        getEdgeResolution),

    # Added keys from Graeme's convention.
    'TwoTheta':(['TWO_THETA'], float),   # _FIXME_ Not really here now...
    'SerialNumber':(['DETECTOR_SN'], str),
    'HeaderSize':(['HEADER_BYTES'], int),
    'EndianType':(['BYTE_ORDER'], endian),
    }

    SpecialRules = {
    
    # No special rules for now
    }
    
    Identifiers = {
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    413:('ESRF','ID14EH2','ADSC Q4'),
    420:('ESRF','ID14EH4','ADSC Q4R'),
    428:('ESRF','ID14EH1','ADSC Q4'),
    444:('ESRF','ID29','ADSC 210'),
    445:('USA?','UNKN','ADSC 210'),
    927:('SOLEIL','PROXIMA1','ADSC 315r')
    }
    
    def getRawHeadDict(self, rawHead):
        _lis = rawHead[2:].split("}")[0].split(";\n")[:-1]
        RawHeadDict = dict([par.split("=") for par in _lis])
        RawHeadDict.update({'MESSAGE':'','TWO_THETA':'0'}) # _FIXME_ Not really here now...
        return RawHeadDict
