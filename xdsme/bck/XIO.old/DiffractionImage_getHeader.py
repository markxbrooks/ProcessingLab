# Extract from the DiffractionImage module written by Graeme Winter
# The DiffractionImage module is distributed in the DNA package (www.dna.ac.uk)
# Used here to ensure compatibility with this DNA module.

import math, os, sys

def GetHeader(imageFileName):
    '''Get back a dictionary with the header information for file filename'''
    
    if not os.environ.has_key('DNAHOME'):
        raise RuntimeError, 'DNAHOME not defined'
    dna = os.environ['DNAHOME']
    sys.path.append(dna + '/scheduler/DiffractionImage/lib')
    #sys.path.append('/home/legrand/prog/DNA/dna-release-rc1/scheduler/DiffractionImage/lib')
    try:
        import DiffractionImage
        d = DiffractionImage.DiffractionImage(imageFileName, 0)
    except:
        raise RuntimeError, 'Error opening image %s' % (imageFileName)

    header = { }
    header['ExposureTime'] = d.getExposureTime()
    header['BeamX'] = d.getBeamX()
    header['BeamY'] = d.getBeamY()
    header['Distance'] = d.getDistance()
    header['Wavelength'] = d.getWavelength()
    header['PixelX'] = d.getPixelX()
    header['PixelY'] = d.getPixelY()
    header['Width'] = d.getWidth()
    header['Height'] = d.getHeight()
    messages = d.getMessage()
    header['Message'] = messages.split(';')
    header['PhiStart'] = d.getPhiStart()
    header['PhiEnd'] = d.getPhiEnd()
    header['PhiWidth'] = d.getPhiEnd() - d.getPhiStart()
    r = 0.5 * d.getPixelX() * d.getWidth()
    A = d.getDistance()

    if A > 0.0:
        theta = math.atan(r / A)
        
        wave = d.getWavelength()
        resolution = wave / math.sin(theta)
        header['EdgeResolution'] = resolution

    else:
        header['EdgeResolution'] = 0.0

    del(d)
    return header

if __name__ == "__main__":
    from pprint import pprint
    pprint(GetHeader(sys.argv[1]))
