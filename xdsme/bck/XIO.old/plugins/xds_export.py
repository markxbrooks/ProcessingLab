
_to_XDS_detectors = {
"mar":       "MAR345",
"mar555":    "MAR345",
"marccd":    "CCDCHESS",
"adsc":      "ADSC",
"raxis":     "RAXIS",
}


#     Header Translator Dictionary.
#     Translate image header entries in a new dictionay
#     newdic['X_RAY_WAVELENGTH'] = float(head['Wavelength'])
#
HTD = {
'X_RAY_WAVELENGTH':(['Wavelength'], float),
'DETECTOR_DISTANCE':(['Distance'], float),
'STARTING_ANGLE':(['PhiStart'], float),
'OSCILLATION_RANGE':(['PhiWidth'], float),
'NX':(['Width'], int),
'NY':(['Height'], int),
'QX':(['PixelX'], float),
'QY':(['PixelY'], float),
'ORGX':(['BeamY','PixelX'], lambda x,y: float(x)/float(y)), # Yes, it is in revers order.
'ORGY':(['BeamX','PixelY'], lambda x,y: float(x)/float(y)), # DoNotAsk!
}
        
#     Collect Translator Dictionary.
#     Translate collect object attributes to a new dictionay
#     newdic['SPOT_RANGE'] = list(collect.imageRanges)
#
CTD = {
'NAME_TEMPLATE_OF_DATA_FRAMES':(['xdsTemplate'], str),
'DATA_RANGE':(['imageNumbers'], lambda x: [x[0],x[-1]]),
'SPOT_RANGE':(['imageRanges'], list),
'BACKGROUND_RANGE':(['imageRanges'], lambda x: [x[0][0],min(x[0][0]+5,x[0][1])]),
'DETECTOR':(['imageType'], lambda x: _to_XDS_detectors[x]),
}
