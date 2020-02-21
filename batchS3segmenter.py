# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 21:57:22 2020

@author: cy101
"""

 # read class probability maps
        listing = []
    for iFile in os.listdir(contoursClassProbPath):
        if fnmatch.fnmatch(iFile,filePrefix + '_Contours*'):
            listing.append(iFile)
    PMfileName = ([''.join(listing)])
    
    nucleiPMListing = []
    for iFile in os.listdir(classProbPath):
        if fnmatch.fnmatch(iFile,filePrefix +'_NucleiPM*'):
            nucleiPMListing.append(iFile)
    PMfileName.append(''.join(nucleiPMListing))
    #probMapSuffix.append(PMfileName[PMfileName.index('_NucleiPM_'):])
    
        nucleiProbMaps = tifffile.imread(classProbPath + os.path.sep + PMfileName[0],key=0)
        
        
 # mask
   TMAmask = tifffile.imread(maskPath + os.path.sep + filePrefix + '_mask.tif')