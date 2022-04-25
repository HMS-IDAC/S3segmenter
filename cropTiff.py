import matplotlib.pyplot as plt
import tifffile
import os
import numpy as np
from scipy.ndimage import *
import scipy.ndimage as ndi
from skimage.morphology import *
from skimage.morphology import extrema
from skimage import morphology
from skimage.measure import regionprops
from skimage.transform import resize
from skimage.filters import gaussian, threshold_otsu, threshold_local
from skimage.feature import peak_local_max
from skimage.color import label2rgb
import skimage.io
from skimage.segmentation import clear_border, watershed, find_boundaries
from scipy.ndimage.filters import uniform_filter
from os.path import *
from os import listdir, makedirs, remove
import pickle
import shutil
import fnmatch
import sys
import argparse
import re
import copy
import datetime
from save_tifffile_pyramid import save_pyramid


if __name__ == '__main__':
    channel = [3]
    channel[:] = [number - 1 for number in channel]
    print(channel)
    
    imagePath = 'V:/cycif-production/114-PCA-ITM/registration/Z238_c10_stitch/Z238_4.ome.tif'
    outputPath = 'V:/cycif-production/114-PCA-ITM/registration/Z238_c10_stitch/cropped'
    fileName = 'Z238_4'
#    subset = np.zeros((20940,42000),dtype=np.uint16)

    metadata = {
        'Pixels': {
            'PhysicalSizeX': '0.65',
            'PhysicalSizeXUnit': '0.65',
            'PhysicalSizeY': 'µm',
            'PhysicalSizeYUnit': 'µm',
        },
        
    }

    append_kwargs = {'bigtiff': True,'metadata': metadata,'append': True}
    save_kwargs = {'bigtiff': True,'metadata': metadata,'append': False}


    
    for iChan in range(40):
        I = tifffile.imread(imagePath,key = iChan)
        subset=I[5527:31017,2400:41450]
        print(iChan+1)
        if iChan == 0:
            skimage.io.imsave(
                    outputPath + os.path.sep + fileName + '.tif',
                    subset, **save_kwargs)
        else:
            skimage.io.imsave(
                    outputPath + os.path.sep + fileName + '.tif',
                    subset, **append_kwargs)




#
#    save_pyramid(
#                subset,
#                outputPath + os.path.sep + fileName + '.ome.tif',
#                pixel_sizes=(0.65,0.65)
#                
#            )         
    
#    subset = np.zeros((12,26400,29630),dtype=np.uint16)
#
#        subset[iChan,:,:]=I[4600:32000,2200:37200]