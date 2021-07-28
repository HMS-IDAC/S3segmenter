import matplotlib.pyplot as plt
import tifffile
import os
import numpy as np
from skimage import io as skio
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
from skimage.io import imsave,imread
from skimage.segmentation import clear_border, watershed, find_boundaries
from scipy.ndimage.filters import uniform_filter
from os.path import *
from os import listdir, makedirs, remove
import pickle
import shutil
import fnmatch
import cv2
import sys
import argparse
import re
import copy
import datetime
from joblib import Parallel, delayed
from rowit import WindowView, crop_with_padding_mask
from save_tifffile_pyramid import save_pyramid
import subprocess
import ome_types


def imshowpair(A,B):
    plt.imshow(A,cmap='Purples')
    plt.imshow(B,cmap='Greens',alpha=0.5)
    plt.show()

    
def imshow(A):
    plt.imshow(A)
    plt.show()
    
def overlayOutline(outline,img):
    img2 = img.copy()
    stacked_img = np.stack((img2,)*3, axis=-1)
    stacked_img[outline > 0] = [65535, 0, 0]
    imshowpair(img2,stacked_img)
    
def normI(I):
    Irs=resize(I,(I.shape[0]//10,I.shape[1]//10) );
    p1 = np.percentile(Irs,10);
    J = I-p1;
    p99 = np.percentile(Irs,99.99);
    J = J/(p99-p1);
    return J

def contour_pm_watershed(
    contour_pm, sigma=2, h=0, tissue_mask=None,
    padding_mask=None, min_area=None, max_area=None
):
    if tissue_mask is None:
        tissue_mask = np.ones_like(contour_pm)
    padded = None
    if padding_mask is not None and np.any(padding_mask == 0):
        contour_pm, padded = crop_with_padding_mask(
            contour_pm, padding_mask, return_mask=True
        )
        tissue_mask = crop_with_padding_mask(
            tissue_mask, padding_mask
        )
    
    maxima = peak_local_max(
        extrema.h_maxima(
            ndi.gaussian_filter(np.invert(contour_pm), sigma=sigma),
            h=h
        ),
        indices=False,
        footprint=np.ones((3, 3))
    )
    maxima = label(maxima).astype(np.int32)
    
    # Passing mask into the watershed function will exclude seeds outside
    # of the mask, which gives fewer and more accurate segments
    maxima = watershed(
        contour_pm, maxima, watershed_line=True, mask=tissue_mask
    ) > 0
    
    if min_area is not None and max_area is not None:
        maxima = label(maxima, connectivity=1).astype(np.int32)
        areas = np.bincount(maxima.ravel())
        size_passed = np.arange(areas.size)[
            np.logical_and(areas > min_area, areas < max_area)
        ]
        maxima *= np.isin(maxima, size_passed)
        np.greater(maxima, 0, out=maxima)

    if padded is None:
        return maxima.astype(np.bool)
    else:
        padded[padded == 1] = maxima.flatten()
        return padded.astype(np.bool)

def S3AreaSegmenter(nucleiPM, images, TMAmask, threshold,fileprefix,outputPath):
    nucleiCenters = nucleiPM[:,:,0]
    TMAmask= (nucleiCenters>np.amax(nucleiCenters)*0.8)*TMAmask
    area = []
    area.append(np.sum(np.sum(TMAmask)))
    for iChan in range(len(images)):
        image_gauss = gaussian(resize(images[iChan,:,:],(int(0.25*images.shape[1]),int(0.25*images.shape[2]))),1)
        if threshold ==-1:
            threshold = threshold_otsu(image_gauss)
        mask=resize(image_gauss>threshold,(images.shape[1],images.shape[2]),order = 0)*TMAmask
        area.append(np.sum(np.sum(mask)))
    np.savetxt(outputPath + os.path.sep + fileprefix + '_area.csv',(np.transpose(np.asarray(area))),fmt='%10.5f')  
    return TMAmask

def getMetadata(path):
    with tifffile.TiffFile(path) as tif:
        metadata=ome_types.from_xml(tif.ome_metadata)
        return metadata.images[0].pixels

def S3NucleiBypass(nucleiPM,nucleiImage,logSigma,TMAmask,nucleiFilter,nucleiRegion):        
        foregroundMask =  nucleiPM
        P = regionprops(foregroundMask, nucleiImage)
        prop_keys = ['mean_intensity', 'label','area']
        def props_of_keys(prop, keys):
            return [prop[k] for k in keys]
        
        mean_ints, labels, areas = np.array(Parallel(n_jobs=6)(delayed(props_of_keys)(prop, prop_keys) 
						for prop in P
						)
		        ).T
        del P
        maxArea = (logSigma[1]**2)*3/4
        minArea = (logSigma[0]**2)*3/4
        passed = np.logical_and(areas > minArea, areas < maxArea)
        
        foregroundMask *= np.isin(foregroundMask, labels[passed])
        np.greater(foregroundMask, 0, out=foregroundMask)
        foregroundMask = label(foregroundMask, connectivity=1).astype(np.int32)
        return foregroundMask

def S3NucleiSegmentationWatershed(nucleiPM,nucleiImage,logSigma,TMAmask,nucleiFilter,nucleiRegion):
    nucleiContours = nucleiPM[:,:,1]
    nucleiCenters = nucleiPM[:,:,0]
    mask = resize(TMAmask,(nucleiImage.shape[0],nucleiImage.shape[1]),order = 0)>0
    if nucleiRegion == 'localThreshold' or nucleiRegion == 'localMax':
        Imax =  peak_local_max(extrema.h_maxima(255-nucleiContours,logSigma[0]),indices=False)
        Imax = label(Imax).astype(np.int32)
        foregroundMask =  watershed(nucleiContours, Imax, watershed_line=True)
        P = regionprops(foregroundMask, np.amax(nucleiCenters) - nucleiCenters - nucleiContours)
        prop_keys = ['mean_intensity', 'label','area']
        def props_of_keys(prop, keys):
            return [prop[k] for k in keys]
        
        mean_ints, labels, areas = np.array(Parallel(n_jobs=6)(delayed(props_of_keys)(prop, prop_keys) 
						for prop in P
						)
		        ).T
        del P
        maxArea = (logSigma[1]**2)*3/4
        minArea = (logSigma[0]**2)*3/4
        passed = np.logical_and.reduce((
            np.logical_and(areas > minArea, areas < maxArea),
            np.less(mean_ints, 50)
            ))
        
        foregroundMask *= np.isin(foregroundMask, labels[passed])
        np.greater(foregroundMask, 0, out=foregroundMask)
        foregroundMask = label(foregroundMask, connectivity=1).astype(np.int32)
    
    else:
        if len(logSigma)==1:
             nucleiDiameter  = [logSigma*0.5, logSigma*1.5]
        else:
             nucleiDiameter = logSigma
        logMask = nucleiCenters > 150
        
        win_view_setting = WindowView(nucleiContours.shape, 2000, 500)
    
        nucleiContours = win_view_setting.window_view_list(nucleiContours)
        padding_mask = win_view_setting.padding_mask()
        mask = win_view_setting.window_view_list(mask)
    
        maxArea = (logSigma[1]**2)*3/4
        minArea = (logSigma[0]**2)*3/4
    
        foregroundMask = np.array(
            Parallel(n_jobs=6)(delayed(contour_pm_watershed)(
                img, sigma=logSigma[1]/30, h=logSigma[1]/30, tissue_mask=tm,
                padding_mask=m, min_area=minArea, max_area=maxArea
            ) for img, tm, m in zip(nucleiContours, mask, padding_mask))
        )
    
        del nucleiContours, mask
    
        foregroundMask = win_view_setting.reconstruct(foregroundMask)
        foregroundMask = label(foregroundMask, connectivity=1).astype(np.int32)
    
        if nucleiFilter == 'IntPM':
            int_img = nucleiCenters
        elif nucleiFilter == 'Int':
            int_img = nucleiImage
    
        print('    ', datetime.datetime.now(), 'regionprops')
        P = regionprops(foregroundMask, int_img)
    
        def props_of_keys(prop, keys):
            return [prop[k] for k in keys]
    
        prop_keys = ['mean_intensity', 'area', 'solidity', 'label']
        mean_ints, areas, solidities, labels = np.array(
            Parallel(n_jobs=6)(delayed(props_of_keys)(prop, prop_keys) 
                for prop in P
            )
        ).T
        del P
    
        MITh = threshold_otsu(mean_ints)
    
        minSolidity = 0.8
    
        passed = np.logical_and.reduce((
            np.greater(mean_ints, MITh),
            np.logical_and(areas > minArea, areas < maxArea),
            np.greater(solidities, minSolidity)
        ))
    
        # set failed mask label to zero
        foregroundMask *= np.isin(foregroundMask, labels[passed])
    
        np.greater(foregroundMask, 0, out=foregroundMask)
        foregroundMask = label(foregroundMask, connectivity=1).astype(np.int32)

    return foregroundMask

def bwmorph(mask,radius):
    mask = np.array(mask,dtype=np.uint8)
    #labels = label(mask)
    background = nucleiMask == 0
    distances, (i, j) = distance_transform_edt(background, return_indices=True)
    cellMask = nucleiMask.copy()
    finalmask = background & (distances <= radius)
    cellMask[finalmask] = nucleiMask[i[finalmask], j[finalmask]]

#    imshowpair(cellMask,mask)
    return cellMask
#    imshow(fg)
#    fg = cv2.dilate(mask,ndimage.generate_binary_structure(2, 2))
#    bg = 1-fg-mask
#    imshowpair(bg,mask)

def S3CytoplasmSegmentation(nucleiMask,cyto,mask,cytoMethod='distanceTransform',radius = 5):
    mask = (nucleiMask + resize(mask,(nucleiMask.shape[0],nucleiMask.shape[1]),order=0))>0
    gdist = distance_transform_edt(1-(nucleiMask>0))
    if cytoMethod == 'distanceTransform':
        mask = np.array(mask,dtype=np.uint32)
        markers= nucleiMask
    elif cytoMethod == 'hybrid':
        cytoBlur = gaussian(cyto,2)
        c1 = uniform_filter(cytoBlur, 3, mode='reflect')
        c2 = uniform_filter(cytoBlur*cytoBlur, 3, mode='reflect')
        grad = np.sqrt(c2 - c1*c1)*np.sqrt(9./8)
        grad[np.isnan(grad)]=0
        gdist= np.sqrt(np.square(grad) + 0.000001*np.amax(grad)/np.amax(gdist)*np.square(gdist))
        bg = binary_erosion(np.invert(mask),morphology.selem.disk(radius, np.uint8))
        markers=nucleiMask.copy()
        markers[bg==1] = np.amax(nucleiMask)+1
        markers = label(markers>0,connectivity=1)
        mask = np.ones(nucleiMask.shape)
        del bg
    elif cytoMethod == 'ring':
#        mask =np.array(bwmorph(nucleiMask,radius)*mask,dtype=np.uint32)>0
        mask =np.array(bwmorph(nucleiMask,radius),dtype=np.uint32)>0
        markers= nucleiMask
    
    cellMask  =clear_border(watershed(gdist,markers,watershed_line=True))
    del gdist, markers, cyto
    cellMask = np.array(cellMask*mask,dtype=np.uint32)
	
    finalCellMask = np.zeros(cellMask.shape,dtype=np.uint32)
    P = regionprops(label(cellMask>0,connectivity=1),nucleiMask>0,cache=False)
    count=0
    for props in P:
         if props.max_intensity>0 :
            count += 1
            yi = props.coords[:, 0]
            xi = props.coords[:, 1]
            finalCellMask[yi, xi] = count
    nucleiMask = np.array(nucleiMask>0,dtype=np.uint32)
    nucleiMask = finalCellMask*nucleiMask
    cytoplasmMask = np.subtract(finalCellMask,nucleiMask)
    return cytoplasmMask,nucleiMask,finalCellMask
    
def exportMasks(mask,image,outputPath,filePrefix,fileName,commit,metadata,saveFig=True,saveMasks = True):
    outputPath =outputPath + os.path.sep + filePrefix
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    metadata_args = dict(
        pixel_sizes=(metadata.physical_size_y,metadata.physical_size_x),
        pixel_size_units=('µm', 'µm'),
        software= 's3segmenter v' + commit
    )
    if saveMasks ==True:
        save_pyramid(
            mask,
            outputPath + os.path.sep + fileName + '.ome.tif',
            channel_names=fileName,
            **metadata_args
        )     
    if saveFig== True:
        mask=np.uint8(mask>0)
        edges = find_boundaries(mask,mode = 'outer')
        stacked_img=np.stack((np.uint16(edges)*np.amax(image),image),axis=0)
        save_pyramid(
            stacked_img,
            outputPath + os.path.sep + fileName + 'Outlines.ome.tif',
            channel_names=[f'{fileName} outlines', 'Segmentation image'],
            **metadata_args
        )
        
def S3punctaDetection(spotChan,mask,sigma,SD):
    Ilog = -gaussian_laplace(np.float32(spotChan),sigma)
    tissueMask = spotChan >0
    fgm=peak_local_max(Ilog, indices=False,footprint=np.ones((3, 3)))*tissueMask
    test=Ilog[fgm==1]
    med = np.median(test)
    mad = np.median(np.absolute(test - med))
    thresh = med + 1.4826*SD*mad
    return (Ilog>thresh)*fgm*(mask>0)
 
     

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument("--imagePath")
    parser.add_argument("--contoursClassProbPath",default ='')
    parser.add_argument("--nucleiClassProbPath",default ='')
    parser.add_argument("--stackProbPath",default ='')
    parser.add_argument("--outputPath")
    parser.add_argument("--dearrayPath")
    parser.add_argument("--maskPath")
    parser.add_argument("--probMapChan",type = int, default = -1)
    parser.add_argument("--mask",choices=['TMA', 'tissue','none'],default = 'tissue')
    parser.add_argument("--crop",choices=['interactiveCrop','autoCrop','noCrop','dearray','plate'], default = 'noCrop')
    parser.add_argument("--cytoMethod",choices=['hybrid','distanceTransform','bwdistanceTransform','ring'],default = 'distanceTransform')
    parser.add_argument("--nucleiFilter",choices=['IntPM','LoG','Int','none'],default = 'IntPM')
    parser.add_argument("--nucleiRegion",choices=['watershedContourDist','watershedContourInt','watershedBWDist','dilation','localThreshold','localMax','bypass','pixellevel'], default = 'watershedContourInt')
    parser.add_argument("--pixelThreshold",type = float, default = -1)
    parser.add_argument("--segmentCytoplasm",choices = ['segmentCytoplasm','ignoreCytoplasm'],default = 'ignoreCytoplasm')
    parser.add_argument("--cytoDilation",type = int, default = 5)
    parser.add_argument("--logSigma",type = int, nargs = '+', default = [3, 60])
    parser.add_argument("--CytoMaskChan",type=int, nargs = '+', default=[2])
    parser.add_argument("--pixelMaskChan",type=int, nargs = '+', default=[2])
    parser.add_argument("--TissueMaskChan",type=int, nargs = '+', default=0)
    parser.add_argument("--detectPuncta",type=int, nargs = '+', default=[0])
    parser.add_argument("--punctaSigma", nargs = '+', type=float, default=[0])
    parser.add_argument("--punctaSD", nargs = '+', type=float, default=[4])
    parser.add_argument("--saveMask",action='store_false')
    parser.add_argument("--saveFig",action='store_false')
    args = parser.parse_args()
    
       
    imagePath = args.imagePath
    outputPath = args.outputPath
    nucleiClassProbPath = args.nucleiClassProbPath
    contoursClassProbPath = args.contoursClassProbPath
    stackProbPath = args.stackProbPath
    maskPath = args.maskPath
    
    metadata = getMetadata(imagePath)
    
    fileName = os.path.basename(imagePath)
    filePrefix = fileName[0:fileName.index('.')]
 
    # convert 1-based indexing to 0-based indexing
    CytoMaskChan = args.CytoMaskChan
    CytoMaskChan[:] = [number - 1 for number in CytoMaskChan]
    pixelMaskChan = args.pixelMaskChan
    pixelMaskChan[:] = [number - 1 for number in pixelMaskChan]
 
    commit = ''#subprocess.check_output(['git', 'describe', '--tags']).decode('ascii').strip()
 
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    
    # get channel used for nuclei segmentation

    if len(contoursClassProbPath)>0:
        legacyMode = 1
        probPrefix = os.path.basename(contoursClassProbPath)
        nucMaskChan = int(probPrefix.split('ContoursPM_')[1].split('.')[0])
    elif len(stackProbPath)>0:
        legacyMode = 0
        probPrefix = os.path.basename(stackProbPath)
        index = re.split('_', stackProbPath)
        index = re.search('.', index[-1]).group(0)
        
        if len(index)==0:
            nucMaskChan = args.probMapChan
        else:
            nucMaskChan  = int(re.sub("\D", "", index))
    else:
        print('NO PROBABILITY MAP PROVIDED')
    if args.probMapChan ==-1:
        if nucMaskChan ==-1:
            sys.exit('INVALID NUCLEI CHANNEL SELECTION. SELECT CHANNEL USING --probMapChan')
        else:
            print('extracting nuclei channel from filename')
    else:
        nucMaskChan = args.probMapChan
    nucMaskChan = nucMaskChan -1 #convert 1-based indexing to 0-based indexing    

    if args.TissueMaskChan==0:
        TissueMaskChan = copy.copy(CytoMaskChan)
        TissueMaskChan.append(nucMaskChan)
    else:
        TissueMaskChan = args.TissueMaskChan[:]
        TissueMaskChan[:] = [number - 1 for number in TissueMaskChan]#convert 1-based indexing to 0-based indexing        
        TissueMaskChan.append(nucMaskChan)

	#crop images if needed
    if args.crop == 'interactiveCrop':
        nucleiCrop = tifffile.imread(imagePath,key = nucMaskChan)
        r=cv2.selectROI(resize(nucleiCrop,(nucleiCrop.shape[0] // 30, nucleiCrop.shape[1] // 30)))
        cv2.destroyWindow('select')
        rect=np.transpose(r)*30
        PMrect= [rect[1], rect[0], rect[3], rect[2]]
        nucleiCrop = nucleiCrop[int(rect[1]):int(rect[1]+rect[3]), int(rect[0]):int(rect[0]+rect[2])]
    elif args.crop == 'noCrop' or args.crop == 'dearray'  or args.crop == 'plate':
        nucleiCrop = tifffile.imread(imagePath,key = nucMaskChan)
        rect = [0, 0, nucleiCrop.shape[0], nucleiCrop.shape[1]]
        PMrect= rect
    elif args.crop == 'autoCrop':
        nucleiCrop = tifffile.imread(imagePath,key = nucMaskChan)
        rect = [np.round(nucleiCrop.shape[0]/3), np.round(nucleiCrop.shape[1]/3),np.round(nucleiCrop.shape[0]/3), np.round(nucleiCrop.shape[1]/3)]
        PMrect= rect
        nucleiCrop = nucleiCrop[int(rect[0]):int(rect[0]+rect[2]), int(rect[1]):int(rect[1]+rect[3])]
        
    if legacyMode ==1:    
        nucleiProbMaps = tifffile.imread(nucleiClassProbPath,key=0)
        nucleiPM = nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
        nucleiProbMaps = tifffile.imread(contoursClassProbPath,key=0)
        PMSize = nucleiProbMaps.shape
        nucleiPM = np.dstack((nucleiPM,nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]))
    else:
        nucleiProbMaps = imread(stackProbPath)
        if len(nucleiProbMaps.shape)==2:
            nucleiPM = nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
        else:
            nucleiPM = nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3]),:]
        PMSize = nucleiProbMaps.shape
      
    # mask the core/tissue
    if args.crop == 'dearray':
        TMAmask = tifffile.imread(maskPath)
    elif args.crop =='plate':
        TMAmask = np.ones(nucleiCrop.shape)
		
    else:
        tissue = np.empty((len(TissueMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)
        count=0
        if args.crop == 'noCrop':
            for iChan in TissueMaskChan:
                tissueCrop =tifffile.imread(imagePath,key=iChan)
                tissue_gauss = gaussian(tissueCrop,1)
                #tissue_gauss[tissue_gauss==0]=np.nan
                tissue[count,:,:] =np.log2(tissue_gauss+1)>threshold_otsu(np.log2(tissue_gauss+1))
                count+=1
        else:
            for iChan in TissueMaskChan:
                tissueCrop = tifffile.imread(imagePath,key=iChan)
                tissue_gauss = gaussian(tissueCrop[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])],1)
                tissue[count,:,:] =  np.log2(tissue_gauss+1)>threshold_otsu(np.log2(tissue_gauss+1))
                count+=1
        TMAmask = np.max(tissue,axis = 0)


        del tissue_gauss, tissue
      
    # nuclei segmentation
    if args.nucleiRegion == 'pixellevel':
        pixelTissue = np.empty((len(pixelMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)
        count=0
        for iChan in pixelMaskChan:
                pixelCrop = tifffile.imread(imagePath,key=iChan)
                pixelTissue[count,:,:] = pixelCrop[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
                count+=1
        nucleiMask = S3AreaSegmenter(nucleiPM, pixelTissue, TMAmask,args.pixelThreshold,filePrefix,outputPath)
    elif args.nucleiRegion == 'bypass':    
        nucleiMask = S3NucleiBypass(nucleiPM,nucleiCrop,args.logSigma,TMAmask,args.nucleiFilter,args.nucleiRegion)
    else:
        nucleiMask = S3NucleiSegmentationWatershed(nucleiPM,nucleiCrop,args.logSigma,TMAmask,args.nucleiFilter,args.nucleiRegion)
    del nucleiPM
    # cytoplasm segmentation
    if args.segmentCytoplasm == 'segmentCytoplasm':
        count =0
        if args.crop == 'noCrop' or args.crop == 'dearray' or args.crop == 'plate':
            cyto=np.empty((len(CytoMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)    
            for iChan in CytoMaskChan:
                cyto[count,:,:] =  tifffile.imread(imagePath, key=iChan)
                count+=1
        elif args.crop =='autoCrop':
            cyto=np.empty((len(CytoMaskChan),int(rect[2]),int(rect[3])),dtype=np.int16)
            for iChan in CytoMaskChan:
                cytoFull= tifffile.imread(imagePath, key=iChan)
                cyto[count,:,:] = cytoFull[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
                count+=1                
        else:
            cyto=np.empty((len(CytoMaskChan),rect[3],rect[2]),dtype=np.int16)
            for iChan in CytoMaskChan:
                cytoFull= tifffile.imread(imagePath, key=iChan)
                cyto[count,:,:] = cytoFull[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
                count+=1
        cyto = np.amax(cyto,axis=0)
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,args.cytoMethod,args.cytoDilation)
        exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,filePrefix,'nuclei',commit,metadata,args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cyto',commit,metadata,args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cell',commit,metadata,args.saveFig,args.saveMask)
  
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,'ring',args.cytoDilation)
        exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,filePrefix,'nucleiRing',commit,metadata,args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cytoRing',commit,metadata,args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cellRing',commit,metadata,args.saveFig,args.saveMask)
        
    elif args.segmentCytoplasm == 'ignoreCytoplasm':
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'nuclei',commit,metadata)
        cellMask = nucleiMask
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'cell',commit,metadata)
        cytoplasmMask = nucleiMask
        
    detectPuncta = args.detectPuncta
    if (np.min(detectPuncta)>0):
        detectPuncta[:] = [number - 1 for number in detectPuncta] #convert 1-based indexing to 0-based indexing    
        if len(detectPuncta) != len(args.punctaSigma):
            args.punctaSigma = args.punctaSigma[0] * np.ones(len(detectPuncta))
 
  
        if len(detectPuncta) != len(args.punctaSD):
            args.punctaSD = args.punctaSD[0] * np.ones(len(detectPuncta))
  
        counter=0
        for iPunctaChan in detectPuncta:
            punctaChan = tifffile.imread(imagePath,key = iPunctaChan)
            punctaChan = punctaChan[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
            spots=S3punctaDetection(punctaChan,cellMask,args.punctaSigma[counter],args.punctaSD[counter])
            cellspotmask = nucleiMask
            P = regionprops(cellspotmask,intensity_image = spots ,cache=False)
            numSpots = []
            for prop in P:
                numSpots.append(np.uint16(np.round(prop.mean_intensity * prop.area)))
            np.savetxt(outputPath + os.path.sep + 'numSpots_chan' + str(iPunctaChan+1) +'.csv',(np.transpose(np.asarray(numSpots))),fmt='%10.5f')    
            edges = 1-(cellMask>0)
            stacked_img=np.stack((np.uint16((spots+edges)>0)*np.amax(punctaChan),punctaChan),axis=0)
            # tifffile.imsave(outputPath + os.path.sep + filePrefix + os.path.sep + 'punctaChan'+str(iPunctaChan+1) + 'Outlines.tif',stacked_img)
            
            
            outputPathPuncta = outputPath + os.path.sep + filePrefix + os.path.sep + 'punctaChan'+str(iPunctaChan+1) + 'Outlines.ome.tif'
            
            metadata_args = dict(
                pixel_sizes=(metadata.physical_size_y, metadata.physical_size_x),
                pixel_size_units=('µm', 'µm'),
                software= 's3segmenter v' + commit
                )
            save_pyramid(
                stacked_img,
                outputPathPuncta,
                channel_names=['puncta outlines', 'image channel'],
                **metadata_args
                )     
            
            counter=counter+1        
        #fix bwdistance watershed