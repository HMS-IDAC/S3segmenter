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
from skimage.measure import label
from skimage.transform import resize
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max
from skimage.color import label2rgb
from skimage.io import imsave
from os.path import *
from os import listdir, makedirs, remove
from sklearn.cluster import KMeans
import pickle
import shutil
import fnmatch
import cv2
import sys
import argparse
import re


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

def S3NucleiSegmentationWatershed(nucleiPM,nucleiImage,logSigma,TMAmask,nucleiFitler,nucleiRegion):
    nucleiContours = nucleiPM[:,:,1]
    nucleiCenters = nucleiPM[:,:,0]
    mask = resize(TMAmask,(nucleiImage.shape[0],nucleiImage.shape[1]),order = 0)>0
    
    if len(logSigma)==1:
         nucleiDiameter  = [logSigma*0.5, logSigma*1.5]
    else:
         nucleiDiameter = logSigma
    logMask = nucleiCenters > 150
   # dist_trans_img = ndi.distance_transform_edt(logMask)
    
    fgm=peak_local_max(extrema.h_maxima(gaussian_filter(255-nucleiContours,logSigma[1]/30),logSigma[1]/30), indices=False,
                                      footprint=np.ones((3, 3)))
    #imshowpair(fgm.astype(np.uint8)*255,nucleiContours)
    _, fgm= cv2.connectedComponents(fgm.astype(np.uint8))
    
    foregroundMask= morphology.watershed(nucleiContours,fgm,watershed_line=True)
   
    allNuclei = ((foregroundMask)*mask)
    P = regionprops(allNuclei,nucleiCenters)
    mean_int = np.array([prop.mean_intensity for prop in P]) 
    #kmeans = KMeans(n_clusters=2).fit(mean_int.reshape(-1,1))
    MITh = threshold_otsu(mean_int.reshape(-1,1))
    X = np.zeros(foregroundMask.shape,dtype=np.uint16)
    count = 0
    maxArea = (logSigma[1]**2)*3/4
    minArea = (logSigma[0]**2)*3/4

    for props in P:
        intensity = props.mean_intensity
        area = props.area
        if intensity >MITh and minArea < area < maxArea:
            count += 1
            yi = props.coords[:, 0]
            xi = props.coords[:, 1]
            X[yi, xi] = count
    
#    idx = [props.label for props in P if props.mean_intensity > MITh and minArea < props.area < maxArea]
#    final_mask =label(np.isin(allNuclei, idx))
    return X
    
#    img2 = nucleiImage.copy()
#    stacked_img = np.stack((img2,)*3, axis=-1)
#    stacked_img[X > 0] = [65535, 0, 0]
#    imshowpair(img2,stacked_img)

def bwmorph(mask,radius):
    mask = np.array(mask,dtype=np.uint8)
    labels = label(mask)
    background = labels == 0
    distances, (i, j) = distance_transform_edt(background, return_indices=True)
    cellMask = labels.copy()
    finalmask = background & (distances <= radius)
    cellMask[finalmask] = labels[i[finalmask], j[finalmask]]

#    imshowpair(cellMask,mask)
    return cellMask
#    imshow(fg)
#    fg = cv2.dilate(mask,ndimage.generate_binary_structure(2, 2))
#    bg = 1-fg-mask
#    imshowpair(bg,mask)

def S3CytoplasmSegmentation(nucleiMask,cyto,mask,cytoMethod='distanceTransform',radius = 5):
    mask = (nucleiMask + resize(mask,(nucleiMask.shape[0],nucleiMask.shape[1]),order=0))>0

#    if cytoMethod == 'bwdistanceTransform':
#    gdist = distance_transform_edt(1-(nucleiMask>0))
#    nucleiMask = np.array(nucleiMask,dtype=np.uint8)
#    markers = label(nucleiMask)
#    cellMask  = watershed(gdist,markers)
#    cellMask = cellMask*mask
    
#    elif cytoMethod == 'distanceTransform':
    gdist = distance_transform_edt(1-(nucleiMask>0))
    nucleiMask = np.array(nucleiMask,dtype=np.uint16)
    mask = np.array(mask,dtype=np.uint16)
    markers = label(nucleiMask)
    cellMask  = watershed(cyto,markers,watershed_line=True)
    cellMask = cellMask*mask
# overlayOutline(np.array(cellMask==0,dtype=np.uint8),nucleiMask)
#    elif cytoMethod == 'ring':
#    cellMask =bwmorph(nucleiMask>0,radius)
    cellMask =  np.array(cellMask,dtype=np.uint16)
    cytoplasmMask = cv2.subtract(cellMask,nucleiMask)
#    imshow(cytoplasmMask)
    return cytoplasmMask,nucleiMask,cellMask
    
def exportMasks(mask,image,outputPath,filePrefix,fileName,saveFig=True,saveMasks = True):
    outputPath =outputPath + os.path.sep + filePrefix
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    if saveMasks ==True:
        kwargs={}
        kwargs['bigtiff'] = True
        kwargs['photometric'] = 'minisblack'
        resolution = np.round(1)
        kwargs['resolution'] = (resolution, resolution, 'cm')
        kwargs['metadata'] = None
        kwargs['description'] = '!!xml!!'
        tifffile.imsave(outputPath + os.path.sep + fileName + 'Mask.tif',mask)
        
    if saveFig== True:
        mask=np.uint8(mask>0)
        edges=cv2.Canny(mask,0,1)
        stacked_img=np.stack((np.uint16(edges)*65535,image),axis=0)
        tifffile.imsave(outputPath + os.path.sep + fileName + 'Outlines.tif',stacked_img)
        
    
        
    # assign nan to tissue mask


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument("--imagePath")
    parser.add_argument("--contoursClassProbPath")
    parser.add_argument("--nucleiClassProbPath")
    parser.add_argument("--outputPath")
    parser.add_argument("--dearrayPath")
    parser.add_argument("--maskPath")
    parser.add_argument("--probMapChan",type = int, default = -1)
    parser.add_argument("--mask",choices=['TMA', 'tissue','none'],default = 'tissue')
    parser.add_argument("--crop",choices=['interactiveCrop','autoCrop','noCrop','dearray'], default = 'noCrop')
    parser.add_argument("--cytoMethod",choices=['distanceTransform','bwdistanceTransform','ring'],default = 'distanceTransform')
    parser.add_argument("--nucleiFilter",choices=['IntPM','LoG','Int','none'],default = 'IntPM')
    parser.add_argument("--nucleiRegion",choices=['watershedContourDist','watershedContourInt','watershedBWDist','dilation'], default = 'watershedContourInt')
    parser.add_argument("--segmentCytoplasm",choices = ['segmentCytoplasm','ignoreCytoplasm'],default = 'segmentCytoplasm')
    parser.add_argument("--cytoDilation",type = int, default = 5)
    parser.add_argument("--logSigma",type = int, nargs = '+', default = [3, 60])
    parser.add_argument("--CytoMaskChan",type=int, nargs = '+', default=[2])
    parser.add_argument("--TissueMaskChan",type=int, nargs = '+', default=[1])
    parser.add_argument("--saveMask",action='store_false')
    parser.add_argument("--saveFig",action='store_false')
    args = parser.parse_args()
    
    # gather filename information
    imagePath = 'Y:/sorger/data/Broad/H1_775/registration/H1_775.ome.tif'
    outputPath = 'Y:/sorger/data/Broad/H1_775/prob_maps/segmentation'
    nucleiClassProbPath = 'Y:/sorger/data/Broad/H1_775/prob_maps/H1_775_ContoursPM_5.tif'
    contoursClassProbPath = 'Y:/sorger/data/Broad/H1_775/prob_maps/H1_775_ContoursPM_5.tif'
    maskPath = 'D:/LSP/cycif/testsets/exemplar-001/dearray/masks/A1_mask.tif'
#    imagePath = args.imagePath
#    outputPath = args.outputPath
#    nucleiClassProbPath = args.nucleiClassProbPath
#    contoursClassProbPath = args.contoursClassProbPath
#    maskPath = args.maskPath
       
    fileName = os.path.basename(imagePath)
    filePrefix = fileName[0:fileName.index('.')]
    
    # get channel used for nuclei segmentation
    if args.probMapChan==-1:
        test = os.path.basename(contoursClassProbPath)
        nucMaskChan = int(test.split('ContoursPM_')[1].split('.')[0])-1
        
    else:
        nucMaskChan = args.probMapChan
    
    #crop images if needed
    if args.crop == 'interactiveCrop':
        nucleiCrop = tifffile.imread(imagePath,key = nucMaskChan)
        r=cv2.selectROI(resize(nucleiCrop,(nucleiCrop.shape[0] // 10, nucleiCrop.shape[1] // 10)))
        cv2.destroyWindow('select')
        rect=np.transpose(r)*10
        PMrect= rect
        nucleiCrop = nucleiCrop[int(rect[1]):int(rect[1]+rect[3]), int(rect[0]):int(rect[0]+rect[2])]
    elif args.crop == 'noCrop' or args.crop == 'dearray' :
        nucleiCrop = tifffile.imread(imagePath,key = nucMaskChan)
        rect = [0, 0, nucleiCrop.shape[0], nucleiCrop.shape[1]]
        PMrect= rect
        
    nucleiProbMaps = tifffile.imread(nucleiClassProbPath,key=0)
    nucleiPM = nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
    nucleiProbMaps = tifffile.imread(contoursClassProbPath,key=0)
    PMSize = nucleiProbMaps.shape
    nucleiPM = np.dstack((nucleiPM,nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]))
        
    # mask the core/tissue
    if args.crop == 'dearray':
        TMAmask = tifffile.imread(maskPath)
    else:
        tissue = np.empty((len(args.TissueMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)
        count=0
        if args.crop == 'noCrop':
            for iChan in args.TissueMaskChan:
                tissue[count,:,:] =tifffile.imread(imagePath,key=iChan)
                count+=1
        else:
            for iChan in args.TissueMaskChan:
                tissue = np.dstack(tissue,normI(tifffile.imread(imagePath,key=iChan)))
        tissueCrop = np.sum(tissue,axis = 0)
        tissue_gauss = gaussian_filter(tissueCrop,1)
        tissue_gauss1 = tissue_gauss.astype(float)
#        tissue_gauss1[tissue_gauss>np.percentile(tissue_gauss,99)]=np.nan
        TMAmask = np.log2(tissue_gauss+1)>threshold_otsu(np.log2(tissue_gauss+1))
        #imshow(TMAmask)
        del tissue_gauss, tissueCrop, tissue, tissue_gauss1
    
    
    # nuclei segmentation
    nucleiMask = S3NucleiSegmentationWatershed(nucleiPM,nucleiCrop,args.logSigma,TMAmask,args.nucleiFilter,args.nucleiRegion)
    del nucleiPM
    
    # cytoplasm segmentation
    if args.segmentCytoplasm == 'segmentCytoplasm':
        count =0
        if args.crop == 'noCrop' or args.crop == 'dearray':
            cyto=np.empty((len(args.CytoMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)    
            for iChan in args.CytoMaskChan:
                cyto[count,:,:] =  tifffile.imread(imagePath, key=iChan)
                count+=1
        else:
            cyto=np.empty((len(args.CytoMaskChan),rect[3],rect[2]),dtype=np.int16)
            for iChan in args.CytoMaskChan:
                cytoFull= tifffile.imread(imagePath, key=iChan)
                cyto[count,:,:] = cytoFull[int(rect[0]):int(rect[0]+rect[2]), int(rect[1]):int(rect[1]+rect[3])]
                count+=1
        cyto = np.amax(cyto,axis=0)
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,args.cytoMethod,args.cytoDilation)
        exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,filePrefix,'nuclei',args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cyto',args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cell',args.saveFig,args.saveMask)
        
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,'ring',args.cytoDilation)
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'nucleiRing',args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cytoRing',args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cellRing',args.saveFig,args.saveMask)
        
    elif args.segmentCytoplasm == 'ignoreCytoplasm':
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'nuclei')
        cellMask = nucleiMask
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'cell')
        cytoplasmMask = nucleiMask
        
        #fix bwdistance watershed
   
                