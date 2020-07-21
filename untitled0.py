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
from skimage.filters import threshold_otsu
from skimage.filters import gaussian
from skimage.feature import peak_local_max
from skimage.color import label2rgb
from skimage.io import imsave, imread
from skimage.segmentation import clear_border
from scipy.ndimage.filters import uniform_filter
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
import copy


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

def S3NucleiSegmentationWatershed(nucleiPM,nucleiImage,logSigma,TMAmask,nucleiFilter,nucleiRegion):
    nucleiContours = nucleiPM[:,:,1]
    nucleiCenters = nucleiPM[:,:,0]
    
    del nucleiPM
    mask = resize(TMAmask,(nucleiImage.shape[0],nucleiImage.shape[1]),order = 0)>0
    nucleiContours = nucleiContours*mask
    if len(logSigma)==1:
         nucleiDiameter  = [logSigma*0.5, logSigma*1.5]
    else:
         nucleiDiameter = logSigma
    logMask = nucleiCenters > 150
    gf=gaussian_filter(255-nucleiContours,logSigma[1]/30)
    gf= extrema.h_maxima(gf,logSigma[1]/30)
    fgm=peak_local_max(gf, indices=False,footprint=np.ones((3, 3)))
    _, fgm= cv2.connectedComponents(fgm.astype(np.uint8))
    foregroundMask= morphology.watershed(nucleiContours,fgm,watershed_line=True)
    del fgm
    allNuclei = ((foregroundMask)*mask)
    del foregroundMask
    if nucleiFilter == 'IntPM':
        P = regionprops(allNuclei,nucleiCenters,cache=False)
    elif nucleiFilter == 'Int':
        P = regionprops(allNuclei,nucleiImage,cache=False)
    mean_int = np.array([prop.mean_intensity for prop in P]) 
    #kmeans = KMeans(n_clusters=2).fit(mean_int.reshape(-1,1))
    MITh = threshold_otsu(mean_int.reshape(-1,1))
    allNuclei = np.zeros(mask.shape,dtype=np.uint32)
    count = 0
    del nucleiImage
    
    for props in P:
        intensity = props.mean_intensity
        if intensity >MITh:
            count += 1
            yi = props.coords[:, 0]
            xi = props.coords[:, 1]
            allNuclei[yi, xi] = count
    
    P = regionprops(allNuclei,cache=False)
    count=0
    maxArea = (logSigma[1]**2)*3/4
    minArea = (logSigma[0]**2)*3/4
    nucleiMask = np.zeros(mask.shape,dtype=np.uint32)
    for props in P:
        area = props.area
        solidity = props.solidity
        if minArea < area < maxArea and solidity >0.8:
            count += 1
            yi = props.coords[:, 0]
            xi = props.coords[:, 1]
            nucleiMask[yi, xi] = count
    return nucleiMask
    
#    img2 = nucleiImage.copy()
#    stacked_img = np.stack((img2,)*3, axis=-1)
#    stacked_img[X > 0] = [65535, 0, 0]
#    imshowpair(img2,stacked_img)

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
        imsave(outputPath + os.path.sep + fileName + 'Mask.tif',mask, plugin="tifffile")
        
    if saveFig== True:
        mask=np.uint8(mask>0)
        edges=cv2.Canny(mask,0,1)
        stacked_img=np.stack((np.uint16(edges)*np.amax(image),image),axis=0)
        tifffile.imsave(outputPath + os.path.sep + fileName + 'Outlines.tif',stacked_img)
        
def S3punctaDetection(spotChan,mask,sigma,SD):
    Ilog = -gaussian_laplace(np.float32(spotChan),sigma)
    fgm=peak_local_max(Ilog, indices=False,footprint=np.ones((3, 3)))
    test=Ilog[fgm==1]
    med = np.median(test)
    mad = np.median(np.absolute(test - med))
    thresh = med + 1.4826*SD*mad
    return (Ilog>thresh)*fgm*(mask>0)
    
#    stacked_img=np.stack((spots*4095,nucleiCrop),axis=0)
#    tifffile.imsave('D:/Seidman/ZeissTest Sets/registration/spottest.tif',stacked_img)
    
    
        
    # assign nan to tissue mask


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
    parser.add_argument("--nucleiRegion",choices=['watershedContourDist','watershedContourInt','watershedBWDist','dilation','localThreshold'], default = 'watershedContourInt')
    parser.add_argument("--segmentCytoplasm",choices = ['segmentCytoplasm','ignoreCytoplasm'],default = 'ignoreCytoplasm')
    parser.add_argument("--cytoDilation",type = int, default = 5)
    parser.add_argument("--logSigma",type = int, nargs = '+', default = [3, 60])
    parser.add_argument("--CytoMaskChan",type=int, nargs = '+', default=[1])
    parser.add_argument("--TissueMaskChan",type=int, nargs = '+', default=-1)
    parser.add_argument("--detectPuncta",type=int, nargs = '+', default=[-1])
    parser.add_argument("--punctaSigma", nargs = '+', default=[1])
    parser.add_argument("--punctaSD", nargs = '+', default=[4])
    parser.add_argument("--saveMask",action='store_false')
    parser.add_argument("--saveFig",action='store_false')
    args = parser.parse_args()
    
    # gather filename information
    #exemplar001
#    imagePath = 'D:/LSP/cycif/testsets/exemplar-001/registration/exemplar-001.ome.tif'
#    outputPath = 'D:/LSP/cycif/testsets/exemplar-001/segmentation'
##    nucleiClassProbPath = 'D:/LSP/cycif/testsets/exemplar-001/probmaps/exemplar-001_NucleiPM_1.tif'
##    contoursClassProbPath = 'D:/LSP/cycif/testsets/exemplar-001/probmaps/exemplar-001_ContoursPM_1.tif'
#    contoursClassProbPath =''
#    stackProbPath = 'D:/LSP/cycif/testsets/exemplar-001/probability_maps/exemplar-001_Probabilities_1.tif'
#    maskPath = 'D:/LSP/cycif/testsets/exemplar-001/dearray/masks/A1_mask.tif'
#    args.cytoMethod = 'hybrid'
#    args.crop = 'autoCrop'
#    args.segmentCytoplasm = 'segmentCytoplasm'
	
    #punctatest
    imagePath = 'D:/Seidman/ZeissTest Sets/registration/13042020_15AP_FAP488_LINC550_DCN647_WGA_40x_1.ome.tif'
    outputPath = 'D:/Seidman/ZeissTest Sets/segmentation'
    nucleiClassProbPath = 'D:/Seidman/ZeissTest Sets/probability-maps/13042020_15AP_FAP488_LINC550_DCN647_WGA_40x_1_NucleiPM_1.tif'
    contoursClassProbPath = 'D:\Seidman\ZeissTest Sets\probability-maps/13042020_15AP_FAP488_LINC550_DCN647_WGA_40x_1_ContoursPM_1.tif'
#    contoursClassProbPath =''
#    stackProbPath = 'D:/LSP/cycif/testsets/exemplar-001/probability_maps/exemplar-001_Probabilities_1.tif'
    maskPath = 'D:/Seidman/ZeissTest Sets/segmentation/13042020_15AP_FAP488_LINC550_DCN647_WGA_40x_1/cellMask.tif'
    args.nucleiRegion = 'localThresh'
    args.crop = 'plate'
    args.logSigma = [6, 300]
    args.segmentCytoplasm = 'ignoreCytoplasm'
    args.detectPuncta = [0,1]
    args.punctaSigma = [1]
    args.punctaSD = [10]
    
    
	    #exemplar002
#    imagePath = 'D:/LSP/cycif/testsets/exemplar-002/dearray/1new.tif'
#    outputPath = 'D:/LSP/cycif/testsets/exemplar-002/segmentation'
#    nucleiClassProbPath =''# D:/LSP/cycif/testsets/exemplar-002/probability-maps/1_NucleiPM_1.tif'
#    contoursClassProbPath = ''#D:/LSP/cycif/testsets/exemplar-002/probability-maps/1_ContoursPM_1.tif'
#    stackProbPath = 'D:/LSP/cycif/testsets/exemplar-002/probability-maps/ilastik_1new_Probabilities_1.tif'
#    maskPath = 'D:/LSP/cycif/testsets/exemplar-002/dearray/masks/1_mask.tif'
#    args.cytoMethod = 'hybrid'
#    args.mask = 'TMA'
#    args.crop = 'dearray'
#    args.segmentCytoplasm = 'ignoreCytoplasm'
		
		#test
#    imagePath = 'Y:\\sorger\\data\\RareCyte\\Zoltan\\CY_200617\\LCN241_076\\dearray\\17.tif'
#    outputPath = 'Y:\\sorger\\data\\RareCyte\\Zoltan\\CY_200617\\LCN241_076\\segmentation'
#    nucleiClassProbPath ='Y:\\sorger\\data\\RareCyte\\Zoltan\\CY_200617\\LCN241_076\\probability-maps\\unmicst\\17_NucleiPM_1.tif'
#    contoursClassProbPath = 'Y:\\sorger\\data\\RareCyte\\Zoltan\\CY_200617\\LCN241_076\\probability-maps\\unmicst\\17_ContoursPM_1.tif'
##    stackProbPath = 'D:/LSP/cycif/testsets/exemplar-002/probability-maps/1_Probabilities_1.tif'
#    maskPath = 'Y:\\sorger\\data\\RareCyte\\Zoltan\\CY_200617\\LCN241_076\\dearray\\masks\\17_mask.tif'
#    args.cytoDilation = 3
##    args.mask = 'TMA'
#    args.crop = 'dearray'
#    args.segmentCytoplasm = 'segmentCytoplasm'	
	
	
	#plate 
#    imagePath = 'Y:/sorger/data/computation/Jeremy/caitlin-ddd-cycif-registered/Plate1/E3_fld_1/registration/E3_fld_1.ome.tif'
#    outputPath = 'Y:/sorger/data/computation/Jeremy/caitlin-ddd-cycif-registered/Plate1/E3_fld_1/segmentation'
#    nucleiClassProbPath = 'Y:/sorger/data/computation/Jeremy/caitlin-ddd-cycif-registered/Plate1/E3_fld_1/prob_maps/E3_fld_1_NucleiPM_1.tif'
#    contoursClassProbPath = 'Y:/sorger/data/computation/Jeremy/caitlin-ddd-cycif-registered/Plate1/E3_fld_1/prob_maps/E3_fld_1_ContoursPM_1.tif'
#    maskPath = 'D:/LSP/cycif/testsets/exemplar-001/dearray/masks/A1_mask.tif'
#    args.crop = 'plate'
#    args.cytoMethod ='hybrid'
        
    #large tissue
#    imagePath =  'Y:/sorger/data/RareCyte/Zoltan/Z174_lung/stitched/7.ome.tif'
#    outputPath = 'D:/LSP/cycif/testsets/exemplar-001/segmentation'
#    nucleiClassProbPath = 'Y:/sorger/data/RareCyte/Zoltan/Z174_lung/unmist/7_NucleiPM_46.tif'
#    contoursClassProbPath = 'Y:/sorger/data/RareCyte/Zoltan/Z174_lung/unmist/7_ContoursPM_46.tif'
#    maskPath = 'D:/LSP/cycif/testsets/exemplar-001/dearray/masks/A1_mask.tif'
#    args.crop = 'interactiveCrop'
    
#    imagePath = args.imagePath
#    outputPath = args.outputPath
#    nucleiClassProbPath = args.nucleiClassProbPath
#    contoursClassProbPath = args.contoursClassProbPath
#    stackProbPath = args.stackProbPath
#    maskPath = args.maskPath
       
    fileName = os.path.basename(imagePath)
    filePrefix = fileName[0:fileName.index('.')]
    
    # get channel used for nuclei segmentation
    if args.probMapChan==-1:
        if len(contoursClassProbPath)>0:
            legacyMode = 1
            probPrefix = os.path.basename(contoursClassProbPath)
            nucMaskChan = int(probPrefix.split('ContoursPM_')[1].split('.')[0])-1
        elif len(stackProbPath)>0:
            legacyMode = 0
            probPrefix = os.path.basename(stackProbPath)
            nucMaskChan = int(probPrefix.split('Probabilities_')[1].split('.')[0])-1
        else:
            print('NO PROBABILITY MAP PROVIDED')
        
    else:
        nucMaskChan = args.probMapChan

    if args.TissueMaskChan==-1:
        TissueMaskChan = copy.copy(args.CytoMaskChan)
        TissueMaskChan.append(nucMaskChan)
    else:
        TissueMaskChan = args.TissueMaskChan[:]
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
        nucleiPM = nucleiProbMaps[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3]),0:2]
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

 #       tissue_gauss = tissueCrop
#        tissue_gauss1 = tissue_gauss.astype(float)
#        tissue_gauss1[tissue_gauss>np.percentile(tissue_gauss,99)]=np.nan
#        TMAmask = np.log2(tissue_gauss+1)>threshold_otsu(np.log2(tissue_gauss+1))
        #imshow(TMAmask)
        del tissue_gauss, tissue

    # nuclei segmentation
    nucleiMask = S3NucleiSegmentationWatershed(nucleiPM,nucleiCrop,args.logSigma,TMAmask,args.nucleiFilter,args.nucleiRegion)
    del nucleiPM
    # cytoplasm segmentation
    if args.segmentCytoplasm == 'segmentCytoplasm':
        count =0
        if args.crop == 'noCrop' or args.crop == 'dearray' or args.crop == 'plate':
            cyto=np.empty((len(args.CytoMaskChan),nucleiCrop.shape[0],nucleiCrop.shape[1]),dtype=np.uint16)    
            for iChan in args.CytoMaskChan:
                cyto[count,:,:] =  skio.imread(imagePath, key=iChan)
                count+=1
        elif args.crop =='autoCrop':
            cyto=np.empty((len(args.CytoMaskChan),int(rect[2]),int(rect[3])),dtype=np.int16)
            for iChan in args.CytoMaskChan:
                cytoFull= skio.imread(imagePath, key=iChan)
                cyto[count,:,:] = cytoFull[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
                count+=1
        else:
            cyto=np.empty((len(args.CytoMaskChan),rect[3],rect[2]),dtype=np.int16)
            for iChan in args.CytoMaskChan:
                cytoFull= skio.imread(imagePath, key=iChan)
                cyto[count,:,:] = cytoFull[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
                count+=1
        cyto = np.amax(cyto,axis=0)
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,args.cytoMethod,args.cytoDilation)
        exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,filePrefix,'nuclei',args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cyto',args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cell',args.saveFig,args.saveMask)
  
        cytoplasmMask,nucleiMaskTemp,cellMask = S3CytoplasmSegmentation(nucleiMask,cyto,TMAmask,'ring',args.cytoDilation)
        exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,filePrefix,'nucleiRing',args.saveFig,args.saveMask)
        exportMasks(cytoplasmMask,cyto,outputPath,filePrefix,'cytoRing',args.saveFig,args.saveMask)
        exportMasks(cellMask,cyto,outputPath,filePrefix,'cellRing',args.saveFig,args.saveMask)
        
    elif args.segmentCytoplasm == 'ignoreCytoplasm':
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'nuclei')
        cellMask = nucleiMask
        exportMasks(nucleiMask,nucleiCrop,outputPath,filePrefix,'cell')
        cytoplasmMask = nucleiMask
        
        
    if (np.min(args.detectPuncta)>-1):
        if len(args.detectPuncta) != len(args.punctaSigma):
            args.punctaSigma = args.punctaSigma[0] * np.ones(len(args.detectPuncta))
 
  
        if len(args.detectPuncta) != len(args.punctaSD):
            args.punctaSD = args.punctaSD[0] * np.ones(len(args.detectPuncta))
  
        counter=0
        for iPunctaChan in args.detectPuncta:
            punctaChan = tifffile.imread(imagePath,key = iPunctaChan)
            punctaChan = punctaChan[int(PMrect[0]):int(PMrect[0]+PMrect[2]), int(PMrect[1]):int(PMrect[1]+PMrect[3])]
            spots=S3punctaDetection(punctaChan,cellMask,args.punctaSigma[counter],args.punctaSD[counter])
            cellspotmask = tifffile.imread(maskPath) #REMOVE THIS LATER
            P = regionprops(cellspotmask,intensity_image = spots ,cache=False)
            numSpots = []
            for prop in P:
                numSpots.append(np.uint16(np.round(prop.mean_intensity * prop.area)))
            np.savetxt(outputPath + 'numSpots.csv',(np.transpose(np.asarray(numSpots))),fmt='%10.5f')    
            edges = 1-(cellMask>0)
            stacked_img=np.stack((np.uint16((spots+edges)>0)*np.amax(punctaChan),punctaChan),axis=0)
            tifffile.imsave(outputPath + os.path.sep + filePrefix + os.path.sep + 'punctaChan'+str(iPunctaChan) + 'Outlines.tif',stacked_img)
            counter=counter+1

#fix bwdistance watershed