function S3segmenterWrapper(testPath,fileName,p)
%this function requires input of nuclei stack range. It assumes that every
%stack beyond that to the end is a cytoplasmic stain. Marker controlled
%watershed based on distance transform of nuclei channel is employed to
%separate nuclei clumps.

% training with nuclei model = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\nucleiTrainingSet\halfsize','sigmas',[3 7 15 30],'sfSigmas',[],'edgeSigmas',1,'ridgeSigmas',1,'radii',[3 5 10 ],'logSigmas',[4 6 8],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei model = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\Real experiment\All Cyles\nucleiTrainingSet10x','sigmas',[1 3 7 9],'sfSigmas',[],'edgeSigmas',2,'ridgeSigmas',2,'radii',[3 5 7 ],'logSigmas',[2 4 6 ],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);

% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\halfsize','sigmas',[1 3 7 ],'sfSigmas',[1 2 3],'logSigmas',[2 4],'nhoodStd',[3 5 11],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\RFTrainingSet','sigmas',[1 3 5 7  ],'sfSigmas',[],'edgeSigmas',[2],'ridgeSigmas',[1 2],'logSigmas',[1 2 4 ],'nhoodStd',[5 7 11 ],'pctMaxNPixelsPerLabel',40,'adjustContrast',0);
% training with nuclei modelCatSF1 = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet','sigmas',[1 3 7 ],'sfSigmas',[1],'logSigmas',[1 2 4 6],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\RFTrainingSet','sigmas',[1 3 5 7 ],'sfSigmas',[1 2],'logSigmas',[1 2 4 6],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',50,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\manualRFTrainingSet','sigmas',[1 2 4 ],'sfSigmas',[1 2],'edgeSigmas',[1 2],'ridgeSigmas',[1 2],'logSigmas',[1 2 4],'nhoodStd',[3 5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);


% override default values because using Docker config file
modelPathName = p.modelPath;


%% extract image properties
imagePath = p.paths.registration;
if nargin < 1 
    if nargin <2 
         testPath = pwd;
    end
     [fileName, testPath] = uigetfile([testPath filesep '*.tif'],'Select file to process');
else
    if ~isequal(testPath(end),filesep)
        testPath = [testPath filesep];
    end
end
filePrefix = fileName(1:strfind(fileName,'.')-1);
paths=regexp(p.paths.samplefolder,filesep,'split');
sampleName = paths{end-1};
rawFileListing = fileName;
metadata =bfGetReader([imagePath rawFileListing]);            
numChan =metadata.getImageCount;

 %% set up output directories
if isequal(p.Docker,'true')
    outputPath = p.dockerParams.outputPath;
else
    outputPath = p.paths.segmentation;
    mkdir(outputPath)
end

name=fileName(1:strfind(fileName,'.')-1);
mkdir([outputPath name])
outputPath = [outputPath name];  

    %% read class probability maps
    
    switch p.ClassProbSource
        case 'RF'
            if isempty(modelPathName)
                [modelFileName, modelPathName] = uigetfile([testPath filesep '*.mat'],'Select RF model');
            end
            % modelFileName = 'modelReducedFeatures.mat';
            % modelPathName = 'D:\LSP\cycif\S3\nucleiTrainingSet\'; 
        %#pragma treeBagger
            load([modelPathName modelFileName])

            nuclei = imread([testPath rawFileListing]);
            [nucleiPM,nucleiCrop]=RFSplitter(nuclei,modelNuc,'split',false);

        case 'unet'
            classProbPath = p.paths.probabilitymaps;
            listing = dir([classProbPath name '_ContoursPM*']);
            PMfileName = listing.name;
            pmI = strfind(PMfileName,'_ContoursPM_');
            probMapSuffix= cellstr(PMfileName(pmI:end));
            nucMaskChan = sscanf(char(probMapSuffix), '_ContoursPM_%d.tif');
            if nucMaskChan >numChan
                nucMaskChan = nucMaskChan - numChan;
            end
               
            nucleiPMListing = dir([classProbPath name '_NucleiPM*']);
            if ~isempty(nucleiPMListing)
                PMfileName = nucleiPMListing.name;
                pmI = strfind(PMfileName,'_NucleiPM_');
                probMapSuffix{2}= PMfileName(pmI:end);
            end
            
        case 'none'
            classProbPath = p.paths.registration;
            probMapSuffix= '.ome.tif';
        
    end 
    
    nucleiCrop = imread([imagePath rawFileListing],nucMaskChan);
    fullResSize = size(nucleiCrop);
    
    
    %% read nuclei channel and crop as necessary
        if exist([classProbPath filePrefix probMapSuffix{1}])
            classProbInfo = imfinfo([classProbPath filePrefix probMapSuffix{1}]);
            
            switch p.crop
                case 'interactiveCrop'
                 [~,rect] = imcrop(nucleiCrop);
                 rect = [round(rect(1)) round(rect(2)) round(rect(3)) round(rect(4))];
                 PMrect= rect; 
                 nucleiCrop = imread([imagePath rawFileListing],nucMaskChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]});
                 
                case 'autoCrop'
                 rect = [round(metadata.getSizeX/3) round(metadata.getSizeY/3) round(metadata.getSizeX/3) round(metadata.getSizeY/3)];
                 PMrect= rect;
                 nucleiCrop = imread([imagePath rawFileListing],nucMaskChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]});
                case 'dearray'
                    TMAcropInfo = imfinfo([p.paths.dearray filesep name '.tif' ]);
                    rawSize = [TMAcropInfo.Height TMAcropInfo.Width];
                    if exist([p.paths.dearray filesep name '_cropCoords.mat'])
                        load([p.paths.dearray filesep name '_cropCoords.mat'])
                        rect = [rect(1) rect(2) rect(3)-rect(1)-1 rect(4)-rect(2)-1];
                        PMrect =  [1 1 rect(3) rect(4)];
                        if (rect(3)+rect(1)) > fullResSize(2) || (rect(4)+rect(2)) > fullResSize(1)
                            rect = [1 1 fullResSize(2)-1 fullResSize(1)-1];
                            PMrect= rect;
                        end
                    else
                        rect = [1 1 rawSize(2)-1 rawSize(1)-1];
                        PMrect= rect;
                    end
                    nucleiCrop = imread([imagePath rawFileListing],nucMaskChan,'PixelRegion',{[rect(2) rect(2)+rect(4)], [rect(1) rect(1)+rect(3)]});
                    fullResSize = size(nucleiCrop);
                case 'noCrop'
                 nucleiCrop = imresize(nucleiCrop,p.resizeFactor);
                 rect = round([1 1 size(nucleiCrop,2) size(nucleiCrop,1)]);
                 PMrect= rect;
            end

            nucleiPM=[];
             for iPM = 1:numel(probMapSuffix)
                nucleiProbMaps = imread([classProbPath filePrefix probMapSuffix{iPM}],1);
                PMSize = size(nucleiProbMaps);
                nucleiProbMaps = imresize(nucleiProbMaps,fullResSize);
                nucleiPM(:,:,iPM) = imcrop(nucleiProbMaps,PMrect);
                
             end
             PMUpsampleFactor = fullResSize(1)/PMSize(1);
        else
            disp([filePrefix contourprobMapSuffix ' not found'])
            return
        end


%% mask the core/tissue

    
if isequal(p.crop,'dearray')
    if exist([p.paths.dearray 'masks' filesep filePrefix '_mask.tif'])
        TMAmask = imread([p.paths.dearray 'masks' filesep filePrefix '_mask.tif'])>0;
    else
        TMAmask = coreSegmenterFigOutput(max(tissue,[],3),'activeContours','true','split','false');
    end
else
    tissue =[];
    if isequal(p.crop,'noCrop')
        for iChan = p.TissueMaskChan
            tissue= cat(3,tissue,normI(double(imresize(imread([imagePath rawFileListing],iChan),p.resizeFactor))));
        end
    else
        for iChan = p.TissueMaskChan
            tissue= cat(3,tissue,normI(double(imresize(imread([imagePath rawFileListing],iChan,...
                'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}),p.resizeFactor))));
        end
    end
    tissueCrop = sum(tissue,3);
    tissue_gauss = imgaussfilt3(imresize(tissueCrop,0.5),1);
    TMAmask=imresize(tissue_gauss>thresholdMinimumError(tissue_gauss,'model','poisson'),size(tissueCrop));
    if sum(sum(TMAmask)) ==0
        return
    end
    
    
    if p.RefineTissueMask >0
        dsfactor =0.02;
        stats=regionprops(TMAmask,'area');
        smallMask = bwareafilt(TMAmask,[0 prctile(cat(1,stats.Area),99.99)]);
        tissue_gauss1 = imresize(tissue_gauss,dsfactor,'nearest');
        bw = activecontour(tissue_gauss1, ones(size(tissue_gauss1)), 500,'Chan-Vese');
        bw = bwareaopen(bw,100);
        dist=bwdist(~bw);
        distMask = imresize(dist>p.RefineTissueMask/50,size(smallMask),'nearest');
        TMAmask = ((smallMask + distMask)>0).*TMAmask;
    end
        
    
end
clear tissue_gauss, clear maxTissue, clear tissueCrop, clear tissue, clear distMask
   %% nuclei segmentation
   if isequal(p.segmentNucleus,'true')
%       [nucleiMask,largestNucleiArea] = S3NucleiSegmentationMultiSize(nucleiPM,nucleiCrop,p.logSigma,'mask',TMAmask,'inferNucCenters',p.inferNucCenters,...
%             'nucleiFilter',p.nucleiFilter,'nucleiRegion',p.nucleiRegion,'resize',p.resizeFactor,'LoGresize',PMSizeFactor);
%         [nucleiMask,largestNucleiArea] = S3NucleiSegmentation(nucleiPM,nucleiCrop,p.logSigma,'mask',TMAmask,'inferNucCenters',p.inferNucCenters,...
%             'nucleiFilter',p.nucleiFilter,'nucleiRegion',p.nucleiRegion,'resize',p.resizeFactor,'LoGresize',PMSizeFactor);
        [nucleiMask,largestNucleiArea] = S3NucleiSegmentationWatershed(nucleiPM,nucleiCrop,p.logSigma,'mask',TMAmask,'inferNucCenters',p.inferNucCenters,...
            'nucleiFilter',p.nucleiFilter,'nucleiRegion',p.nucleiRegion,'resize',1);
   else
       listing = dir([p.paths.segmentation filePrefix filesep '*nucleiMask.tif']);
       if isempty(listing)
           disp(['No nuclei segmentation mask found!'])  
           return
       else
           nucleiMask = imread([listing.folder filesep listing.name]);
           largestNucleiArea=100000;
       end
   end
   
   if sum(sum(nucleiMask))==0
       return
   else       
    disp(['Segmented Nuclei'])  
   end

clear nucleiPM
clear nuclei

    %% cytoplasm segmentation
   switch p.segmentCytoplasm
       case 'segmentCytoplasm'
            cyto =[];
            for iChan = p.CytoMaskChan
                if isequal(p.crop,'noCrop')
                   cyto= cat(3,cyto,normI(imresize(double(imread([imagePath rawFileListing],iChan)),p.resizeFactor)));
                else
                    cyto= cat(3,cyto,normI(double(imread([imagePath rawFileListing],iChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}))));
                end
            end
            
            cyto=max(cyto,[],3);
           
            % load random forest model if 
            if isequal(p.cytoMethod,'RF')
                %#pragma treeBagger
                load('D:\LSP\cycif\S3\cytoTrainingSet\modelContours1.mat')
            else
            modelCat=[];
            end
%             cyto=S3tileReturn(cyto);
%             tissue = S3tileReturn(tissue);
            [cytoplasmMask,nucleiMaskTemp,cellMask]=S3CytoplasmSegmentation(imdilate(nucleiMask,strel('disk',3)),cyto,modelCat,'mask',TMAmask,...
                'cytoMethod',p.cytoMethod,'resize',1,'sizeFilter',largestNucleiArea,'upSample',p.upSample);
            
            exportMasks(imerode(nucleiMask>0,strel('disk',2)),nucleiCrop,outputPath,'nuclei',p.saveFig,p.saveMasks)
            exportMasks(cytoplasmMask,cyto,outputPath,'cyto',p.saveFig,p.saveMasks)
%             exportMasks(cellMask,cyto,outputPath,'cell',p.saveFig,p.saveMasks)

%             [cytoplasmMaskRing,nucleiMaskRing,cellMaskRing]=S3CytoplasmSegmentation(imdilate(nucleiMask,strel('disk',3)),cyto,modelCat,'mask',TMAmask,...
%                 'cytoMethod','ring','resize',1,'sizeFilter',largestNucleiArea,'upSample',p.upSample);
%             exportMasks(imerode(nucleiMaskRing>0,strel('disk',3)),nucleiCrop,outputPath,'nucleiRing',p.saveFig,p.saveMasks)
%             exportMasks(cytoplasmMaskRing,cyto,outputPath,'cytoRing',p.saveFig,p.saveMasks)
%             exportMasks(cellMaskRing,cyto,outputPath,'cellRing',p.saveFig,p.saveMasks)
%             
        case 'loadMask'
            listing = dir([p.paths.segmentation filePrefix filesep '*cytoplasmmask.tif']);
            if isempty(listing)
                disp(['No cytoplasm segmentation mask found!'])  
                return
            else
                cytoplasmMask = imread([listing.folder filesep listing.name]);
                nucleiMask = cast(nucleiMask>0,class(cytoplasmMask)).*cytoplasmMask;
            end
                        
        case 'ignoreCytoplasm'
            nucleiMask = bwlabel(nucleiMask);
            exportMasks(nucleiMask,nucleiCrop,outputPath,'nuclei',p.saveFig,p.saveMasks)
            cellMask=nucleiMask;
            exportMasks(cellMask,nucleiCrop,outputPath,'cell',p.saveFig,p.saveMasks)
            cytoplasmMask = nucleiMask;
    
   end
   disp(['Segmented Cytoplasm'])
  
%  imshowpair(bwperim(cytoplasmMask),imadjust(cyto))
  
    %% measureFeatures
    if isequal(p.measureFeatures,'true')
        registrationFileList = dir([p.paths.registration filesep]);
        mainFile=[];
        for iFolder = 1:length(registrationFileList)
        fName =  registrationFileList(iFolder).name;
            if ~isfolder(fName) && ~contains(fName,'..') && ~isequal(fName,'.')
                mainFile{end+1} = fName;
            end
        end
        [cellID, meanIntRegionTable, medianIntRegionTable,meanStdTable,meanEntTable, haralickTable,majorAxisTable, minorAxisTable, ...
                    eccTable,solidityTable,areaTable, meanLawTable,centroidCellTable,headMeanInt, headMedianInt,headStd, headEnt,headHaralick,headLaw,headShape]=...
        S3MeasureFeatures(cat(3,nucleiMaskTemp,cytoplasmMask),p.paths,fileName,'MedianIntensity',p.MedianIntensity,...
                    'Docker',p.Docker,'crop',rect,'chanRange',p.chanRange,'channelNames',p.channelNames);
        
        if ~isempty(meanIntRegionTable)
                 writetable(array2table([cellID, meanIntRegionTable, medianIntRegionTable,areaTable,centroidCellTable, majorAxisTable, minorAxisTable, eccTable, ...
               solidityTable, meanStdTable,meanEntTable,haralickTable,meanLawTable],...
               'VariableNames',lower(regexprep(['cellid',headMeanInt,headMedianInt,headShape,headStd, headEnt,headHaralick,headLaw], '[()/\-\s+]', ''))),...
                            [p.paths.analysis filesep name '_Features.txt'],'Delimiter','\t')
        end
        
        
        disp(['Measured all features'])
    end
    
    save([outputPath filesep 'segParams.mat'],'p')
    close all
    disp(['Completed ' fileName])
 
  
end


function exportMasks(mask,image,outputPath,fileNameNuclei,saveFig,saveMasks)
   image = im2double(image);
   if isequal(saveFig,'true')
    t = Tiff([outputPath filesep fileNameNuclei 'Outlines.tif'],'w8');
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(t,'Compression',Tiff.Compression.None);
    setTag(t,'BitsPerSample',16);
    setTag(t,'SamplesPerPixel',2);
    setTag(t,'ImageLength',size(mask,1));
    setTag(t,'ImageWidth',size(mask,2));
    setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Separate);
    setTag(t,'RowsPerStrip', 1);    
    write(t,cat(3,uint16(bwperim(mask))*65535,image*65535));
    close(t);
%     TiffWrite(cat(3,uint16(bwperim(mask))*65535,image*65535),[outputPath filesep fileNameNuclei 'Outlines.tif'])
   end
   if isequal(saveMasks,'true')
    t = Tiff([outputPath filesep fileNameNuclei 'Mask.tif'],'w8');
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(t,'Compression',Tiff.Compression.None);
    setTag(t,'BitsPerSample',32);
    setTag(t,'SamplesPerPixel',1);
    setTag(t,'ImageLength',size(mask,1));
    setTag(t,'ImageWidth',size(mask,2));
    setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    setTag(t,'RowsPerStrip', 1);
    write(t,uint32(mask));
    close(t);
%     TiffWrite(mask,[outputPath filesep fileNameNuclei 'Mask.tif'])
   end
end

