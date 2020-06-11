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
sampleName =paths{end-1};
rawFileListing = dir([p.paths.registration sampleName '*.' p.paths.fileExt]);
metadata =bfGetReader([imagePath rawFileListing(1).name]);            
p.numChan =metadata.getImageCount;

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
    
%     switch p.ClassProbSource
%         case 'RF'
%             if isempty(modelPathName)
%                 [modelFileName, modelPathName] = uigetfile([testPath filesep '*.mat'],'Select RF model');
%             end
%             % modelFileName = 'modelReducedFeatures.mat';
%             % modelPathName = 'D:\LSP\cycif\S3\nucleiTrainingSet\'; 
%         %#pragma treeBagger
%             load([modelPathName modelFileName])
% 
%             nuclei = imread([testPath rawFileListing(1).name]);
%             [nucleiPM,nucleiCrop]=RFSplitter(nuclei,modelNuc,'split',false);
% 
%         case 'unet' 
            classProbPath = p.paths.probabilitymaps;
            PMlisting=dir([classProbPath filesep '*_Probabilities*']);
            if numel(PMlisting)>0 % if using a single file containing all PMs
                singleChannelPM=false;
                PMfileName = PMlisting.name;
                pmI = strfind(PMfileName,'_Probabilities');
                nucMaskChan =str2num(PMfileName(strfind(PMfileName,'Probabilities')+14:strfind(PMfileName,'.tif')-1));
                if isempty(nucMaskChan)
                    nucMaskChan = p.NucMaskChan;
                end
                
                if nucMaskChan >p.numChan
                    nucMaskChan = nucMaskChan - p.numChan;
                end
                nucleiprobMapSuffix= cellstr(PMfileName(pmI:end));
                nucleiprobMapSuffix{2}= PMfileName(pmI:end);
                nucleiprobMapSuffix{3}= PMfileName(pmI:end);
            else
                singleChannelPM = true;
                listing = dir([classProbPath name '_ContoursPM*']);
                PMfileName = listing.name;
                pmI = strfind(PMfileName,'_NucleiPM_');
                nucleiprobMapSuffix= cellstr(PMfileName(pmI:end));
                nucMaskChan = sscanf(char(nucleiprobMapSuffix), '_NucleiPM_%d.tif');
                if nucMaskChan >p.numChan
                    nucMaskChan = nucMaskChan - p.numChan;
                end

                nucleiPMListing = dir([classProbPath name '_ContoursPM*']);
                if ~isempty(nucleiPMListing)
                    PMfileName = nucleiPMListing.name;
                    pmI = strfind(PMfileName,'_ContoursPM_');
                    nucleiprobMapSuffix{2}= PMfileName(pmI:end);
                end
                nucleiprobMapSuffix{3}= '_backgroundPM_';
            end
            
            if isequal(p.ClassProbSource,'ilastik')
                 nucleiprobMapSuffix{4}= PMfileName(pmI:end);
            else
                cytoPMListing = dir([classProbPath name '_CytoPM*']);
                if ~isempty(cytoPMListing)
                    PMfileName = cytoPMListing.name;
                    pmI = strfind(PMfileName,'_CytoPM_');
                    cytoprobMapSuffix= PMfileName(pmI:end);
                    p.CytoMaskChan = sscanf(char(cytoprobMapSuffix), '_CytoPM_%d.tif');
                end
            end

%         case 'none'
%             classProbPath = p.paths.registration;
%             probMapSuffix= '.ome.tif';
%         
%     end 
    

    
    
    %% read nuclei channel and crop as necessary
        if exist([classProbPath filePrefix nucleiprobMapSuffix{1}])
            classProbInfo = imfinfo([classProbPath filePrefix nucleiprobMapSuffix{1}]);
            
            switch p.crop
                case 'interactiveCrop'
                 nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan);
                 fullResSize = size(nucleiCrop);
                 [~,rect] = imcrop(nucleiCrop);
                 rect = [round(rect(1)) round(rect(2)) round(rect(3)) round(rect(4))];
                 PMrect= rect; 
                 nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]});
                 
                case 'autoCrop'
                 rect = [round(metadata.getSizeX/3) round(metadata.getSizeY/3) round(metadata.getSizeX/3) round(metadata.getSizeY/3)];
                 PMrect= rect;
                 nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]});
                case 'dearray'
                    TMAcropInfo = imfinfo([p.paths.dearray filesep name '.tif' ]);
                    rawSize = [TMAcropInfo.Height TMAcropInfo.Width];
                    if exist([p.paths.dearray filesep name '_cropCoords.mat'])
                        load([p.paths.dearray filesep name '_cropCoords.mat'])
                        rect = [rect(1) rect(2) rect(3)-rect(1)-1 rect(4)-rect(2)-1];
                        PMrect =  [1 1 rect(3) rect(4)];
                    else
                        rect = [1 1 rawSize(2)-1 rawSize(1)-1];
                        PMrect= rect;
                    end
                    nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan,'PixelRegion',{[rect(2) rect(2)+rect(4)], [rect(1) rect(1)+rect(3)]});
                    fullResSize = size(nucleiCrop);
                case 'noCrop' 
                 nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan);
                 fullResSize = size(nucleiCrop);
                 nucleiCrop = imresize(nucleiCrop,p.resizeFactor);
                 rect = round([1 1 size(nucleiCrop,2) size(nucleiCrop,1)]);
                 PMrect= rect;
                 case 'plate' 
                 nucleiCrop = imread([imagePath rawFileListing(1).name],nucMaskChan);
                 fullResSize = size(nucleiCrop);
                 nucleiCrop = imresize(nucleiCrop,p.resizeFactor);
                 rect = round([1 1 size(nucleiCrop,2) size(nucleiCrop,1)]);
                 PMrect= rect;
            end

            %read nucleiPMs and cytoPM (if any)
            nucleiPM=[];
            cytoPM=[];
            if singleChannelPM==true
                 for iPM = 1:numel(nucleiprobMapSuffix)
                    nucleiProbMaps = imread([classProbPath filePrefix nucleiprobMapSuffix{iPM}],1);
                    PMSize = size(nucleiProbMaps);
                    nucleiProbMaps = imresize(nucleiProbMaps,fullResSize);
                    nucleiPM(:,:,iPM) = imcrop(nucleiProbMaps,PMrect);
                 end
                    cytoProbMaps = imread([classProbPath filePrefix nucleiprobMapSuffix],1);
                    PMSize = size(cytoProbMaps);
                    cytoProbMaps = imresize(cytoProbMaps,fullResSize);
                    nucleiPM(:,:,4) = imcrop(cytoProbMaps,PMrect);
                    del cytoProbMaps
                 
            else
                nucleiProbMaps = imread([classProbPath filePrefix nucleiprobMapSuffix{1}]);
                if isequal(p.ClassProbSource,'unet') || isequal(p.ClassProbSource,'ilastik')
                    if exist('cytoprobMapSuffix','var')==1 
                        nucleiProbMaps =cat(3,nucleiProbMaps,imread([classProbPath filePrefix cytoprobMapSuffix]));
                    end
                    if numel(nucleiprobMapSuffix)==4
                        nucleiProbMaps =cat(3,nucleiProbMaps,nucleiProbMaps(:,:,p.probMapOrder(4)));
                    end
                end
                    for iPM = 1:size(nucleiProbMaps,3)
                        nucleiPM(:,:,iPM) = imcrop(imresize(nucleiProbMaps(:,:,iPM),fullResSize),PMrect);
                    end
                
                    PMSize = size(nucleiProbMaps);
            end
            
            
            
             
             PMUpsampleFactor = fullResSize(1)/PMSize(1);
        else
            disp([filePrefix contourprobMapSuffix ' not found'])
            return
        end
p.rect= rect;

%% mask the core/tissue
if isempty(p.TissueMaskChan)
    p.TissueMaskChan = [nucMaskChan p.CytoMaskChan];
end

if isequal(p.crop,'dearray')
    if exist([p.paths.dearray 'masks' filesep filePrefix '_mask.tif'])
        TMAmask = imread([p.paths.dearray 'masks' filesep filePrefix '_mask.tif'])>0;
    else
        TMAmask = coreSegmenterFigOutput(max(tissue,[],3),'activeContours','true','split','false');
    end
elseif isequal(p.crop,'plate')
    TMAmask = ones(size(nucleiCrop));
    
else
    tissue =[];
    if isequal(p.crop,'noCrop')
        for iChan = p.TissueMaskChan
            tissueCrop= normI(double(imresize(imread([imagePath rawFileListing(1).name],iChan),p.resizeFactor)));
            tissue_gauss = imgaussfilt3(tissueCrop,1);
            tissue_gauss(tissue_gauss==0)=NaN; %remove outlier contributions to automatic threshold calculation
            tissue=cat(3,tissue,tissue_gauss>thresholdMinimumError(tissue_gauss,'model','poisson'));
        end
    else
        for iChan = p.TissueMaskChan
%             tissue= cat(3,tissue,normI(double(imresize(imread([imagePath rawFileListing(1).name],iChan,...
%                 'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}),p.resizeFactor))));
              tissueCrop= double(imread([imagePath rawFileListing(1).name],iChan,...
                'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}));
               tissue_gauss = imgaussfilt3(tissueCrop,1);
               tissue_gauss(tissue_gauss==0)=NaN; %remove outlier contributions to automatic threshold calculation
               tissue=cat(3,tissue,tissue_gauss>thresholdMinimumError(tissue_gauss,'model','poisson'));
               tissue=cat(3,tissue,tissue_gauss>thresholdOtsu(tissue_gauss));
        end
    end
        TMAmask = max(tissue,[],3);
%     tissueCrop = max(tissue,[],3);
%     tissue_gauss = imgaussfilt3(tissueCrop,1);
%     tissue_gauss(tissue_gauss==0)=NaN; %remove outlier contributions to automatic threshold calculation
%     TMAmask=imresize(tissue_gauss>thresholdMinimumError(tissue_gauss,'model','poisson'),size(tissueCrop));
   
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
       disp(['No nuclei found. Consider increasing nuclei size range or choose a different tissue mask channel'])  
       return
   else       
    disp(['Segmented Nuclei'])  
   end


clear nuclei

    %% cytoplasm segmentation
   switch p.segmentCytoplasm
       case 'segmentCytoplasm'
            if (isequal(p.ClassProbSource,'ilastik') || isequal(p.ClassProbSource,'unet')) && size(nucleiPM,3)>3 %use cytoPM
                    p.cytoPM = nucleiPM(:,:,4);
                    cyto = nucleiPM(:,:,4);
                    if isequal(p.ClassProbSource,'ilastik')
                        TMAmask = nucleiPM(:,:,3)<thresholdOtsu(nucleiPM(:,:,3));
                    end
            else 
                p.cytoPM=[];
                cyto =[];
                for iChan = p.CytoMaskChan
                    if isequal(p.crop,'noCrop')
                       cyto= cat(3,cyto,normI(imresize(double(imread([imagePath rawFileListing(1).name],iChan)),p.resizeFactor)));
                    elseif isequal(p.crop,'dearray')
                        cyto= cat(3,cyto,normI(imresize(double(imread([p.paths.dearray filesep fileName],iChan)),p.resizeFactor)));
                    else
                        cyto= cat(3,cyto,normI(double(imread([imagePath rawFileListing(1).name],iChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}))));
                    end
                end
                cyto=max(cyto,[],3);
            end

            [cytoplasmMask,nucleiMaskTemp,cellMask]=S3CytoplasmSegmentation(nucleiMask,cyto,'mask',TMAmask,...
                'cytoMethod',p.cytoMethod,'resize',1,'sizeFilter',largestNucleiArea,'upSample',p.upSample,...
                'cytoDilation',p.cytoDilation,'cytoPM',p.cytoPM);
            exportMasks(nucleiMaskTemp,nucleiCrop,outputPath,'nuclei',p.saveFig,p.saveMasks)
            exportMasks(cytoplasmMask,cyto,outputPath,'cyto',p.saveFig,p.saveMasks)
            exportMasks(cellMask,cyto,outputPath,'cell',p.saveFig,p.saveMasks)
            
            [cytoplasmMaskRing,nucleiMaskRing,cellMaskRing]=S3CytoplasmSegmentation(nucleiMask,cyto,'mask',TMAmask,...
                'cytoMethod','ring','resize',1,'sizeFilter',largestNucleiArea,'upSample',p.upSample,...
                'cytoDilation',p.cytoDilation,'cytoPM',p.cytoPM);
            exportMasks(nucleiMaskRing,nucleiCrop,outputPath,'nucleiRing',p.saveFig,p.saveMasks)
            exportMasks(cytoplasmMaskRing,cyto,outputPath,'cytoRing',p.saveFig,p.saveMasks)
            exportMasks(cellMaskRing,cyto,outputPath,'cellRing',p.saveFig,p.saveMasks)
            
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
  
%  clear nucleiPM
 
  %% detect puncta
if (min(p.detectPuncta)>0) && (numel(p.detectPuncta) <= p.numChan)
      if numel(p.detectPuncta) ~= numel(p.punctaSigma)
          p.punctaSigma = p.punctaSigma(1) * ones(1,numel(p.detectPuncta));
      end
      
      if numel(p.detectPuncta) ~= numel(p.punctaSD)
          p.punctaSD = p.punctaSD(1) * ones(1,numel(p.detectPuncta));
      end
    counter=1;
    for iPunctaChan = p.detectPuncta
        punctaChan =imread([imagePath rawFileListing(1).name],iPunctaChan,'PixelRegion',{[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]}); 
        spots=S3punctaDetection(punctaChan,p.punctaSigma(counter),p.punctaSD(counter));
        exportMasks(spots.*cellMask,punctaChan,outputPath,['punctaChan' int2str(iPunctaChan)],'false','true')
        exportMasks(bwperim(cellMask>0)+(spots.*cellMask),punctaChan,outputPath,['punctaChan' int2str(iPunctaChan)],'true','false')
        exportMasks(spots.*cellMaskRing,punctaChan,outputPath,['punctaRingChan' int2str(iPunctaChan)],'false','true')
        exportMasks(bwperim(cellMaskRing>0)+(spots.*cellMaskRing),punctaChan,outputPath,['punctaRingChan' int2str(iPunctaChan)],'true','false')
        counter=counter+1;
    end
end
  
    
      
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
        S3MeasureFeatures(cat(3,nucleiMaskTemp,cytoplasmMask),p.paths,char(mainFile),'MedianIntensity',p.MedianIntensity,...
                    'Docker',p.Docker,'crop',rect,'chanRange',p.chanRange);
                
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

