function O2batchS3segmenterWrapperR(mainPath,varargin)


ip = inputParser;
ip.addParamValue('HPC','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('fileNum',1,@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('ClassProbSource','unet',@(x)(ismember(x,{'RF','unet','none'})));
ip.addParamValue('NucMaskChan',[2],@(x)(numel(x) > 0 & all(x > 0 )));  % deprecated. Channel number implied from prob map filename
ip.addParamValue('CytoMaskChan',[2],@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('TissueMaskChan',[3],@(x)(numel(x) > 0 & all(x > 0 ))); 
ip.addParamValue('RefineTissueMask',[0],@(x)(numel(x) > 0 & all(x > 0 ))); 
ip.addParamValue('mask','tissue',@(x)(ismember(x,{'TMA','tissue','none'}))); % set to true if sample is TMA cores
ip.addParamValue('crop','noCrop',@(x)(ismember(x,{'interactiveCrop','autoCrop','dearray','noCrop'})));
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','bwdistanceTransform','ring'})));
ip.addParamValue('MedianIntensity','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveFig','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveMasks','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('segmentNucleus','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiFilter','IntPM',@(x)(ismember(x,{'LoG','Int','IntPM','none'})));
ip.addParamValue('measureFeatures','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','watershedBWDist','dilation'})));
ip.addParamValue('segmentCytoplasm','segmentCytoplasm',@(x)(ismember(x,{'segmentCytoplasm','loadMask','ignoreCytoplasm'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('inferNucCenters','UNet',@(x)(ismember(x,{'UNet','RF','Int'})));
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resizeFactor',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('logSigma',[2.5],@(x)(numel(x) >0 & all(x > 0 )));
ip.addParamValue('chanRange',[0],@(x)(numel(x) >0 & all(x > 0 ))); %channels for measuring features. If 0, assume all channels.
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('dockerParams',0,@(x)(numel(x)==1));
ip.parse(varargin{:});          
p = ip.Results;  

if isequal(p.Docker,'true')
    mainPath = p.dockerParams.parentPath;
    p.modelPath = [p.dockerParams.modelPath filesep];
    p.modelCatPath = [p.dockerParams.modelCatPath filesep];
    
else
    if nargin<1 
        mainPath = uigetdir(pwd,'Select a folder with image(s) to process');
    end
    p.modelPath ='';
    p.modelCatPath ='';
end

paths.metadata = ['metadata' ];
paths.dearray = ['dearray' ];
paths.probabilitymaps= ['prob_maps3Class'];
paths.segmentation = ['segmentation3Class'];
paths.analysis = ['analysis'];
paths.registration = [''];

if isequal(p.crop,'dearray')
    FileExt = 'tif';
    searchPath = ['dearray'];
    
else
    FileExt = 'tif';
    searchPath = [''];
end

%% read file names
    sampleFolderList = dir([mainPath filesep '*' searchPath]);
    finalSampleFolderList = [];
    
    for iFolder = 1:length(sampleFolderList)
        if isfolder(sampleFolderList(iFolder).folder) && ~isequal(sampleFolderList(iFolder).name,'..') ...
                && ~isequal(sampleFolderList(iFolder).name,'.') && ~isequal(sampleFolderList(iFolder).name,'metadata')
            finalSampleFolderList{end+1} = sampleFolderList(iFolder).name;
            
        end
    end
    disp (['Found ' num2str(length(finalSampleFolderList)) ' samples(s)!'])
    
    if isequal(p.HPC,'false')
        folderNumStart = 1;
        folderNumEnd = length(finalSampleFolderList);
    else 
        folderNumStart=p.fileNum;
        folderNumEnd = p.fileNum;
    end

    
%% import metadata file
    metadataListing= dir([mainPath filesep 'metadata' filesep '*channel_metadata*']);
    M = readtable([metadataListing(1).folder filesep metadataListing(1).name]);
    channelNames = [M.Properties.VariableNames;table2cell(M)];
%% analysis    
    for iFolder = folderNumStart:folderNumEnd
        
        listing = dir([mainPath filesep finalSampleFolderList{iFolder} filesep '*.tif']);
        for iFile = 1: numel(listing)
            tic
            subpaths =paths;
            subp=p;
            subpaths.samplefolder = [mainPath filesep finalSampleFolderList{iFolder} filesep ]; 
            subpaths.registration = [subpaths.samplefolder paths.registration ];
            subpaths.metadata = [mainPath filesep paths.metadata filesep];
            subpaths.dearray = [subpaths.samplefolder paths.dearray filesep];
            subpaths.probabilitymaps= [subpaths.samplefolder paths.probabilitymaps filesep];
            subpaths.segmentation = [fileparts(mainPath) filesep 'Train_segmentation' filesep finalSampleFolderList{iFolder} filesep];
            
            subpaths.analysis = [fileparts(mainPath) filesep 'Train_analysis' filesep finalSampleFolderList{iFolder} filesep];
            subp.channelNames= channelNames;
            subpaths.fileExt = FileExt;
            subp.paths = subpaths;
            disp (['Processing ' listing(iFile).name ' in folder ' subpaths.samplefolder])
            kagglesegmenterWrapper(subpaths.samplefolder,listing(iFile).name,subp);
            toc
        end
    end