function O2batchS3segmenterWrapperR(mainPath,varargin)


ip = inputParser;
ip.addParamValue('HPC','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('fileNum',1,@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('ClassProbSource','unet',@(x)(ismember(x,{'RF','ilastik','unet','none'})));
ip.addParamValue('NucMaskChan',[1],@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('CytoMaskChan',[2],@(x)(numel(x) > 0 & all(x > -1 )));  
ip.addParamValue('TissueMaskChan',[],@(x)(isnumeric(x))); 
ip.addParamValue('RefineTissueMask',[0],@(x)(numel(x) > 0 & all(x > 0 ))); 
ip.addParamValue('cytoDilation',5,@(x)(numel(x) > 0 & all(x > 0 ))); 
ip.addParamValue('mask','tissue',@(x)(ismember(x,{'TMA','tissue','none'}))); % set to true if sample is TMA cores
ip.addParamValue('crop','noCrop',@(x)(ismember(x,{'interactiveCrop','autoCrop','dearray','noCrop'})));
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'hybrid','ilastik','distanceTransform','bwdistanceTransform','ring','UNet'})));
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
ip.addParamValue('detectPuncta',0,@(x)(numel(x) >0 & all(x > -1 )));
ip.addParamValue('punctaSigma',1,@(x)(numel(x) >0 & all(x > 0 )));
ip.addParamValue('punctaSD',4,@(x)(numel(x) >0 & all(x > 0 )));
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
if isequal(p.ClassProbSource,'ilastik')
    paths.probabilitymaps= ['prob_maps_ilastik'];
    p.probMapOrder = [1 2 3 2];
else
    paths.probabilitymaps= ['prob_maps'];
    p.probMapOrder = [1 2 3 4];
end
paths.segmentation = ['segmentation'];
paths.analysis = ['analysis'];
paths.registration = ['registration'];

if isequal(p.crop,'dearray')
    FileExt = 'tif';
    searchPath = ['dearray'];
else
    FileExt = 'tif';
    searchPath = ['registration'];
end

%% read file names
    sampleFolderList = dir([mainPath filesep '**' filesep searchPath]);
    finalSampleFolderList = [];
    finalSampleFileList =[];
    for iFolder = 1:length(sampleFolderList)
        fName = [sampleFolderList(iFolder).folder filesep sampleFolderList(iFolder).name];
        if isfolder(sampleFolderList(iFolder).folder) && ~isequal(sampleFolderList(iFolder).name,'..') ...
            && ~isequal(sampleFolderList(iFolder).name,'.') && ~contains(sampleFolderList(iFolder).name,'TMA_MAP') && contains(sampleFolderList(iFolder).name,FileExt)
            finalSampleFolderList{end+1} = fileparts(sampleFolderList(iFolder).folder);
            finalSampleFileList{end+1} = sampleFolderList(iFolder).name;
        end
    end
    disp (['Found ' num2str(length(finalSampleFolderList)) ' samples(s)!'])
    
    if isequal(p.HPC,'false')
        fileNumStart = 1;
        fileNumEnd = length(finalSampleFolderList);
    else 
        fileNumStart=p.fileNum;
        fileNumEnd = p.fileNum;
    end

    for iFile = fileNumStart:fileNumEnd
        tic
        subpaths =paths;
        subp=p;
        subpaths.samplefolder = [finalSampleFolderList{iFile} filesep ]; 
        subpaths.registration = [subpaths.samplefolder paths.registration filesep];
        subpaths.metadata = [subpaths.samplefolder paths.metadata filesep];
        subpaths.dearray = [subpaths.samplefolder paths.dearray filesep];
        subpaths.probabilitymaps= [subpaths.samplefolder paths.probabilitymaps filesep];
        subpaths.segmentation = [subpaths.samplefolder paths.segmentation filesep];
        subpaths.analysis = [subpaths.samplefolder paths.analysis filesep];
        subpaths.fileExt = FileExt;
        subp.paths = subpaths;
        disp (['Processing ' finalSampleFileList{iFile}])
        S3segmenterWrapper(subpaths.samplefolder,finalSampleFileList{iFile},subp);
        toc
    end