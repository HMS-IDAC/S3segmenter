function [totalClassProbs,nucleiCrop]=RFSplitter(nuclei,modelNuc,varargin)
% model = pixelClassifierTrain('D:\LSP\cycif\S3\nucleiTrainingSet','sigmas',[1 2 4],'radii',[10 20],'cfSigma',4,'logSigmas',[2 4],'edgeSigmas',[1 3] ,'ridgeSigmas',[3],'nhoodEntropy',[5 11],'nhoodStd',[5 11 17],'pctMaxNPixelsPerLabel',20,'adjustContrast',false);
% modelContour = pixelClassifierTrain('D:\LSP\cycif\S3\cytoTrainingSet','sigmas',[2 4],'radii',[5 20],'cfSigma',2,'sfSigmas',[1 2],'logSigmas',[1 2 4],'edgeSigmas',[1 3] ,'ridgeSigmas',[1 3],'nhoodEntropy',[5 11],'nhoodStd',[],'pctMaxNPixelsPerLabel',50,'adjustContrast',false);
        

ip = inputParser;
ip.addParamValue('split',false,@islogical);
ip.parse(varargin{:});          
p = ip.Results;  

% get a folder with single channels
numClass=numel(modelNuc.treeBag.ClassNames);
totalClassProbs=zeros(size(nuclei,1),size(nuclei,2),2);

if p.split
    for i = 1:8
        for j = 1:8
            nucleiCrop = nuclei((1+size(nuclei,1)/16*(i-1)):(size(nuclei,1)/16*i),(1+size(nuclei,2)/16*(j-1)):(size(nuclei,2)/16*j));

             F = pcImageFeatures(double(nucleiCrop)/65535,modelNuc.sigmas,modelNuc.offsets,modelNuc.osSigma,modelNuc.radii,modelNuc.cfSigma,...
                modelNuc.logSigmas,modelNuc.sfSigmas,modelNuc.ridgeSigmas,modelNuc.ridgenangs,modelNuc.edgeSigmas,modelNuc.edgenangs,...
                modelNuc.nhoodEntropy,modelNuc.nhoodStd);
            [~,classProbs] = imClassify(F,modelNuc.treeBag,100);

    %         totalClassProbs((1+size(nuclei,1)/8*(i-1)):(size(nuclei,1)/8*i),(1+size(nuclei,2)/8*(j-1)):(size(nuclei,2)/8*j),1:2)=classProbs(:,:,2:3);
    totalClassProbs=classProbs(:,:,2:3);
    %         tiffwriteimj(uint8(finalNuclei), [testPath filesep int2str(i) '_' int2str(j) '.tif'])
        end
    end

else
    i=3;
    j=3;
    nucleiCrop = nuclei((1+size(nuclei,1)/16*(i-1)):(size(nuclei,1)/16*i),(1+size(nuclei,2)/16*(j-1)):(size(nuclei,2)/16*j));
    
    F = pcImageFeatures(double(nucleiCrop)/65535,modelNuc.sigmas,modelNuc.offsets,modelNuc.osSigma,modelNuc.radii,modelNuc.cfSigma,...
                modelNuc.logSigmas,modelNuc.sfSigmas,modelNuc.ridgeSigmas,modelNuc.ridgenangs,modelNuc.edgeSigmas,modelNuc.edgenangs,...
                modelNuc.nhoodEntropy,modelNuc.nhoodStd);
            [~,classProbs] = imClassify(F,modelNuc.treeBag,100);

    %         totalClassProbs((1+size(nuclei,1)/8*(i-1)):(size(nuclei,1)/8*i),(1+size(nuclei,2)/8*(j-1)):(size(nuclei,2)/8*j),1:2)=classProbs(:,:,2:3);
    totalClassProbs=classProbs(:,:,2:3);
    
end


