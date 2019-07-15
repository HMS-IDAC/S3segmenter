function S3MeasureFeatures(regions,paths,fileName,varargin)   
ip = inputParser;
ip.addParamValue('MedianIntensity','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('crop',[],@(x)(numel(x) == 4 & all(x > 0 )));  
ip.addParamValue('chanRange',[0],@(x)(numel(x) >0 & all(x > -1 )));
ip.addParamValue('channelNames',[]);
ip.parse(varargin{:});          
p = ip.Results; 

if isempty(p.channelNames)
    listing = dir([paths.metadata filesep '*channel_metadata*']);
    if isempty(listing)
        disp('No metadata file found!')
        return
    else
    %     [~,~,channelNames] = xlsread([listing(1).folder filesep listing(1).name]);
          M = readtable([listing(1).folder filesep listing(1).name]);
          channelNames = [M.Properties.VariableNames;table2cell(M)];
    end
else 
    channelNames = p.channelNames;
end

metadata =bfGetReader([paths.registration fileName]);
numChan= metadata.getImageCount;

if p.chanRange==0
    p.chanRange(1)=1;
    p.chanRange(2)=numChan;
end


mkdir(paths.analysis)
name=extractBefore(fileName,'.');

[cellID, meanIntRegionTable, medianIntRegionTable,meanStdTable,meanEntTable, haralickTable,majorAxisTable, minorAxisTable, eccTable,solidityTable,areaTable, meanLawTable,centroidCellTable,...
                       headMeanInt, headMedianInt,headStd, headEnt,headHaralick,headLaw,headShape] =...
                       S3measureIntShapeTextureFeatures([paths.registration fileName] ,regions,channelNames,[p.chanRange(1) p.chanRange(2)],...
                       'MedianIntensity',p.MedianIntensity,'crop',p.crop);
                   
if ~isempty(meanIntRegionTable)
     writetable(array2table([cellID, meanIntRegionTable, medianIntRegionTable,areaTable,centroidCellTable, majorAxisTable, minorAxisTable, eccTable, ...
   solidityTable, meanStdTable,meanEntTable,haralickTable,meanLawTable],...
   'VariableNames',lower(regexprep(['cellid',headMeanInt,headMedianInt,headShape,headStd, headEnt,headHaralick,headLaw], '[()/\-\s+]', ''))),...
                [paths.analysis filesep name '_Features.txt'],'Delimiter','\t')
end
     