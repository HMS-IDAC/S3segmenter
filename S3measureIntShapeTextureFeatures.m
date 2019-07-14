function  [cellID, meanIntRegionTable, medianIntRegionTable,meanStdTable,meanEntTable, haralickTable,majorAxisTable, minorAxisTable, eccTable,solidityTable,areaTable, meanLawTable,centroidCellTable,...
    headMeanInt, headMedianInt,headStd, headEnt,headHaralick,headLaw,headShape] =  S3measureIntShapeTextureFeatures(filePath,region,channelNames,chanRange,varargin)
%If region is 3D, assume to be different compartments in 2D
% assume the 1st slice of the variable region is nuclei, 2nd slice is
% cytoplasm

ip = inputParser;
ip.addParamValue('MedianIntensity','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('Haralick','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('laws','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('crop',[],@(x)(numel(x) == 4 & all(x > 0 )));  
ip.parse(varargin{:});          
p = ip.Results; 

metadata =bfGetReader(filePath);
% if isempty(p.crop)
%     p.crop =[1 1 metadata.getSizeX metadata.getSizeY];
% end
numChan= metadata.getImageCount;


numRegions = size(region,3);
haralick =0;
medianIntensity=0;

shapeTable = [];                        
meanIntRegionTable = [];
medianAllCellsAllRegionsAllChan = [];
meanStdTable = [];
meanEntTable = [];
centroidCellTable = [];
haralickAllCellsAllRegionsAllChan = [];
meanLawTable=[];


regionTypes = {'nuc', 'cyto','plasmem','nucring','cytosol'};
harTypes = {'ang2ndmom' 'contr' 'corr' 'var' 'invdiffmom' 'sumave' 'sumvar' 'sument'...
    'ent' 'diffvar' 'diffent' 'infmeascor1' 'infmeascor2' 'maxcorcoeff'};    
lawTypes = {'LE', 'LR','ES','SS','RR','LS','EE','ER','SR'};
lawWindows = [5 9];

headShape = [];
headMeanInt = [];
headStd =[];
headEnt = [];
headMedianInt =[];
headHaralick = [];
headLaw =[];
   
%% shape characteristics
                        shapeStats{1}=regionprops(region(:,:,1),'Centroid','Area','MajorAxisLength','MinorAxisLength','Eccentricity','Solidity');
                        shapeStats{2}=regionprops(region(:,:,2),'Centroid','Area','MajorAxisLength','MinorAxisLength','Eccentricity','Solidity');
                        
                        for iChan = chanRange(1):chanRange(2)
                            medianAllCellsAllRegions = [];
                            haralickAllCellsAllRegions = [];
                            if isempty(p.crop)
                                image = imresize(imread(filePath,iChan),[size(region,1) size(region,2)]);
                            else
                                image = imresize(imread(filePath,iChan,'PixelRegion',{[p.crop(2),p.crop(2)+p.crop(4)],[p.crop(1),p.crop(1)+p.crop(3)]}),[size(region,1) size(region,2)]);
                            end
                            for iRegion = 1:size(region,3)
                                regionMask = region(:,:,iRegion);
%% mean intensity and std and entropy for each sub region
                                regionStats=regionprops(regionMask,image,'MeanIntensity','PixelValues');
                                meanIntRegionTable= [meanIntRegionTable padarray(cat(1,regionStats.MeanIntensity),double(max(region(:))-max(regionMask(:))),'post')];
                                headMeanInt = cat(2,headMeanInt,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                    char(channelNames(iChan+1,4)) '_int_mean_' regionTypes{iRegion}]});
                                    
                                stdImage = stdfilt(image,ones(3,3));
                                stdStats=regionprops(regionMask,stdImage,'MeanIntensity');
                                meanStdTable = [meanStdTable padarray(cat(1,stdStats.MeanIntensity),double(max(region(:))-max(regionMask(:))),'post')];
                                headStd = cat(2,headStd,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                    char(channelNames(iChan+1,4)) '_txt_standev_' regionTypes{iRegion}]});

                                entImage = entropyfilt(image,ones(3,3));
                                entStats=regionprops(regionMask,entImage,'MeanIntensity');
                                meanEntTable = [meanEntTable padarray(cat(1,entStats.MeanIntensity),double(max(region(:))-max(regionMask(:))),'post')];
                                headEnt = cat(2,headEnt,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                    char(channelNames(iChan+1,4)) '_txt_ent_' regionTypes{iRegion}]});
                                
                                if isequal(p.laws,'true')
                                    lawFeatures = calculateLawFeatures(image,lawWindows);
                                    for iWindow = 1: size(lawFeatures,4)
                                        for iLawFeatures = 1:size(lawFeatures,3)
                                            lawStats=regionprops(regionMask,lawFeatures(:,:,iLawFeatures,iWindow),'MeanIntensity');
                                            meanLawTable = [meanLawTable padarray(cat(1,lawStats.MeanIntensity),double(max(region(:))-max(regionMask(:))),'post')];
                                            headLaw = cat(2,headLaw,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                                char(channelNames(iChan+1,4)) '_txt_laws' lawTypes{iLawFeatures} num2str(lawWindows(iWindow)) '_' regionTypes{iRegion}]});
                                        end
                                    end
                                else
                                    headLawTable =[];
                                    meanLawTable=[];
                                end

                                                                
                                haralickAllCells = [];
                                medianAllCells = [];
                                 
                                for iCell = 1: numel(shapeStats{iRegion})

    %% haralick features for each sub region
                                    if isequal(p.Haralick,'true')
                                        cell = image;
                                        cell(~(regionMask==iCell))=NaN;
                                        haralickAllCells = [haralickAllCells ; computeHaralick(cell)'];
                                    else
                                        haralickAllCells = [haralickAllCells ; zeros(1,14)];
                                    end

    %% median intensity for each sub region
                                    if isequal(p.MedianIntensity,'true')
                                        medianAllCells(end+1) = median(regionStats(iCell).PixelValues);
                                    else
                                        medianAllCells(end+1) = 0;
                                    end
                                end
                                if ~isempty(medianAllCellsAllRegions)
                                   medianAllCells = padarray(medianAllCells,[0 double(size(medianAllCellsAllRegions,1)-numel(medianAllCells))],'post');
                                end
                                medianAllCellsAllRegions = [medianAllCellsAllRegions medianAllCells'];
                                headMedianInt = cat(2,headMedianInt,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                    char(channelNames(iChan+1,4)) '_int_med_' regionTypes{iRegion}]});

                                if ~isempty( haralickAllCellsAllRegions)
                                   haralickAllCells = padarray(haralickAllCells,double(size(haralickAllCellsAllRegions,1)-size(haralickAllCells,1)),'post');
                                end
                                haralickAllCellsAllRegions = [haralickAllCellsAllRegions  haralickAllCells];
                                 for i = 1:14
                                    headHaralick = cat(2,headHaralick,{[char(channelNames(iChan+1,3)) '_' num2str(channelNames{iChan+1,1}) '_' ...
                                        char(channelNames(iChan+1,4)) '_txt_Har' harTypes{i} '_' regionTypes{iRegion}]});
                                 end
                            
                            end
                            medianAllCellsAllRegionsAllChan = [medianAllCellsAllRegionsAllChan medianAllCellsAllRegions];
                            haralickAllCellsAllRegionsAllChan = [haralickAllCellsAllRegionsAllChan haralickAllCellsAllRegions];
                            disp(['Measured features for channel ' num2str(iChan)])
                        end
                        cellID = [1:numel(shapeStats{1})]';
                        majorAxisTable = [cat(1,shapeStats{1}.MajorAxisLength) padarray(cat(1,shapeStats{2}.MajorAxisLength),double(numel(shapeStats{1})-numel(shapeStats{2})),'post') ];
                        minorAxisTable = [cat(1,shapeStats{1}.MinorAxisLength) padarray(cat(1,shapeStats{2}.MinorAxisLength),double(numel(shapeStats{1})-numel(shapeStats{2})),'post') ];
                        eccTable = [cat(1,shapeStats{1}.Eccentricity) padarray(cat(1,shapeStats{2}.Eccentricity),double(numel(shapeStats{1})-numel(shapeStats{2})),'post') ];
                        solidityTable = [cat(1,shapeStats{1}.Solidity) padarray(cat(1,shapeStats{2}.Solidity),double(numel(shapeStats{1})-numel(shapeStats{2})),'post') ];
                        areaTable = [cat(1,shapeStats{1}.Area) padarray(cat(1,shapeStats{2}.Area),double(numel(shapeStats{1})-numel(shapeStats{2})),'post') ];
                        centroidCellTable = cat(1,shapeStats{1}.Centroid);       
                        medianIntRegionTable = medianAllCellsAllRegionsAllChan;
                        haralickTable = haralickAllCellsAllRegionsAllChan;
                        headShape = {'mor_area_nuc' 'mor_area_cyto' 'mor_cellposx_nuc' 'mor_cellposy_nuc' 'mor_majoraxis_nuc' 'mor_majoraxis_cyto' ...
                            'mor_minoraxis_nuc' 'mor_minoraxis_cyto' 'mor_Eccen_nuc' 'mor_Eccen_cyto' 'mor_solidity_nuc' 'mor_solildity_cyto'};
                        
   
    function x= computeHaralick(cell)                        
    glcm = graycomatrix(cell, 'offset', [0 1], 'GrayLimits',[],'Symmetric', true);
    x = haralickTextureFeatures(glcm,[2 3 4 7 9]);
    end

    function lawFeatures = calculateLawFeatures(I,lawWindows)
        lawFeatures = [];
        for iwindow = 1:numel(lawWindows)
            lawFeatures = cat(4,lawFeatures,reshape(cell2mat(laws(I,lawWindows(iwindow))),[size(I,1) size(I,2) 9]));
        end
    end

end