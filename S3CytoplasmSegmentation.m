function [cytoplasmMask,nucleiMask,finalCellMask] = S3CytoplasmSegmentation(nucleiMask,cyto,modelCat,varargin)

ip = inputParser;
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','bwdistanceTransform','ring'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resize',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('sizeFilter',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('mask', [], @(x) isnumeric(x) || islogical(x));
ip.parse(varargin{:});          
p = ip.Results;  
        
%% mask
if ~isempty(p.mask)
    mask = (nucleiMask + imresize(p.mask,size(nucleiMask)))>0;
else
    mask = ones(size(nucleiMask));
end

%% cytoplasm segmentation methods
        
        switch p.cytoMethod
            case 'RF'
          F = pcImageFeatures(imresize(double(cyto)/65335,p.resize,'bilinear'),modelCat.sigmas,modelCat.offsets,...
              modelCat.osSigma,modelCat.radii,modelCat.cfSigma,modelCat.logSigmas,modelCat.sfSigmas,...
              modelCat.ridgeSigmas,modelCat.ridgenangs,modelCat.edgeSigmas,modelCat.edgenangs,modelCat.nhoodEntropy,...
              modelCat.nhoodStd);
             [imL,catClassProbs] = imClassify(F,modelCat.treeBag,100);
              contours = imresize(catClassProbs(:,:,2),2);
              bgm =imresize(bwmorph( imgaussfilt3(cyto,2)<100,'thin',Inf),2);
              cytograd= imimposemin(imresize(contours,[size(nucleiMask,1) size(nucleiMask,2)]),bgm|nucleiMask);
              cellMask= watershed(cytograd);
                
            case 'distanceTransform'
                gdist = graydist(cyto,bwmorph(nucleiMask>0,'shrink','Inf'));
                cytograd = imimposemin(gdist, (imerode(1-mask,strel('disk',5))) | nucleiMask>0 );
                cellMask=watershed(cytograd);
                clear cytograd, clear gdist
                cellMask = cellMask.*cast(mask,class(cellMask));
                
            case 'bwdistanceTransform'
                gdist = bwdist(nucleiMask>0);
                cytograd = imimposemin(gdist, nucleiMask>0 );
                cellMask=watershed(cytograd);
                cellMask = cellMask.*cast(mask,class(cellMask));
            case 'contours'
                
                contours = normalize(steerableDetector(im2double(cyto),2,1.5));

            case 'ring'
                cellMask = bwlabel(bwmorph(nucleiMask>0,'thicken',9));
                mask = ones(size(cellMask));
        end

         

%% filter based on tissue 
            stats=regionprops(cellMask,mask,'MeanIntensity','Area');
            idx = find([stats.MeanIntensity] > 0.05 );
            tissueCellMask = bwlabel(mask.*ismember(cellMask,idx));
            
%% filter based on a nuclei
          stats=regionprops(tissueCellMask,nucleiMask>0,'MaxIntensity','Area');
          clear cyto, clear mask
          idx = find([stats.MaxIntensity] > 0 & [stats.Area]>5  & [stats.Area]<prctile(cat(1,stats.Area),99.9));
          finalCellMask = bwlabel(ismember(tissueCellMask,idx));
          clear tissueCellMask
          nucleiMask = cast(nucleiMask>0,class(finalCellMask)).*finalCellMask; 
          cytoplasmMask = finalCellMask - nucleiMask;
          