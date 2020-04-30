function [cytoplasmMask,nucleiMask,finalCellMask] = S3CytoplasmSegmentation(nucleiMask,cyto,varargin)

ip = inputParser;
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'hybrid','ilastik','distanceTransform','bwdistanceTransform','ring','UNet'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resize',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('cytoDilation',5,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('sizeFilter',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParamValue('cytoPM', [], @(x) isnumeric(x));
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
                         
            case 'hybrid'
                cytoBlur = imgaussfilt(cyto,2);
%                 grad = imgradient(cytoBlur);
                grad=stdfilt(cytoBlur,true(3,3));  
                Idist= bwdist(nucleiMask>0);
                Idist= sqrt(grad.^2 + 0.00001*max(grad(:))/max(Idist(:))*Idist.^2);
                cytograd = imimposemin(Idist,nucleiMask>0 | (imerode(1-mask,strel('disk',3)))  );  
                cellMask = watershed(cytograd);
                mask=ones(size(cellMask));
%                   cellMask = cellMask.*cast(mask,class(cellMask));
%                   stats=regionprops(cellMask,nucleiMask>0,'MaxIntensity','Area');
%                   idx = find([stats.MaxIntensity]>0 );
%                   finalCellMask = ismember(cellMask,idx);
%                     imshowpair(finalCellMask,imadjust(cyto))
%                 % get cells with cytoplasm
%                 cytoMarkerMask = bwareaopen(imfill(mask,'holes'),200);
%                 
%                 % get cells without cytoplasm
%                 stats=regionprops(nucleiMask,cytoMarkerMask,'MaxIntensity');
%                 idx = find([stats.MaxIntensity]==0 );
%                 nocytoNuclei = ismember(nucleiMask,idx);
%                 Itest  = bwmorph(nocytoNuclei,'shrink','Inf');
%                 nocytoNucleicyto = bwdist(Itest | cytoMarkerMask);
%                 nocytoNucleicyto(cytoMarkerMask)=0;
%                 nocytoNucleicyto(nocytoNucleicyto>10)=0;
%                 
%                 
%                 cytoNuclei = (nucleiMask>0) - (nocytoNuclei>0);
%                 Imin  = bwmorph(cytoNuclei,'shrink','Inf');
%                 gdist = graydist(cyto+ nocytoNucleicyto,Imin|Itest).*cytoMarkerMask;
%                 cytograd = imimposemin(gdist, 1-(cytoMarkerMask | nocytoNucleicyto>0)|nucleiMask>0);
%                 cellMask=watershed(cytograd);
%                 mask = ones(size(cellMask,1),size(cellMask,2));
% %                 stats = regionprops(cellMask,nucleiMask,'MaxIntensity','Area');
% %                 idx = find([stats.MaxIntensity]>0 & [stats.Area]<5000);
% %                 finalCellMask = ismember(cellMask,idx);
                
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
                cellMask = bwlabel(bwmorph(nucleiMask>0,'thicken',p.cytoDilation));
                mask = ones(size(cellMask));
                
            case 'ilastik'
                gdist = bwdist(nucleiMask>0);
                cytograd = imimposemin(gdist, nucleiMask>0 );
                cellMask=watershed(cytograd);
                cellMask = cellMask.*cast(mask,class(cellMask));
            case 'UNet'
                % split cytoplasm
                cellMembrane = p.cytoPM < thresholdOtsu(imresize(p.cytoPM,0.25));
                Idist = -bwdist(~cellMembrane);
                Imin = imregionalmin(imhmin(Idist,10));
                cytograd = imimposemin(Idist, Imin );
                cellMask=watershed(cytograd);
                
                % remove all nuclei that touch membranes
                stats = regionprops(nucleiMask,cellMask==0,'MaxIntensity');
                idx = find([stats.MaxIntensity]<1 );
                nucleiMask = ismember(bwlabel(nucleiMask),idx);
                
                mask = cellMask>0;
                
                % find cytoplasm with leftover nuclei
%                 stats = regionprops(cellMask,nucleiMaskTest>0,'MaxIntensity');
%                 idx = find([stats.MaxIntensity]>0);
%                 finalCellMask = bwlabel(ismember(cellMask,idx));
%                 nucleiMask = finalCellMask.*(nucleiMaskTest>0);
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
          