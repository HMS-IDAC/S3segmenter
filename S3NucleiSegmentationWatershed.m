function [nucleiMask,largestNucleusArea] = S3NucleiSegmentationWatershed(NucleiPM,nucleiImage,logSigma,varargin)

ip = inputParser;
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','watershedBWDist','dilation'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('inferNucCenters','UNet',@(x)(ismember(x,{'UNet','RF','Int'})));
ip.addParamValue('resize',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('PMUpSampleFactor',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParamValue('nucleiFilter','Int',@(x)(ismember(x,{'LoG','Int','IntPM','none'})));
ip.parse(varargin{:});          
p = ip.Results;  


%% preprocess
if size(NucleiPM,3)>2
    nucleiCentersResized = imresize(NucleiPM(:,:,1),p.resize);
    nucleiContoursResized = imresize(NucleiPM(:,:,2),p.resize);
else
    nucleiContoursResized = imresize(NucleiPM(:,:,1),p.resize);
    nucleiCentersResized = max(nucleiContoursResized(:))-nucleiContoursResized;
end

if ~isequal(p.resize,1)
    nucleiContoursResized = imresize(NucleiPM(:,:,1),p.resize);
    nucleiImageResized = imresize(nucleiImage,size(nucleiContoursResized));
    
else
    nucleiImageResized = nucleiImage;
   
end
clear NucleiPM, clear nucleiImage

%% mask
if isempty(p.mask) || sum(sum(p.mask))==0
    %% sample background for nuclei intensity filter for TMA cores
%     mask = imresize(p.mask,size(nucleiImageResized));
%     MITh=multithresh(nucleiImageResized,2);
%     MITh=MITh(1);
% else
%     mask = ones(size(nucleiImageResized));
%     MITh=multithresh(nucleiImageResized,2);
%     MITh=MITh(1);
    p.mask = ones(size(nucleiImageResized));
else
    p.mask = imresize(p.mask,size(nucleiImageResized));
end


%% use nuclei contour class prob maps from UNet or RF
%     if isequal(p.inferNucCenters,'RF') || ~isempty(nucleiCentersResized) %RF
%         nucleiClass = nucleiCentersResized;
%     elseif isequal(p.inferNucCenters,'UNet') % UNet
%         nucleiClass =  nucleiContoursResized;
%     else
%         nucleiClass = nucleiImageResized;
%     end
    
   %% markers based on log filter on classProbs 3
     if isequal(p.useGPUArray,'true')
      logfgm=  imregionalmax(gather(imgaussfilt3(gpuArray(filterLoG(nucleiContoursResized,logSigma)),2))); %3.5, 2
     else
        if numel(logSigma)==1
         nucleiDiameter  = [logSigma*0.5 logSigma*1.5];
        else
         nucleiDiameter = logSigma;
        end
%         numLoGScales = ceil((nucleiDiameter(2)-nucleiDiameter(1))/3);
%         logmask = (nucleiImageResized>thresholdMinimumError(nucleiImageResized,'model','poisson'));
        logmask = nucleiCentersResized>0.6*max(nucleiCentersResized(:));
        [logfgm,centers] = filterMultiScaleMultiDirDConstrLoG(max(nucleiContoursResized(:)) - nucleiContoursResized,logmask,'globalThreshold',nucleiDiameter(2),nucleiImageResized);
     end
%      figure,imshowpair(logfgm,nucleiCentersResized)
%         logfgm=logfgm.*p.mask;
%         masktest=(nucleiCentersResized+nucleiContoursResized)>thresholdOtsu(nucleiCentersResized+nucleiContoursResized);
%         bg=(1-masktest).*p.mask;
%         bg = imerode(bg,strel('disk',3));
%% apply watershed transform
switch p.nucleiRegion 
    case 'watershedContourInt'
      gradmag2= imimposemin(nucleiContoursResized,imresize(logfgm,size(nucleiContoursResized),'nearest'));
      foregroundMask= watershed(gradmag2);
      clear gradmag2
    case 'dilation'
      foregroundMask = imdilate(logfgm,strel('square',2));

    case 'watershedContourDist'
        
        gdist = graydist(nucleiContoursResized,imresize(logfgm,size(nucleiContoursResized),'nearest'));
        cytograd= imimposemin(gdist,imresize(logfgm,size(gdist)));
        foregroundMask=watershed(cytograd);
    case 'watershedBWDist'
%     fgth=imtophat(nucleiCentersResized,strel('disk',7));
    fgMask = nucleiCentersResized>thresholdOtsu(nucleiCentersResized);%imclose(fgth>50,strel('disk',3));
    IDist = -bwdist(~fgMask);
    logfgm = imregionalmin(imhmin(IDist,1));
    gdist = graydist(nucleiContoursResized,logfgm);
    
    cytograd= imimposemin(IDist,logfgm);
    foregroundMask=watershed(cytograd);
    foreground =nucleiCentersResized>thresholdMinimumError(nucleiCentersResized,'model','poisson');
    foregroundMask = foregroundMask.*cast(foreground,class(foregroundMask));
%     logfgm = (logfgm.*foregroundMask)>0; 
 end

    %% process mask
   allNuclei = bwlabel((foregroundMask>0).*p.mask);
   
   if isequal(p.inferNucCenters,'false')
       nucleiCentersResized = imresize(nucleiCentersResized,size(nucleiContoursResized));
    stats=regionprops(allNuclei,nucleiCentersResized+nucleiContoursResized,'MeanIntensity','MinIntensity','Area');   
    idx = find([stats.MeanIntensity] > 0.75 & [stats.Area] > 20 & [stats.Area] < median(cat(1,stats.Area))*5 );
   else
%     stats=regionprops(allNuclei,mask,'MeanIntensity');
%     idx = find([stats.MeanIntensity] > 0.25);
%     nucleiTissueMask = bwlabel(ismember(allNuclei,idx));   
      
    if isequal(p.nucleiFilter,'LoG')
        statsFilt=regionprops(allNuclei,(double(nucleiCentersResized))/255.*double(nucleiImageResized),'MaxIntensity');
        MITh = median(cat(1,statsFilt.MaxIntensity))+0.5*mad(cat(1,statsFilt.MaxIntensity));
    elseif isequal(p.nucleiFilter,'Int')    
        statsFilt=regionprops(logfgm.*allNuclei,nucleiImageResized,'MeanIntensity');
        MITh = thresholdMinimumError(nucleiImageResized,'model','poisson');
    elseif isequal(p.nucleiFilter,'IntPM')
        statsFilt=regionprops(allNuclei,nucleiCentersResized,'MeanIntensity');
        MITh = thresholdMinimumError(cat(1,statsFilt.MeanIntensity),'model','poisson');
    else
        statsFilt=regionprops(allNuclei,nucleiImageResized,'MeanIntensity');
        MITh =0;
    end

    idx = find([statsFilt.MeanIntensity]>MITh );
    allNuclei = bwlabel(ismember(allNuclei,idx));
    stats=regionprops(allNuclei,'Area','Solidity');
    idx = find( [stats.Area]>(nucleiDiameter(1)^2)*3/4 & [stats.Area] < (nucleiDiameter(2)^2)*3/4  & [stats.Solidity]>0.8);
   end

   nucleiMask = ismember(allNuclei,idx);
   statsNM=regionprops(imresize(nucleiMask,[size(allNuclei,1) size(allNuclei,2)]),'Area');
   largestNucleusArea=prctile([statsNM.Area],95);


