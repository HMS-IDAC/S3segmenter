function [imMultiscaleLoGResponse,centers] = filterMultiScaleMultiDirDConstrLoG(image,mask,seedMethod,downSample,raw)


switch seedMethod
    case 'LoG'
        origimageSize = size(image);
        sigmax = [ 2.5 4 6] * downSample;
        sigmay = [ 2.5 4 6] * downSample;

        % testSizes=1:2:9;
        % for isigma = 1:numel(sigmax)
        %     maxTest=[];
        %     for iSize= testSizes
        %         testI = zeros(31,31);
        %         testI(16-(iSize-1)/2:16+(iSize-1)/2,16-(iSize-1)/2:16+(iSize-1)/2)=64;
        %         gLoG= multiScaleAngleLoGKernel(testI,sigmax(isigma), sigmay(isigma),1);
        %         gLoG= fft2(gLoG);
        %         testLog = real(ifft2(fft2(testI).*gLoG));
        %         
        %         maxTest = cat(2, maxTest,testLog(16,16));
        %     end
        %     [nucleiLogResponse(isigma),nucleiLogThInd]=max(maxTest(:));
        %     nucleiRadius(isigma) = testSizes(nucleiLogThInd);
        % 
        % end


        mask = (max(image(:))-image)>90;
        dsimage = imresize(image,downSample);
        % % 
        nang = 4; 
        DT = bwdist(~mask);
        DT=2.*log2(DT);
        DT(DT==-Inf)=0;
        LoGResponse=[];
        clear mask, clear image
        % % LoG 
        I = fft2(max(dsimage(:))-dsimage); 
        % seedMask = ones(size(I));
        % seeds = 0.*seedMask;

        for iang = 1:nang
            for isigma = 1:numel(sigmay)
                gLoG= multiScaleAngleLoGKernel(dsimage,sigmax(isigma), sigmay(isigma),iang);
                gLoG= fft2(gLoG);
                LoGResponse = cat(3,LoGResponse,real(ifft2(I.*gLoG)));
        %         newSeeds =  seedMask.*((maxLoGResponse.*LoGResponse)>0.5*nucleiLogResponse(isigma));
        %         newSeedMask = (bwdist(newSeeds)>nucleiRadius(isigma));
        %         seeds = seeds + newSeeds;
        %         seedMask = seedMask.*newSeedMask;

        %         maskResponse = cat(3, maskResponse, real(ifft2(I.*gLoG)));%.*(DT>max([sigmax(isigma) sigmay(isigma)])));
            end
        end

    case 'globalThreshold'
        imMultiscaleLoGResponse=imregionalmax(imhmax(imgaussfilt(image,downSample/30),downSample/30));
%         imMultiscaleLoGResponse=imregionalmax(imgaussfilt(uint8(imextendedmax(image,10)),1));
        mask = imfill(mask,'holes');
        centers = mask.*imfill(bwmorph(imMultiscaleLoGResponse,'shrink','Inf'),'holes');

%         invertImage =max(image(:))-image;
%         Ilog = imresize(filterLoG(imresize(invertImage,0.25),1),size(mask));
%         Imax = imregionalmax(Ilog).*(~mask);
%         imMultiscaleLoGResponse=imMultiscaleLoGResponse+Imax;

    case 'gaussian'

% imMultiscaleLoGResponse = imresize(max(imregionalmax(LoGResponse),[],3),origimageSize);
% imshowpair(imresize(dsimage,origimageSize),imdilate(imMultiscaleLoGResponse,strel('square',2)))
%  disp ('done')

% 
% % gaussian
% image=max(image(:))-image;
% sigmax = [ 1 4];
% sigmay = [ 1 4];
% maskResponse=[];
% for iang = 1:nang
%     for isigma = 1:numel(sigmay)
%         ang = pi/nang*(iang-1);
%         gkernel = anisotropicGaussD(sigmax(isigma),sigmay(isigma),ang);
%         imCurLoGResponse = conv2(max(dsimage(:))-dsimage,gkernel,'same');
%         maskResponse = cat(3, maskResponse,imCurLoGResponse.*(DT>max([sigmax(isigma) sigmay(isigma)])));
% %                 maskResponse = cat(3, maskResponse,imCurLoGResponse);
%     end
% end
% imMultiscaleLoGResponse = imresize(imregionalmax(imgaussfilt3(max(maskResponse,[],3),1.5)),size(dsimage));
%  figure,imshowpair(imresize(dsimage,1),imdilate(imMultiscaleLoGResponse,strel('square',4)))
%    disp ('done')
end
