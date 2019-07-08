function [F,featNames] = pcImageFeatures(stack,sigmas,offsets,osSigma,radii,cfSigma,logSigmas,sfSigmas,ridgeSigmas,nRidgeAngs,edgeSigmas,nEdgeAngs,entropyNhoods,stdNHoods)

% tic

F = [];
featIndex = 0;
featNames = {};
for iChan = 1:size(stack,3)
I=stack(:,:,iChan);

    derivNames = {'d0','dx','dy','dxx','dxy','dyy','hessEV1','hessEV2'};
    for sigma = sigmas
        for i = 1:length(derivNames)
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('sigma%d',sigma) derivNames{i} ' Chan' int2str(iChan)];
        end
        D = zeros(size(I,1),size(I,2),8);
        [D(:,:,1),D(:,:,2),D(:,:,3),D(:,:,4),D(:,:,5),D(:,:,6),D(:,:,7),D(:,:,8)] = derivatives(I,sigma);
        F = cat(3,F,D);
        featIndex = featIndex+1;
        featNames{featIndex} = [sprintf('sigma%dedges',sigma) ' Chan' int2str(iChan)];
        F = cat(3,F,sqrt(D(:,:,2).^2+D(:,:,3).^2)); % edges
    end

    if ~isempty(offsets)
        J = filterGauss2D(I,osSigma);
        for r = offsets
            aIndex = 0;
            for a = 0:pi/4:2*pi-pi/4
                aIndex = aIndex+1;
                v = r*[cos(a) sin(a)];
                T = imtranslate(J,v,'OutputView','same');
                F = cat(3,F,T);
                featIndex = featIndex+1;
                featNames{featIndex} = [sprintf('offset%da%d',r,aIndex) ' Chan' int2str(iChan)];
            end
        end
    end

    if ~isempty(radii)
        for r = radii
            [C1,C2] = circlikl(imresize(I,0.5),r,cfSigma,8,0.001);
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('circfeat%dC1',r) ' Chan' int2str(iChan)];
            F = cat(3,F,imresize(C1,2));
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('circfeat%dC2',r) ' Chan' int2str(iChan)];
            F = cat(3,F,imresize(C2,2));
        end
    end

    if ~isempty(logSigmas)
        for sigma = logSigmas
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('sigma%dlog',sigma) ' Chan' int2str(iChan)];
            F = cat(3,F,filterLoG(I,sigma));
        end
    end

    if ~isempty(sfSigmas)
        for sigma = sfSigmas
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('sigma%dsteer',sigma) ' Chan' int2str(iChan)];
            F = cat(3,F,steerableDetector(I,4,sigma));
        end
    end
    
    if ~isempty(ridgeSigmas)
        for sigma = ridgeSigmas
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('sigma%dRidge',sigma) ' Chan' int2str(iChan)];
            F = cat(3,F,ridgelikl(I,nRidgeAngs,sigma));
        end
    end
    
    if ~isempty(edgeSigmas)
        for sigma = edgeSigmas
            featIndex = featIndex+1;
            featNames{featIndex} = [sprintf('sigma%dEdge',sigma) ' Chan' int2str(iChan)];
            F = cat(3,F,edgeliklRF(I,nEdgeAngs,sigma));
        end
    end

    if ~isempty(entropyNhoods)
       for nhood = entropyNhoods
          featIndex = featIndex+1;
          featNames{featIndex} = [sprintf('nHood%dEntropyFilt',nhood) ' Chan' int2str(iChan)];
          F = cat(3,F,entropyfilt(I,true(nhood)));
       end
    end


    if ~isempty(stdNHoods)
       for nhood = stdNHoods
          featIndex = featIndex+1;
          featNames{featIndex} = [sprintf('nHood%dStdFilt',nhood) ' Chan' int2str(iChan)];
          F = cat(3,F,stdfilt(I,true(nhood)));
       end
    end
    % toc
end

end