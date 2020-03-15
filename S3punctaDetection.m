function spots = S3punctaDetection(image,sigma,alpha)
    Ilog = filterLoG(image,sigma);
    Imax = imregionalmax(Ilog);
    threshold = median(Ilog(Imax))+ alpha*1.4826*mad(Ilog(Imax),1);
    spots = Imax.*(Ilog>threshold);
end