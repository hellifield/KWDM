function [ fusionImage ] = fuseImages( irImage, visImage, tform )
    rOriginal = imref2d(size(irImage));
    recovered = imwarp(visImage,tform,'OutputView',rOriginal);
    fusionImage = imfuse(irImage, recovered);
end

