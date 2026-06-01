function [fx, fy] = makeFrequencyGrid(imageHeight, imageWidth, pixelSizeNm)
%MAKEFREQUENCYGRID Return centered frequency coordinates in cycles/nm.

fxAxis = ((1:imageWidth) - floor(imageWidth/2) - 1) ./ (imageWidth * pixelSizeNm);
fyAxis = ((1:imageHeight) - floor(imageHeight/2) - 1) ./ (imageHeight * pixelSizeNm);
[fx, fy] = meshgrid(fxAxis, fyAxis);
end
