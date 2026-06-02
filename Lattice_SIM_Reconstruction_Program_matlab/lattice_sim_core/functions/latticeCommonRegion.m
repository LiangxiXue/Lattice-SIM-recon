function [newb0, newb1] = latticeCommonRegion(band0, band1, otfValues, carrierPixel, dist, weightLimit, divideByOtf)
%LATTICECOMMONREGION Select HiFi-style common support for Lattice bands.

if nargin < 5 || isempty(dist)
    dist = 0.15;
end
if nargin < 6 || isempty(weightLimit)
    weightLimit = 1e-4;
end
if nargin < 7
    divideByOtf = true;
end

[h, w] = size(band0);
cnt = [h/2 + 1, w/2 + 1];

weight0 = abs(otfValues);
weight1 = abs(otfValues);
wt0 = abs(shiftOtfByCarrier(otfValues, carrierPixel));
wt1 = abs(shiftOtfByCarrier(otfValues, -carrierPixel));

if max(weight0(:)) > 0
    weightThreshold = max(weight0(:)) * weightLimit;
else
    weightThreshold = weightLimit;
end

mask0 = abs(weight0) < weightThreshold | abs(wt0) < weightThreshold;
mask1 = abs(weight1) < weightThreshold | abs(wt1) < weightThreshold;

newb0 = band0;
newb1 = band1;
newb0(mask0) = 0;
newb1(mask1) = 0;

weight0(weight0 == 0) = 1;
weight1(weight1 == 0) = 1;
if divideByOtf
    newb0 = newb0 ./ weight0;
    newb1 = newb1 ./ weight1;
end

[x, y] = meshgrid(1:w, 1:h);
radius = sqrt((y - cnt(1)).^2 + (x - cnt(2)).^2);
carrierRadius = sqrt(sum(carrierPixel .^ 2));
hasStructuredOtf = max(abs(otfValues(:)) - min(abs(otfValues(:)))) > eps;
if carrierRadius > eps && hasStructuredOtf
    ratio = radius ./ carrierRadius;
    annulusMask = ratio < dist | ratio > (1 - dist);
    newb0(annulusMask) = 0;

    shiftedMask = circshift(annulusMask, [round(carrierPixel(2)), round(carrierPixel(1))]);
    newb1(shiftedMask) = 0;
end
end

function shiftedOtf = shiftOtfByCarrier(otfValues, carrierPixel)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);
end
