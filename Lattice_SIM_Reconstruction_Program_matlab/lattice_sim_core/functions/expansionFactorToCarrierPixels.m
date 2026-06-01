function [ksPixel, ktPixel, carrierMagnitudePixels] = expansionFactorToCarrierPixels( ...
    expansionFactor, imageSize, pixelSizeNm, emissionWavelengthNm, NA)
%EXPANSIONFACTORTOCARRIERPIXELS Convert SIM expansion factor to carrier pixels.

if ~isnumeric(expansionFactor) || ~isscalar(expansionFactor) ...
        || ~isfinite(expansionFactor) || expansionFactor <= 1 || expansionFactor > 2
    error('LatticeSIM:InvalidExpansionFactor', ...
        'Expansion factor must be a finite scalar greater than 1 and no larger than 2.');
end

if ~isnumeric(imageSize) || numel(imageSize) ~= 2 || any(imageSize <= 0)
    error('LatticeSIM:InvalidExpansionFactor', 'imageSize must be a two-element positive vector.');
end
if ~isnumeric(pixelSizeNm) || ~isscalar(pixelSizeNm) || pixelSizeNm <= 0 ...
        || ~isnumeric(emissionWavelengthNm) || ~isscalar(emissionWavelengthNm) || emissionWavelengthNm <= 0 ...
        || ~isnumeric(NA) || ~isscalar(NA) || NA <= 0
    error('LatticeSIM:InvalidExpansionFactor', ...
        'pixelSizeNm, emissionWavelengthNm, and NA must be positive scalars.');
end

cutoffCyclesPerNm = 2 * NA / emissionWavelengthNm;
carrierCyclesPerNm = (expansionFactor - 1) * cutoffCyclesPerNm;

carrierMagnitudePixels = carrierCyclesPerNm * imageSize(2) * pixelSizeNm;
carrierMagnitudePixelsY = carrierCyclesPerNm * imageSize(1) * pixelSizeNm;

ksPixel = [round(carrierMagnitudePixels), 0];
ktPixel = [0, round(carrierMagnitudePixelsY)];
end
