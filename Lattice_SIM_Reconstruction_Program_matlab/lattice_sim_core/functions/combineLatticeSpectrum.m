function [SIM, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params)
%COMBINELATTICESPECTRUM Combine Lattice bands with OTF/Wiener damping.

[h, w] = size(bands.C0);
outputSize = [2*h, 2*w];
otfDoubleParams = params;
otfDoubleParams.pixelSizeNm = params.pixelSizeNm / 2;
otfDouble = buildLatticeOTF(outputSize(1), outputSize(2), otfDoubleParams);
bandDomain = getBandDomain(bands);

freq0 = placeSpectrumAtCenter(bandToFrequency(bands.C0, bandDomain), outputSize);
otf0 = otfDouble.values;
otf0Taper = double(abs(otf0) > 0);
otf0Mask = otf0Taper > 0;
fftDirectlyCombined = freq0 .* conj(otf0);
den = abs(otf0).^2;
blendDenominator = den;
bandWeightSum = den;
attenuatedBandWeightSum = den;
bandCoverageCount = double(otf0Mask);
bandContributionMasks = {otf0Mask};
bandTaperMasks = {otf0Taper};
sidebandPhaseReference = zeros(1, 4);
sidebandPhaseMagnitude = zeros(1, 4);
overlapPhaseMasks = cell(1, 4);

components = {
    bands.CsPlus,  carriers.ksPixel
    bands.CsMinus, -carriers.ksPixel
    bands.CtPlus,  carriers.ktPixel
    bands.CtMinus, -carriers.ktPixel
};

for idx = 1:size(components, 1)
    component = components{idx, 1};
    carrierPixel = components{idx, 2};
    freq = placeSpectrumAtCenter(bandToFrequency(component, bandDomain), outputSize);
    freq = shiftSpectrumOnCanvas(freq, -carrierPixel);
    shiftedOtf = shiftOtfByCarrier(otfDouble.values, carrierPixel);
    otfTaper = double(abs(shiftedOtf) > 0);
    otfMask = otfTaper > 0;
    phaseMask = overlapPhaseMask(freq0, freq, otf0Taper, otfTaper);
    phaseReference = sum(freq(phaseMask) .* conj(freq0(phaseMask)), 'all');
    if abs(phaseReference) <= eps
        centerRow = floor(outputSize(1)/2) + 1;
        centerCol = floor(outputSize(2)/2) + 1;
        phaseReference = freq(centerRow, centerCol) * conj(freq0(centerRow, centerCol));
    end
    freq = freq .* exp(-1i * angle(phaseReference));
    fftDirectlyCombined = fftDirectlyCombined + freq .* conj(shiftedOtf);
    bandWeight = abs(shiftedOtf).^2;
    den = den + bandWeight;
    blendDenominator = den;
    bandWeightSum = bandWeightSum + bandWeight;
    attenuatedBandWeightSum = attenuatedBandWeightSum + bandWeight;
    bandCoverageCount = bandCoverageCount + double(otfMask);
    bandContributionMasks{end + 1} = otfMask;
    bandTaperMasks{end + 1} = otfTaper;
    sidebandPhaseReference(idx) = angle(phaseReference);
    sidebandPhaseMagnitude(idx) = abs(phaseReference);
    overlapPhaseMasks{idx} = phaseMask;
end

wienerFilter = 1 ./ (blendDenominator + params.wiener);
combinedSpectrum = fftDirectlyCombined .* wienerFilter;

blendMax = max(blendDenominator(:));
supportMask = blendDenominator > blendMax * params.supportThreshold;
reliabilityMask = blendDenominator > blendMax * params.reliabilityThreshold;
finalSupportMask = supportMask;
finalConfidenceMask = double(supportMask);
combinedSpectrum = combinedSpectrum .* finalConfidenceMask;
[combinedSpectrum, apodizationMask] = applyLatticeApodization(combinedSpectrum, supportMask, params);

SIM = real(FFT2D(combinedSpectrum, true));

if any(~isfinite(SIM(:)))
    error('LatticeSIM:InvalidReconstruction', 'Reconstructed SIM image contains NaN or Inf values.');
end

diagnostics.simSpectrum = combinedSpectrum;
diagnostics.fftDirectlyCombined = fftDirectlyCombined;
diagnostics.wienerFilter = wienerFilter;
diagnostics.wienerDenominator = blendDenominator;
diagnostics.physicalOtfDenominator = den;
diagnostics.blendDenominator = blendDenominator;
diagnostics.physicalOtfValues = otfDouble.values;
diagnostics.bandWeightSum = bandWeightSum;
diagnostics.attenuatedBandWeightSum = attenuatedBandWeightSum;
diagnostics.bandCoverageCount = bandCoverageCount;
diagnostics.overlapTransitionMap = abs(del2(double(bandCoverageCount)));
diagnostics.supportMask = supportMask;
diagnostics.finalSupportMask = finalSupportMask;
diagnostics.reliabilityMask = reliabilityMask;
diagnostics.finalConfidenceMask = finalConfidenceMask;
diagnostics.apodizationMask = apodizationMask;
diagnostics.bandContributionMasks = bandContributionMasks;
diagnostics.bandTaperMasks = bandTaperMasks;
diagnostics.sidebandPhaseReference = sidebandPhaseReference;
diagnostics.sidebandPhaseMagnitude = sidebandPhaseMagnitude;
diagnostics.overlapPhaseMasks = overlapPhaseMasks;
diagnostics.outputImageMode = 'real';
diagnostics.bandDomain = bandDomain;
diagnostics.fusionMode = 'hifi-style-direct-wiener';
diagnostics.outputOtfPixelSizeNm = otfDoubleParams.pixelSizeNm;
diagnostics.modulationCompensationMode = 'separation-matrix';
end

function domain = getBandDomain(bands)
if isfield(bands, 'domain')
    domain = char(bands.domain);
else
    domain = 'space';
end
end

function frequency = bandToFrequency(component, domain)
if strcmp(domain, 'frequency')
    frequency = component;
else
    frequency = FFT2D(component, false);
end
end

function taper = smoothOtfTaper(otfValues, params)
maxOtf = max(otfValues(:));
low = maxOtf * params.otfTaperLow;
high = maxOtf * params.otfTaperHigh;
t = (otfValues - low) ./ (high - low + eps);
t = min(max(t, 0), 1);
taper = t .* t .* (3 - 2*t);
end

function mask = overlapPhaseMask(freq0, freq, taper0, taper1)
reliable = taper0 > 0.5 & taper1 > 0.5;
magnitude0 = abs(freq0);
magnitude1 = abs(freq);
reliable = reliable & magnitude0 > max(magnitude0(:)) * 1e-5 & ...
    magnitude1 > max(magnitude1(:)) * 1e-5;

[h, w] = size(freq0);
[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
radius = hypot(x, y);
reliable(radius < max(2, min(h, w) * 0.02)) = false;

if nnz(reliable) > 8
    mask = reliable;
    return;
end

mask = taper0 > 0.1 & taper1 > 0.1;
if nnz(mask) <= 8
    mask = taper0 > 0 & taper1 > 0;
end
end

function confidence = smoothConfidenceMask(denominator, denMax, params)
low = denMax * params.supportThreshold;
high = denMax * params.reliabilityThreshold;
t = (denominator - low) ./ (high - low + eps);
t = min(max(t, 0), 1);
confidence = t .* t .* (3 - 2*t);
end

function output = placeSpectrumAtCenter(spectrum, outputSize)
[h, w] = size(spectrum);
output = zeros(outputSize);
rowStart = floor((outputSize(1) - h) / 2) + 1;
colStart = floor((outputSize(2) - w) / 2) + 1;
output(rowStart:rowStart+h-1, colStart:colStart+w-1) = spectrum;
end

function shifted = shiftSpectrumOnCanvas(spectrum, shiftPixel)
[h, w] = size(spectrum);
[x, y] = meshgrid(0:w-1, 0:h-1);
image = FFT2D(spectrum, true);
phaseRamp = exp(2i*pi * (shiftPixel(1) * x / w + shiftPixel(2) * y / h));
shifted = FFT2D(image .* phaseRamp, false);
end

function shiftedOtf = shiftOtfByCarrier(otfValues, carrierPixel)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);
end
