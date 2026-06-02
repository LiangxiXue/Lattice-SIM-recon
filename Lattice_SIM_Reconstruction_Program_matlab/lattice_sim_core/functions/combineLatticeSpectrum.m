function [SIM, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params)
%COMBINELATTICESPECTRUM Combine Lattice bands with OTF/Wiener damping.

[h, w] = size(bands.C0);
outputSize = [2*h, 2*w];
otfDouble = buildLatticeOTF(outputSize(1), outputSize(2), params);
bandDomain = getBandDomain(bands);

freq0 = placeSpectrumAtCenter(bandToFrequency(bands.C0, bandDomain), outputSize);
otf0 = otfDouble.values;
otf0Taper = smoothOtfTaper(abs(otf0), params);
otf0Mask = otf0Taper > 0;
attenuation0 = otfAttenuationMask(otf0, params);
centralBlend = otf0Taper .* attenuation0;
num = freq0 .* otf0Taper .* conj(otf0) .* attenuation0;
den = abs(otf0).^2;
blendDenominator = abs(otf0).^2 .* centralBlend;
bandWeightSum = abs(otf0).^2 .* otf0Taper;
attenuatedBandWeightSum = blendDenominator;
bandCoverageCount = double(otf0Mask);
bandContributionMasks = {otf0Mask};
bandTaperMasks = {otf0Taper};
sidebandPhaseReference = zeros(1, 4);
sidebandPhaseMagnitude = zeros(1, 4);
overlapPhaseMasks = cell(1, 4);

components = {
    bands.CsPlus,  carriers.ksPixel, params.modulationS
    bands.CsMinus, -carriers.ksPixel, params.modulationS
    bands.CtPlus,  carriers.ktPixel, params.modulationT
    bands.CtMinus, -carriers.ktPixel, params.modulationT
};

for idx = 1:size(components, 1)
    component = components{idx, 1};
    carrierPixel = components{idx, 2};
    sidebandAmplitude = components{idx, 3} / 2;
    freq = placeSpectrumAtCenter(bandToFrequency(component, bandDomain), outputSize);
    freq = shiftSpectrumOnCanvas(freq, -carrierPixel);
    shiftedOtf = shiftOtfByCarrier(otfDouble.values, -carrierPixel);
    otfTaper = smoothOtfTaper(abs(shiftedOtf), params);
    otfMask = otfTaper > 0;
    phaseMask = overlapPhaseMask(freq0, freq, otf0Taper, otfTaper);
    phaseReference = sum(freq(phaseMask) .* conj(freq0(phaseMask)), 'all');
    if abs(phaseReference) <= eps
        centerRow = floor(outputSize(1)/2) + 1;
        centerCol = floor(outputSize(2)/2) + 1;
        phaseReference = freq(centerRow, centerCol) * conj(freq0(centerRow, centerCol));
    end
    freq = freq .* exp(-1i * angle(phaseReference));
    attenuation = otfAttenuationMask(shiftedOtf, params);
    bandBlend = otfTaper .* attenuation;
    num = num + freq .* otfTaper .* conj(shiftedOtf) .* attenuation ./ sidebandAmplitude;
    den = den + abs(shiftedOtf).^2;
    blendDenominator = blendDenominator + abs(shiftedOtf).^2 .* bandBlend;
    bandWeight = abs(shiftedOtf).^2 .* otfTaper;
    bandWeightSum = bandWeightSum + bandWeight;
    attenuatedBandWeightSum = attenuatedBandWeightSum + abs(shiftedOtf).^2 .* bandBlend;
    bandCoverageCount = bandCoverageCount + double(otfMask);
    bandContributionMasks{end + 1} = otfMask;
    bandTaperMasks{end + 1} = otfTaper;
    sidebandPhaseReference(idx) = angle(phaseReference);
    sidebandPhaseMagnitude(idx) = abs(phaseReference);
    overlapPhaseMasks{idx} = phaseMask;
end

combinedSpectrum = num ./ (blendDenominator + params.wiener);

blendMax = max(blendDenominator(:));
supportMask = blendDenominator > blendMax * params.supportThreshold;
reliabilityMask = blendDenominator > blendMax * params.reliabilityThreshold;
finalSupportMask = supportMask;
finalConfidenceMask = smoothConfidenceMask(blendDenominator, blendMax, params);
combinedSpectrum = combinedSpectrum .* finalConfidenceMask;
[combinedSpectrum, apodizationMask] = applyLatticeApodization(combinedSpectrum, finalConfidenceMask > 0, params);

SIM = real(ifft2c(combinedSpectrum));

if any(~isfinite(SIM(:)))
    error('LatticeSIM:InvalidReconstruction', 'Reconstructed SIM image contains NaN or Inf values.');
end

diagnostics.simSpectrum = combinedSpectrum;
diagnostics.wienerDenominator = blendDenominator;
diagnostics.physicalOtfDenominator = den;
diagnostics.blendDenominator = blendDenominator;
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
    frequency = fft2c(component);
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

function attenuation = otfAttenuationMask(otfValues, params)
attenuation = ones(size(otfValues));
if ~isfield(params, 'otfAttenuationEnabled') || ~params.otfAttenuationEnabled
    return;
end

[~, maxIdx] = max(abs(otfValues(:)));
[centerRow, centerCol] = ind2sub(size(otfValues), maxIdx);
[h, w] = size(otfValues);
[x, y] = meshgrid(1:w, 1:h);
radius = hypot((x - centerCol) * 2 / w, (y - centerRow) * 2 / h);
attenuation = 1 - params.otfAttenuationStrength .* ...
    exp(-(radius .^ 2) ./ ((0.5 * params.otfAttenuationFwhm) ^ 2));
attenuation = min(max(attenuation, 0), 1);
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
image = ifft2c(spectrum);
phaseRamp = exp(2i*pi * (shiftPixel(1) * x / w + shiftPixel(2) * y / h));
shifted = fft2c(image .* phaseRamp);
end

function shiftedOtf = shiftOtfByCarrier(otfValues, carrierPixel)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);
end
