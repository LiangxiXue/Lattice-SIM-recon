function [SIM, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params)
%COMBINELATTICESPECTRUM Combine Lattice bands with OTF/Wiener damping.

[h, w] = size(bands.C0);
outputSize = [2*h, 2*w];
otfDoubleParams = params;
otfDoubleParams.pixelSizeNm = params.pixelSizeNm / 2;
otfDouble = buildLatticeOTF(outputSize(1), outputSize(2), otfDoubleParams);
effectiveOtfValues = getEffectiveOtfValues(otfDouble);
nonAttenuatedOtfValues = getNonAttenuatedOtfValues(otfDouble);
bandDomain = getBandDomain(bands);

freq0 = placeSpectrumAtCenter(bandToFrequency(bands.C0, bandDomain), outputSize);
[peakPixel, centroidPixel] = spectrumPositionPixels(freq0);
shiftedBandNames = {'C0', 'CsPlus', 'CsMinus', 'CtPlus', 'CtMinus'};
shiftedBandSpectra = {freq0};
shiftedBandCarrierPixels = zeros(5, 2);
shiftedBandShiftPixels = zeros(5, 2);
shiftedBandPeakPixels = zeros(5, 2);
shiftedBandCentroidPixels = zeros(5, 2);
shiftedBandPeakPixels(1, :) = peakPixel;
shiftedBandCentroidPixels(1, :) = centroidPixel;
otf0 = effectiveOtfValues;
otf0Taper = double(abs(otf0) > 0);
otf0Mask = otf0Taper > 0;
fusionNumeratorContributions = cell(1, 5);
fusionNumeratorContributions{1} = freq0 .* conj(otf0);
fftDirectlyCombined = fusionNumeratorContributions{1};
den = abs(otf0).^2;
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
    shiftedOtf = shiftOtfByCarrier(effectiveOtfValues, carrierPixel);
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
    [peakPixel, centroidPixel] = spectrumPositionPixels(freq);
    componentSlot = idx + 1;
    shiftedBandSpectra{componentSlot} = freq;
    shiftedBandCarrierPixels(componentSlot, :) = carrierPixel;
    shiftedBandShiftPixels(componentSlot, :) = -carrierPixel;
    shiftedBandPeakPixels(componentSlot, :) = peakPixel;
    shiftedBandCentroidPixels(componentSlot, :) = centroidPixel;
    contribution = freq .* conj(shiftedOtf);
    fusionNumeratorContributions{componentSlot} = contribution;
    fftDirectlyCombined = fftDirectlyCombined + contribution;
    bandWeight = abs(shiftedOtf).^2;
    den = den + bandWeight;
    bandWeightSum = bandWeightSum + bandWeight;
    attenuatedBandWeightSum = attenuatedBandWeightSum + bandWeight;
    bandCoverageCount = bandCoverageCount + double(otfMask);
    bandContributionMasks{end + 1} = otfMask;
    bandTaperMasks{end + 1} = otfTaper;
    sidebandPhaseReference(idx) = angle(phaseReference);
    sidebandPhaseMagnitude(idx) = abs(phaseReference);
    overlapPhaseMasks{idx} = phaseMask;
end

[wienerDenominatorW1Base, wienerDenominatorW2Base, denominatorDiagnostics] = ...
    buildHifiStyleDenominators(otfDouble, carriers, params);
blendDenominator = wienerDenominatorW1Base;
wienerFilter = 1 ./ (blendDenominator + params.wiener);

blendMax = max(blendDenominator(:));
supportMask = blendDenominator > blendMax * params.supportThreshold;
reliabilityMask = blendDenominator > blendMax * params.reliabilityThreshold;
finalConfidenceMask = smoothConfidenceMask(blendDenominator, blendMax, params);
finalSupportMask = finalConfidenceMask > 0;

fusionMode = char(params.fusionMode);
hifiMask = double(supportMask);
hifiOtf = otfDouble.values;
wienerDenominatorW1 = wienerDenominatorW1Base;
wienerDenominatorW2 = wienerDenominatorW2Base;
wienerW1 = [];
wienerW2 = [];
spectrumAfterW1 = [];
spectrumAfterW2 = [];

switch fusionMode
    case 'single-step'
        combinedSpectrum = fftDirectlyCombined .* wienerFilter;
        combinedSpectrum = combinedSpectrum .* finalConfidenceMask;
        [combinedSpectrum, apodizationMask] = applyLatticeApodization(combinedSpectrum, supportMask, params);
    case 'hifi-two-step'
        hifiOtf = makeExtendedHifiOtf(supportMask);
        apodizationMask = makeLatticeApodizationMask(supportMask, params);
        wienerDenominatorW1 = wienerDenominatorW1Base * params.hifiDenominatorScaleW1;
        wienerDenominatorW2 = wienerDenominatorW2Base * params.hifiDenominatorScaleW2;
        wienerW1 = hifiOtf ./ (wienerDenominatorW1 + params.wienerW1 ^ 2);
        spectrumAfterW1 = fftDirectlyCombined .* wienerW1 .* hifiMask;
        wienerW2 = apodizationMask ./ (wienerDenominatorW2 + params.wienerW2 ^ 2);
        spectrumAfterW2 = spectrumAfterW1 .* wienerW2 .* hifiMask;
        combinedSpectrum = spectrumAfterW2;
    otherwise
        error('LatticeSIM:InvalidFusionMode', ...
            'fusionMode must be "single-step" or "hifi-two-step".');
end

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
diagnostics.effectiveOtfValues = effectiveOtfValues;
diagnostics.nonAttenuatedOtfValues = nonAttenuatedOtfValues;
diagnostics.reconstructionFilter = otfDouble.reconstructionFilter;
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
diagnostics.shiftedBandNames = shiftedBandNames;
diagnostics.shiftedBandSpectra = shiftedBandSpectra;
diagnostics.shiftedBandCarrierPixels = shiftedBandCarrierPixels;
diagnostics.shiftedBandShiftPixels = shiftedBandShiftPixels;
diagnostics.shiftedBandPeakPixels = shiftedBandPeakPixels;
diagnostics.shiftedBandCentroidPixels = shiftedBandCentroidPixels;
diagnostics.fusionNumeratorContributions = fusionNumeratorContributions;
diagnostics.outputImageMode = 'real';
diagnostics.bandDomain = bandDomain;
diagnostics.fusionMode = fusionMode;
diagnostics.outputOtfPixelSizeNm = otfDoubleParams.pixelSizeNm;
diagnostics.modulationCompensationMode = 'separation-matrix';
diagnostics.spectrumAfterW1 = spectrumAfterW1;
diagnostics.spectrumAfterW2 = spectrumAfterW2;
diagnostics.wienerW1 = wienerW1;
diagnostics.wienerW2 = wienerW2;
diagnostics.wienerDenominatorW1 = wienerDenominatorW1;
diagnostics.wienerDenominatorW2 = wienerDenominatorW2;
diagnostics.wienerDenominatorW1Base = wienerDenominatorW1Base;
diagnostics.wienerDenominatorW2Base = wienerDenominatorW2Base;
diagnostics.denominatorBandNames = denominatorDiagnostics.bandNames;
diagnostics.denominatorBaseOtfValues = denominatorDiagnostics.baseOtfValues;
diagnostics.denominatorAttenuationMasksW1 = denominatorDiagnostics.attenuationMasksW1;
diagnostics.denominatorAttenuationMasksW2 = denominatorDiagnostics.attenuationMasksW2;
diagnostics.denominatorW1Contributions = denominatorDiagnostics.w1Contributions;
diagnostics.denominatorW2Contributions = denominatorDiagnostics.w2Contributions;
diagnostics.hifiMask = hifiMask;
diagnostics.hifiOtf = hifiOtf;
end

function values = getEffectiveOtfValues(otf)
if isfield(otf, 'reconstructionFilter')
    values = otf.values .* otf.reconstructionFilter;
else
    values = otf.values;
end
end

function values = getNonAttenuatedOtfValues(otf)
values = otf.values;
if isfield(otf, 'empiricalDampingMask')
    values = values .* otf.empiricalDampingMask;
end
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

function hifiOtf = makeExtendedHifiOtf(supportMask)
[h, w] = size(supportMask);
[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
radius = hypot(x, y);
supportRadius = max(radius(supportMask), [], 'all');
cutoffRadius = max(1, supportRadius + 1);
rho = radius ./ cutoffRadius;
hifiOtf = zeros(h, w);
inside = rho <= 1;
hifiOtf(inside) = (2 / pi) * (acos(rho(inside)) - ...
    rho(inside) .* sqrt(1 - rho(inside) .^ 2));
hifiOtf = hifiOtf .* double(supportMask);
end

function apo = makeLatticeApodizationMask(supportMask, params)
apo = double(supportMask);
if ~isfield(params, 'apodizationEnabled') || ~params.apodizationEnabled
    return;
end

if isfield(params, 'apodizationMode') && strcmp(char(params.apodizationMode), 'radial-gaussian')
    [h, w] = size(supportMask);
    [x, y] = meshgrid(1:w, 1:h);
    center = [floor(h/2) + 1, floor(w/2) + 1];
    radius = hypot((x - center(2)) * 2 / w, (y - center(1)) * 2 / h);
    apo = exp(-0.5 * (radius ./ params.apodizationRadius .* sqrt(2*log(2))).^2);
    apo = apo .* double(supportMask);
elseif exist('bwdist', 'file') == 2
    dist = bwdist(~supportMask);
    if max(dist(:)) > 0
        apo = (dist ./ max(dist(:))) .^ params.apodizationStrength;
    end
end
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

function [w1Base, w2Base, diagnostics] = buildHifiStyleDenominators(otf, carriers, params)
baseOtf = getNonAttenuatedOtfValues(otf);
centerWeight = getParam(params, 'hifiCenterDenominatorWeight', 0.5);
w1SideScale = getParam(params, 'hifiW1SidebandAttenuationScale', 1.0);
w2CenterScale = getParam(params, 'hifiW2CenterAttenuationScale', 1.05);
w2SideScale = getParam(params, 'hifiW2SidebandAttenuationScale', 1.15);

componentNames = {'C0+', 'C0-', 'CsPlus', 'CsMinus', 'CtPlus', 'CtMinus'};
componentCarriers = {
    [0, 0]
    [0, 0]
    carriers.ksPixel
    -carriers.ksPixel
    carriers.ktPixel
    -carriers.ktPixel
};
componentWeights = [centerWeight, centerWeight, 1, 1, 1, 1];

centerAttenuationW1 = ones(size(baseOtf));
centerAttenuationW2 = hifiAttenuationMaskForOtf(otf, params, w2CenterScale);
sideAttenuationW1 = hifiAttenuationMaskForOtf(otf, params, w1SideScale);
sideAttenuationW2 = hifiAttenuationMaskForOtf(otf, params, w2SideScale);

w1Base = zeros(size(baseOtf));
w2Base = zeros(size(baseOtf));
diagnostics.bandNames = componentNames;
diagnostics.baseOtfValues = cell(1, numel(componentNames));
diagnostics.attenuationMasksW1 = cell(1, numel(componentNames));
diagnostics.attenuationMasksW2 = cell(1, numel(componentNames));
diagnostics.w1Contributions = cell(1, numel(componentNames));
diagnostics.w2Contributions = cell(1, numel(componentNames));

for idx = 1:numel(componentNames)
    carrierPixel = componentCarriers{idx};
    shiftedBaseOtf = shiftOtfByCarrier(baseOtf, carrierPixel);
    if idx <= 2
        attenuationW1 = centerAttenuationW1;
        attenuationW2 = centerAttenuationW2;
    else
        attenuationW1 = shiftOtfByCarrier(sideAttenuationW1, carrierPixel);
        attenuationW2 = shiftOtfByCarrier(sideAttenuationW2, carrierPixel);
    end

    basePower = abs(shiftedBaseOtf).^2;
    contributionW1 = componentWeights(idx) .* basePower .* attenuationW1;
    contributionW2 = componentWeights(idx) .* basePower .* attenuationW2;
    w1Base = w1Base + contributionW1;
    w2Base = w2Base + contributionW2;

    diagnostics.baseOtfValues{idx} = shiftedBaseOtf;
    diagnostics.attenuationMasksW1{idx} = attenuationW1;
    diagnostics.attenuationMasksW2{idx} = attenuationW2;
    diagnostics.w1Contributions{idx} = contributionW1;
    diagnostics.w2Contributions{idx} = contributionW2;
end
end

function mask = hifiAttenuationMaskForOtf(otf, params, strengthScale)
mask = ones(size(otf.values));
if ~isfield(params, 'otfAttenuationEnabled') || ~params.otfAttenuationEnabled || ...
        params.otfAttenuationStrength <= 0
    return;
end

radius = hypot(otf.fxCyclesPerNm, otf.fyCyclesPerNm);
fwhmCyclesPerNm = params.otfAttenuationFwhm * 1e-3;
strength = params.otfAttenuationStrength ./ strengthScale;
mask = 1 - strength .* exp(-(radius .^ 2) ./ ((0.5 * fwhmCyclesPerNm) ^ 2));
mask = min(max(mask, 0), 1);
end

function value = getParam(params, name, defaultValue)
if isfield(params, name)
    value = params.(name);
else
    value = defaultValue;
end
end

function [peakPixel, centroidPixel] = spectrumPositionPixels(spectrum)
[h, w] = size(spectrum);
[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
magnitude = abs(spectrum);
[~, maxIdx] = max(magnitude(:));
peakPixel = [x(maxIdx), y(maxIdx)];

weights = log10(1 + magnitude);
weightSum = sum(weights(:));
if weightSum > 0
    centroidPixel = [sum(x(:) .* weights(:)) / weightSum, ...
        sum(y(:) .* weights(:)) / weightSum];
else
    centroidPixel = [NaN, NaN];
end
end
