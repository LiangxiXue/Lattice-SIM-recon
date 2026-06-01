function summary = run_component_separation_debug()
%RUN_COMPONENT_SEPARATION_DEBUG Isolate whether five Lattice-SIM bands demodulate correctly.
%
% This script intentionally stops before final frequency-spectrum stitching.
% It uses the same testpat.tiff-driven simulation path as
% simulate_testpat_lattice_sim.m, compares separateLatticeBands output with
% analytically known simulated components, and displays the five separated
% component spectra.

debugDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(debugDir);
functionsDir = fullfile(coreDir, 'functions');
addpath(functionsDir);

[objectImage, simParams] = loadTestpatSimulationInput(coreDir);

cases = {
    makeCase('ideal_no_otf_no_noise_no_norm', false, 0, 0, false, false)
    makeCase('otf_no_noise_no_norm', true, 0, 0, false, false)
    makeCase('current_script_like_with_norm', true, 0.05, 0.01, true, true)
};

summary = struct([]);
for caseIdx = 1:numel(cases)
    caseSpec = cases{caseIdx};

    currentSimParams = simParams;
    currentSimParams.useOTF = caseSpec.useOTF;
    currentSimParams.noiseLevel = caseSpec.noiseLevel;
    currentSimParams.phaseErrorStd = caseSpec.phaseErrorStd;

    [rawStack, truth] = simulateLatticeSIMExperiment(objectImage, currentSimParams);

    reconParams = defaultLatticeSIMParams();
    reconParams.pixelSizeNm = currentSimParams.pixelSizeNm;
    reconParams.emissionWavelengthNm = currentSimParams.emissionWavelengthNm;
    reconParams.NA = currentSimParams.NA;
    reconParams.modulationS = currentSimParams.modulationS;
    reconParams.modulationT = currentSimParams.modulationT;
    reconParams.normalizeFrames = caseSpec.normalizeFrames;

    [stackForSeparation, ~] = normalizeSIMFrames(rawStack, reconParams);
    bands = separateLatticeBands(stackForSeparation, reconParams);
    expectedBands = nominalBandsFromSimulationTruth(truth);
    metrics = compareBands(bands, expectedBands);

    try
        carriers = estimateLatticeCarrier(bands, reconParams);
        metrics.ksErrorPixels = norm(carriers.ksPixel - truth.ksPixel);
        metrics.ktErrorPixels = norm(carriers.ktPixel - truth.ktPixel);
        metrics.estimatedKsPixel = carriers.ksPixel;
        metrics.estimatedKtPixel = carriers.ktPixel;
    catch carrierError
        fprintf('Carrier estimation failed: %s\n', carrierError.message);
        metrics.ksErrorPixels = NaN;
        metrics.ktErrorPixels = NaN;
        metrics.estimatedKsPixel = [NaN, NaN];
        metrics.estimatedKtPixel = [NaN, NaN];
    end
    metrics.truthKsPixel = truth.ksPixel;
    metrics.truthKtPixel = truth.ktPixel;

    summary(caseIdx).name = caseSpec.name; %#ok<AGROW>
    summary(caseIdx).metrics = metrics; %#ok<AGROW>

    fprintf('\n[%s]\n', caseSpec.name);
    printMetrics(metrics);
    if caseSpec.useOTF && ~caseSpec.normalizeFrames && caseSpec.noiseLevel == 0 && caseSpec.phaseErrorStd == 0
        diagnoseCombineVariants(bands, truth, reconParams);
    end
    if caseSpec.showSpectrum
        showSeparatedSpectra(bands, caseSpec.name, metrics);
    end
end
end

function caseSpec = makeCase(name, useOTF, noiseLevel, phaseErrorStd, normalizeFrames, showSpectrum)
caseSpec.name = name;
caseSpec.useOTF = useOTF;
caseSpec.noiseLevel = noiseLevel;
caseSpec.phaseErrorStd = phaseErrorStd;
caseSpec.normalizeFrames = normalizeFrames;
caseSpec.showSpectrum = showSpectrum;
end

function [objectImage, simParams] = loadTestpatSimulationInput(coreDir)
inputPath = fullfile(coreDir, 'testpat.tiff');
if exist(inputPath, 'file') ~= 2
    error('LatticeSIM:MissingInputImage', 'Expected test image not found: %s', inputPath);
end

objectImageFull = double(imread(inputPath));
cropSize = 256;
if size(objectImageFull, 1) < cropSize || size(objectImageFull, 2) < cropSize
    error('LatticeSIM:InvalidInputImage', 'testpat.tiff must be at least %d x %d pixels.', cropSize, cropSize);
end

rowStart = floor((size(objectImageFull, 1) - cropSize) / 2) + 1;
colStart = floor((size(objectImageFull, 2) - cropSize) / 2) + 1;
objectImage = objectImageFull(rowStart:rowStart+cropSize-1, colStart:colStart+cropSize-1);

simParams.imageSize = size(objectImage);
simParams.expansionFactor = 1.6;
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 97.5;
simParams.emissionWavelengthNm = 561;
simParams.NA = 1.42;
simParams.randomSeed = 7;
[simParams.ksPixel, simParams.ktPixel] = expansionFactorToCarrierPixels( ...
    simParams.expansionFactor, simParams.imageSize, simParams.pixelSizeNm, ...
    simParams.emissionWavelengthNm, simParams.NA);
end

function expectedBands = nominalBandsFromSimulationTruth(truth)
[h, w] = size(truth.object);
[x, y] = meshgrid(0:w-1, 0:h-1);
phaseS = 2*pi*(truth.ksPixel(1) * x / w + truth.ksPixel(2) * y / h);
phaseT = 2*pi*(truth.ktPixel(1) * x / w + truth.ktPixel(2) * y / h);

params = truth.params;
expectedBands.C0 = params.meanIllumination .* truth.object;
expectedBands.CsPlus = (params.modulationS / 2) .* truth.object .* exp(1i * phaseS);
expectedBands.CsMinus = (params.modulationS / 2) .* truth.object .* exp(-1i * phaseS);
expectedBands.CtPlus = (params.modulationT / 2) .* truth.object .* exp(1i * phaseT);
expectedBands.CtMinus = (params.modulationT / 2) .* truth.object .* exp(-1i * phaseT);

if params.useOTF
    expectedBands.C0 = applySimulationOtf(expectedBands.C0, truth.otf);
    expectedBands.CsPlus = applySimulationOtf(expectedBands.CsPlus, truth.otf);
    expectedBands.CsMinus = applySimulationOtf(expectedBands.CsMinus, truth.otf);
    expectedBands.CtPlus = applySimulationOtf(expectedBands.CtPlus, truth.otf);
    expectedBands.CtMinus = applySimulationOtf(expectedBands.CtMinus, truth.otf);
end
end

function image = applySimulationOtf(image, otf)
image = ifft2(fft2(image) .* fftshift(otf));
end

function metrics = compareBands(bands, expectedBands)
names = {'C0', 'CsPlus', 'CsMinus', 'CtPlus', 'CtMinus'};
for idx = 1:numel(names)
    name = names{idx};
    actual = bands.(name);
    expected = expectedBands.(name);
    diff = actual - expected;
    denom = max(norm(expected(:)), eps);
    metrics.(name).maxAbsError = max(abs(diff(:)));
    metrics.(name).relativeL2Error = norm(diff(:)) / denom;
    metrics.(name).actualNorm = norm(actual(:));
    metrics.(name).expectedNorm = norm(expected(:));
end

metrics.maxBandRelativeL2Error = max(cellfun(@(name) metrics.(name).relativeL2Error, names));
metrics.maxBandAbsError = max(cellfun(@(name) metrics.(name).maxAbsError, names));
end

function showSeparatedSpectra(bands, caseName, metrics)
names = {'C0', 'CsPlus', 'CsMinus', 'CtPlus', 'CtMinus'};
displayNames = {'C0', 'Cs+', 'Cs-', 'Ct+', 'Ct-'};

figure('Name', ['Separated spectra: ' caseName], 'Color', 'w');
for idx = 1:numel(names)
    name = names{idx};
    subplot(1, numel(names), idx);
    spectrumView = log1p(abs(centeredFft2(bands.(name))));
    imagesc(spectrumView);
    axis image off;
    colormap(gca, 'hot');
    title(sprintf('%s\\nrelL2 %.2e', displayNames{idx}, ...
        metrics.(name).relativeL2Error), 'Interpreter', 'none');
end
sgtitle(sprintf('%s | ks err %.3g px, kt err %.3g px', ...
    caseName, metrics.ksErrorPixels, metrics.ktErrorPixels), 'Interpreter', 'none');
drawnow;
end

function spectrum = centeredFft2(image)
spectrum = fftshift(fft2(ifftshift(image)));
end

function diagnoseCombineVariants(bands, truth, params)
fprintf('\n[combine variant check against ideal object spectrum]\n');

otf = buildLatticeOTF(size(bands.C0, 1), size(bands.C0, 2), params);
truthCarriers.ksPixel = truth.ksPixel;
truthCarriers.ktPixel = truth.ktPixel;
truthCarriers.ksRadPerPixel = [2*pi*truth.ksPixel(1)/size(bands.C0, 2), ...
    2*pi*truth.ksPixel(2)/size(bands.C0, 1)];
truthCarriers.ktRadPerPixel = [2*pi*truth.ktPixel(1)/size(bands.C0, 2), ...
    2*pi*truth.ktPixel(2)/size(bands.C0, 1)];

variants = {
    makeCombineVariant('current_code_equivalent', -1, 1, false, 0.5)
    makeCombineVariant('old_divide_by_modulation', -1, 1, false, 1)
    makeCombineVariant('opposite_phase_ramp', 1, 1, false, 0.5)
    makeCombineVariant('opposite_otf_shift', -1, -1, false, 0.5)
    makeCombineVariant('conjugate_sidebands', -1, 1, true, 0.5)
};

objectSpectrum = centeredFft2(truth.object);
for idx = 1:numel(variants)
    variant = variants{idx};
    [combinedSpectrum, denominator] = combineVariantSpectrum(bands, truthCarriers, otf.values, params, variant);
    metrics = compareSpectrumToTruth(combinedSpectrum, objectSpectrum, denominator);
    fprintf('%-25s relL2=% .6e  rawRelL2=% .6e  phaseRMSE=% .6e  supportCorr=% .6f\n', ...
        variant.name, metrics.relativeL2Error, metrics.rawRelativeL2Error, ...
        metrics.phaseRmseRad, metrics.supportCorrelation);
end

defaultModulationParams = params;
defaultModulationParams.modulationS = 0.7;
defaultModulationParams.modulationT = 0.7;
defaultVariant = makeCombineVariant('main_script_modulation_0p7', -1, 1, false, 1);
[combinedSpectrum, denominator] = combineVariantSpectrum(bands, truthCarriers, otf.values, ...
    defaultModulationParams, defaultVariant);
metrics = compareSpectrumToTruth(combinedSpectrum, objectSpectrum, denominator);
fprintf('%-25s relL2=% .6e  rawRelL2=% .6e  phaseRMSE=% .6e  supportCorr=% .6f\n', ...
    defaultVariant.name, metrics.relativeL2Error, metrics.rawRelativeL2Error, ...
    metrics.phaseRmseRad, metrics.supportCorrelation);
end

function variant = makeCombineVariant(name, phaseRampSign, otfShiftSign, conjugateSidebands, modulationFactor)
variant.name = name;
variant.phaseRampSign = phaseRampSign;
variant.otfShiftSign = otfShiftSign;
variant.conjugateSidebands = conjugateSidebands;
variant.modulationFactor = modulationFactor;
end

function [combinedSpectrum, den] = combineVariantSpectrum(bands, carriers, otfValues, params, variant)
[h, w] = size(bands.C0);
[x, y] = meshgrid(0:w-1, 0:h-1);

freq0 = centeredFft2(bands.C0);
num = freq0 .* conj(otfValues);
den = abs(otfValues).^2;

components = {
    bands.CsPlus,  carriers.ksRadPerPixel,  carriers.ksPixel, params.modulationS
    bands.CsMinus, -carriers.ksRadPerPixel, -carriers.ksPixel, params.modulationS
    bands.CtPlus,  carriers.ktRadPerPixel,  carriers.ktPixel, params.modulationT
    bands.CtMinus, -carriers.ktRadPerPixel, -carriers.ktPixel, params.modulationT
};

for componentIdx = 1:size(components, 1)
    component = components{componentIdx, 1};
    carrierRad = components{componentIdx, 2};
    carrierPixel = components{componentIdx, 3};
    modulation = components{componentIdx, 4} * variant.modulationFactor;

    if variant.conjugateSidebands
        component = conj(component);
    end

    phaseRamp = exp(1i * variant.phaseRampSign * ...
        (carrierRad(1) * x + carrierRad(2) * y));
    centered = component .* phaseRamp;
    freq = centeredFft2(centered);
    shiftedOtf = shiftOtfForVariant(otfValues, carrierPixel, variant.otfShiftSign);

    centerRow = floor(h/2) + 1;
    centerCol = floor(w/2) + 1;
    phaseReference = freq(centerRow, centerCol) / (freq0(centerRow, centerCol) + eps);
    freq = freq .* exp(-1i * angle(phaseReference));

    num = num + freq .* conj(shiftedOtf) ./ modulation;
    den = den + abs(shiftedOtf).^2;
end

combinedSpectrum = num ./ (den + params.wiener);
end

function shiftedOtf = shiftOtfForVariant(otfValues, carrierPixel, shiftSign)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + shiftSign * carrierPixel(1), ...
    rowGrid + shiftSign * carrierPixel(2), 'linear', 0);
end

function metrics = compareSpectrumToTruth(combinedSpectrum, objectSpectrum, denominator)
supportMask = denominator > 0;
truthMagnitude = abs(objectSpectrum);
supportMask = supportMask & truthMagnitude > max(truthMagnitude(:)) * 1e-8;
truthOnSupport = objectSpectrum(supportMask);
combinedOnSupport = combinedSpectrum(supportMask);

scale = (combinedOnSupport(:)' * truthOnSupport(:)) / ...
    (combinedOnSupport(:)' * combinedOnSupport(:) + eps);
scaledCombined = combinedSpectrum .* scale;
diff = scaledCombined(supportMask) - truthOnSupport;
rawDiff = combinedOnSupport - truthOnSupport;

metrics.relativeL2Error = norm(diff(:)) / max(norm(truthOnSupport(:)), eps);
metrics.rawRelativeL2Error = norm(rawDiff(:)) / max(norm(truthOnSupport(:)), eps);
phaseDiff = angle(scaledCombined(supportMask) .* conj(truthOnSupport));
metrics.phaseRmseRad = sqrt(mean(phaseDiff(:).^2));
metrics.supportCorrelation = abs(truthOnSupport(:)' * combinedOnSupport(:)) / ...
    max(norm(truthOnSupport(:)) * norm(combinedOnSupport(:)), eps);
end

function printMetrics(metrics)
names = {'C0', 'CsPlus', 'CsMinus', 'CtPlus', 'CtMinus'};
for idx = 1:numel(names)
    name = names{idx};
    fprintf('%-8s maxAbs=% .6e  relL2=% .6e\n', name, ...
        metrics.(name).maxAbsError, metrics.(name).relativeL2Error);
end
fprintf('ks truth=[%.3f %.3f], estimated=[%.3f %.3f], error=%.6g px\n', ...
    metrics.truthKsPixel(1), metrics.truthKsPixel(2), ...
    metrics.estimatedKsPixel(1), metrics.estimatedKsPixel(2), metrics.ksErrorPixels);
fprintf('kt truth=[%.3f %.3f], estimated=[%.3f %.3f], error=%.6g px\n', ...
    metrics.truthKtPixel(1), metrics.truthKtPixel(2), ...
    metrics.estimatedKtPixel(1), metrics.estimatedKtPixel(2), metrics.ktErrorPixels);
end
