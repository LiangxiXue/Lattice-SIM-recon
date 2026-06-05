scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
coreDir = fileparts(fileparts(scriptDir));
functionsDir = fullfile(coreDir, 'functions');
addpath(genpath(functionsDir));

inputPath = fullfile(coreDir, 'roi.tif');
resultPath = fullfile(coreDir, 'output', 'diagnostics', 'result.mat');

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 610;
params.NA = 1.42;
params.phasePairs = [
    0,       0
    2*pi/3, 0
    4*pi/3, 0
    0,       2*pi/3
    2*pi/3, 4*pi/3
];
params.phaseMatrix = [];
params.normalizeFrames = true;
params.preprocessingMode = "hifi-rl-fft";
params.deconvolutionIterations = 5;
params.separationFunction = "separateLatticeBandsFrequency";
params.separationInputDomain = "frequency";
params.enableLatticeParameterEstimation = true;
params.latticeCarrierRefinementIterations = 3;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCorrelationOverlap = 0.15;
params.carrierSearchMode = "unconstrained";
params.carrierAxisToleranceDeg = 15;
params.wiener = 0.04;
params.apodizationEnabled = true;
params.apodizationMode = "radial-gaussian";
params.apodizationStrength = 0.4;
params.apodizationRadius = 0.5;
params.supportThreshold = 1e-4;
params.reliabilityThreshold = 1e-3;
params.otfClipThreshold = 1e-4;
params.otfTaperLow = 1e-4;
params.otfTaperHigh = 5e-3;
params.otfAttenuationEnabled = true;
params.otfAttenuationStrength = 0.15;
params.otfAttenuationFwhm = 0.25;
params.notchScale = 0.5;
params.outputScaleMode = "none";
params = defaultLatticeSIMParams(params);

loaded = load(resultPath);
result = loaded.result;
diagnostics = result.diagnostics;
phaseSearch = diagnostics.phaseSearch;
carriers = diagnostics.latticeEstimation.diagnostics.phaseOnlyEstimate.carriers;

fprintf('INPUT=%s\n', inputPath);
fprintf('saved raw modulation S=%g T=%g; protected S=%g T=%g\n', ...
    diagnostics.rawEstimatedModulationS, diagnostics.rawEstimatedModulationT, ...
    diagnostics.protectedModulationS, diagnostics.protectedModulationT);
fprintf('saved phase S=%g T=%g; carriers ks=[%g,%g] kt=[%g,%g]\n', ...
    phaseSearch.phaseOffsetS, phaseSearch.phaseOffsetT, ...
    carriers.ksPixel(1), carriers.ksPixel(2), carriers.ktPixel(1), carriers.ktPixel(2));

[rawStack, ~] = readLatticeSIMInput(inputPath);
rawStack = double(rawStack);
fprintf('raw frame means:');
fprintf(' %.6g', squeeze(mean(mean(rawStack, 1), 2)));
fprintf('\n');
fprintf('raw temporal RMS/mean=%g\n', temporalRms(rawStack));

stackFFT = preprocessLikeReconstruction(rawStack, params);
otf = buildLatticeOTF(size(stackFFT, 1), size(stackFFT, 2), params);

axis = [0.02 0.05 0.08 0.10 0.15 0.20 0.35 0.50 0.70 1.00];
fprintf('\nDiagonal modulation score using current objective:\n');
fprintf('m\tfullScore\tmaskedScore\tSenergy/C0\tTenergy/C0\tSls\tTls\n');
for idx = 1:numel(axis)
    m = axis(idx);
    searchParams = params;
    searchParams.phaseOffsetS = phaseSearch.phaseOffsetS;
    searchParams.phaseOffsetT = phaseSearch.phaseOffsetT;
    searchParams.estimatedModulationS = m;
    searchParams.estimatedModulationT = m;
    bands = separateLatticeBandsFrequency(stackFFT, searchParams);
    fullScore = currentPhaseConsistencyScore(bands, carriers);
    maskedScore = maskedPhaseConsistencyScore(bands, carriers, otf.values);
    [sEnergy, tEnergy] = shiftedSidebandEnergyRatios(bands, carriers);
    [sLs, tLs] = shiftedLeastSquaresScale(bands, carriers, otf.values);
    fprintf('%.3f\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n', ...
        m, fullScore, maskedScore, sEnergy, tEnergy, sLs, tLs);
end

fprintf('\nOne-axis current objective, other axis fixed at 0.05:\n');
fprintf('m\tSscore(T=.05)\tTscore(S=.05)\n');
for idx = 1:numel(axis)
    m = axis(idx);
    sScore = oneAxisScore(stackFFT, params, phaseSearch, carriers, m, 0.05, 'S');
    tScore = oneAxisScore(stackFFT, params, phaseSearch, carriers, 0.05, m, 'T');
    fprintf('%.3f\t%.6g\t%.6g\n', m, sScore, tScore);
end

function stackFFT = preprocessLikeReconstruction(rawStack, params)
otf = buildLatticeOTF(size(rawStack, 1), size(rawStack, 2), params);
psf = abs(otf2psf(ifftshift(otf.values)));
windowedStack = importImages(rawStack);
deconvolvedStack = deconvlucy(windowedStack, psf, params.deconvolutionIterations);
stackFFT = zeros(size(rawStack));
for frameIdx = 1:size(rawStack, 3)
    stackFFT(:, :, frameIdx) = FFT2D(deconvolvedStack(:, :, frameIdx), false);
end
end

function value = temporalRms(rawStack)
frameMean = mean(rawStack, 3);
temporal = rawStack - frameMean;
value = sqrt(mean(temporal(:) .^ 2)) ./ mean(frameMean(:));
end

function score = oneAxisScore(stackFFT, params, phaseSearch, carriers, modulationS, modulationT, axisName)
searchParams = params;
searchParams.phaseOffsetS = phaseSearch.phaseOffsetS;
searchParams.phaseOffsetT = phaseSearch.phaseOffsetT;
searchParams.estimatedModulationS = modulationS;
searchParams.estimatedModulationT = modulationT;
bands = separateLatticeBandsFrequency(stackFFT, searchParams);
if axisName == 'S'
    score = sidebandPairScore(bands.C0, bands.CsPlus, bands.CsMinus, carriers.ksPixel);
else
    score = sidebandPairScore(bands.C0, bands.CtPlus, bands.CtMinus, carriers.ktPixel);
end
end

function score = currentPhaseConsistencyScore(bands, carriers)
score = sidebandPairScore(bands.C0, bands.CsPlus, bands.CsMinus, carriers.ksPixel) + ...
    sidebandPairScore(bands.C0, bands.CtPlus, bands.CtMinus, carriers.ktPixel);
end

function score = sidebandPairScore(centerBand, plusBand, minusBand, carrierPixel)
score = normalizedSidebandDifference(centerBand, plusBand, carrierPixel) + ...
    normalizedSidebandDifference(centerBand, minusBand, -carrierPixel);
end

function score = normalizedSidebandDifference(centerBand, sideBand, carrierPixel)
b0 = FFT2D(centerBand, true);
b1 = FFT2D(sideBand, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));
energy0 = sum(abs(b0) .^ 2, 'all');
score = -sum(abs(b1 - b0) .^ 2, 'all') ./ energy0;
end

function score = maskedPhaseConsistencyScore(bands, carriers, otfValues)
score = maskedSidebandPairScore(bands.C0, bands.CsPlus, bands.CsMinus, otfValues, carriers.ksPixel) + ...
    maskedSidebandPairScore(bands.C0, bands.CtPlus, bands.CtMinus, otfValues, carriers.ktPixel);
end

function score = maskedSidebandPairScore(centerBand, plusBand, minusBand, otfValues, carrierPixel)
score = maskedSidebandDifference(centerBand, plusBand, otfValues, carrierPixel) + ...
    maskedSidebandDifference(centerBand, minusBand, otfValues, -carrierPixel);
end

function score = maskedSidebandDifference(centerBand, sideBand, otfValues, carrierPixel)
b0 = FFT2D(centerBand, true);
b1 = FFT2D(sideBand, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));
support0 = abs(otfValues) > 1e-3;
support1 = abs(latticeFourierShift(otfValues, -carrierPixel(1), -carrierPixel(2))) > 1e-3;
mask = support0 & support1;
if ~any(mask(:))
    score = -inf;
    return;
end
energy0 = sum(abs(b0(mask)) .^ 2, 'all');
score = -sum(abs(b1(mask) - b0(mask)) .^ 2, 'all') ./ energy0;
end

function [sRatio, tRatio] = shiftedSidebandEnergyRatios(bands, carriers)
b0 = FFT2D(bands.C0, true);
sPlus = shiftToCenter(bands.CsPlus, carriers.ksPixel);
sMinus = shiftToCenter(bands.CsMinus, -carriers.ksPixel);
tPlus = shiftToCenter(bands.CtPlus, carriers.ktPixel);
tMinus = shiftToCenter(bands.CtMinus, -carriers.ktPixel);
cEnergy = sum(abs(b0) .^ 2, 'all');
sRatio = mean([sum(abs(sPlus) .^ 2, 'all'), sum(abs(sMinus) .^ 2, 'all')]) ./ cEnergy;
tRatio = mean([sum(abs(tPlus) .^ 2, 'all'), sum(abs(tMinus) .^ 2, 'all')]) ./ cEnergy;
end

function shifted = shiftToCenter(band, carrierPixel)
shifted = FFT2D(band, true);
shifted = latticeFourierShift(shifted, -carrierPixel(1), -carrierPixel(2));
end

function [sScale, tScale] = shiftedLeastSquaresScale(bands, carriers, otfValues)
sScale = mean([
    leastSquaresScale(bands.C0, bands.CsPlus, otfValues, carriers.ksPixel), ...
    leastSquaresScale(bands.C0, bands.CsMinus, otfValues, -carriers.ksPixel)]);
tScale = mean([
    leastSquaresScale(bands.C0, bands.CtPlus, otfValues, carriers.ktPixel), ...
    leastSquaresScale(bands.C0, bands.CtMinus, otfValues, -carriers.ktPixel)]);
end

function scale = leastSquaresScale(centerBand, sideBand, otfValues, carrierPixel)
b0 = FFT2D(centerBand, true);
b1 = FFT2D(sideBand, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));
support0 = abs(otfValues) > 1e-3;
support1 = abs(latticeFourierShift(otfValues, -carrierPixel(1), -carrierPixel(2))) > 1e-3;
mask = support0 & support1;
scale = real(sum(conj(b0(mask)) .* b1(mask), 'all') ./ sum(abs(b0(mask)) .^ 2, 'all'));
end
