scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
coreDir = fileparts(fileparts(scriptDir));
functionsDir = fullfile(coreDir, 'functions');
addpath(genpath(functionsDir));

rng(7);

inputPath = fullfile(coreDir, 'roi.tif');
resultPath = fullfile(coreDir, 'output', 'diagnostics', 'result.mat');

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 610;
params.NA = 1.42;
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;
params.carrierSearchMode = "unconstrained";
params.latticeCorrelationOverlap = 0.15;

defaultPhasePairs = [
    0,       0
    2*pi/3, 0
    4*pi/3, 0
    0,       2*pi/3
    2*pi/3, 4*pi/3
];

loaded = load(resultPath);
diagnostics = loaded.result.diagnostics;
carriers = diagnostics.latticeEstimation.diagnostics.phaseOnlyEstimate.carriers;

[rawStack, ~] = readLatticeSIMInput(inputPath);
rawStack = double(rawStack);
otf = buildLatticeOTF(size(rawStack, 1), size(rawStack, 2), params);
psf = abs(otf2psf(ifftshift(otf.values)));
deconvolvedStack = deconvlucy(importImages(rawStack), psf, 5);
stackFFT = zeros(size(rawStack));
for frameIdx = 1:size(rawStack, 3)
    stackFFT(:, :, frameIdx) = FFT2D(deconvolvedStack(:, :, frameIdx), false);
end

baselineParams = params;
baselineParams.phasePairs = defaultPhasePairs;
baselineParams.phaseOffsetS = diagnostics.phaseOffsetS;
baselineParams.phaseOffsetT = diagnostics.phaseOffsetT;
[baselineFastS, baselineFastT] = fastEffectiveModulation(stackFFT, baselineParams, carriers, otf.values);
[baselineSlowS, baselineSlowT] = slowEffectiveModulation(stackFFT, baselineParams, carriers, otf.values);

fprintf('Baseline default phasePairs + saved offsets:\n');
fprintf('  fast S=%g T=%g mean=%g\n', baselineFastS, baselineFastT, mean([baselineFastS, baselineFastT]));
fprintf('  slow S=%g T=%g mean=%g\n', baselineSlowS, baselineSlowT, mean([baselineSlowS, baselineSlowT]));
fprintf('  saved phaseOnly modulation S=%g T=%g\n', ...
    diagnostics.phaseOnlyEstimate.modulationS, diagnostics.phaseOnlyEstimate.modulationT);

numLocal = 300;
numRandom = 300;
topCount = 12;
condLimit = 20;
rows = zeros(numLocal + numRandom + 1, 15);
labels = strings(numLocal + numRandom + 1, 1);

rowIdx = 1;
[rows(rowIdx, :), labels(rowIdx)] = scoreCandidate("baseline", defaultPhasePairs, ...
    diagnostics.phaseOffsetS, diagnostics.phaseOffsetT, stackFFT, params, carriers, otf.values);

for idx = 1:numLocal
    jitter = randn(size(defaultPhasePairs)) * 0.45;
    jitter(1, :) = 0;
    phasePairs = wrapPhase(defaultPhasePairs + jitter);
    rowIdx = rowIdx + 1;
    [rows(rowIdx, :), labels(rowIdx)] = scoreCandidate("local", phasePairs, 0, 0, ...
        stackFFT, params, carriers, otf.values);
end

for idx = 1:numRandom
    phasePairs = wrapPhase((rand(5, 2) * 2 - 1) * pi);
    phasePairs = wrapPhase(phasePairs - phasePairs(1, :));
    rowIdx = rowIdx + 1;
    [rows(rowIdx, :), labels(rowIdx)] = scoreCandidate("random", phasePairs, 0, 0, ...
        stackFFT, params, carriers, otf.values);
end

rows = rows(1:rowIdx, :);
labels = labels(1:rowIdx);
valid = isfinite(rows(:, 1)) & rows(:, 4) <= condLimit;
validRows = rows(valid, :);
validLabels = labels(valid);
[validRows, sortIdx] = sortrows(validRows, -1);
validLabels = validLabels(sortIdx);

fprintf('\nFast-screened candidates: total=%d valid(cond<=%g)=%d\n', ...
    rowIdx, condLimit, size(validRows, 1));
printGroupSummary("local", rows, labels);
printGroupSummary("random", rows, labels);
printStrictCondSummary(5, rows, labels);
printStrictCondSummary(10, rows, labels);
fprintf('Top candidates by fast mean effective modulation:\n');
fprintf('rank\tlabel\tfastMean\tfastS\tfastT\tcond\tphasePairs(row-major)\n');
for idx = 1:min(topCount, size(validRows, 1))
    fprintf('%d\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t', ...
        idx, validLabels(idx), validRows(idx, 1), validRows(idx, 2), validRows(idx, 3), validRows(idx, 4));
    fprintf('%.4f ', validRows(idx, 6:15));
    fprintf('\n');
end

fprintf('\nSlow OTF-overlap validation for top candidates:\n');
fprintf('rank\tlabel\tfastMean\tslowMean\tslowS\tslowT\tcond\n');
for idx = 1:min(topCount, size(validRows, 1))
    phasePairs = reshape(validRows(idx, 6:15), 5, 2);
    candidateParams = params;
    candidateParams.phasePairs = phasePairs;
    candidateParams.phaseOffsetS = 0;
    candidateParams.phaseOffsetT = 0;
    [slowS, slowT] = slowEffectiveModulation(stackFFT, candidateParams, carriers, otf.values);
    fprintf('%d\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n', ...
        idx, validLabels(idx), validRows(idx, 1), mean([slowS, slowT]), slowS, slowT, validRows(idx, 4));
    fprintf('  phasePairs = [\n');
    for pairIdx = 1:size(phasePairs, 1)
        fprintf('    %.6f %.6f\n', phasePairs(pairIdx, 1), phasePairs(pairIdx, 2));
    end
    fprintf('  ];\n');
end

function [row, label] = scoreCandidate(label, phasePairs, phaseOffsetS, phaseOffsetT, ...
        stackFFT, params, carriers, otfValues)
candidateParams = params;
candidateParams.phasePairs = phasePairs;
candidateParams.phaseOffsetS = phaseOffsetS;
candidateParams.phaseOffsetT = phaseOffsetT;
candidateParams.estimatedModulationS = 1;
candidateParams.estimatedModulationT = 1;
W = makeLatticePhaseMatrix(candidateParams);
matrixCond = cond(W);
try
    [fastS, fastT] = fastEffectiveModulation(stackFFT, candidateParams, carriers, otfValues);
    fastMean = mean([fastS, fastT]);
catch
    fastS = NaN;
    fastT = NaN;
    fastMean = -inf;
end
row = [fastMean, fastS, fastT, matrixCond, rand(), reshape(phasePairs, 1, [])];
end

function printGroupSummary(groupName, rows, labels)
mask = labels == groupName & isfinite(rows(:, 1)) & rows(:, 4) <= 20;
if ~any(mask)
    fprintf('Best %s candidate: none\n', groupName);
    return;
end
groupRows = rows(mask, :);
groupRows = sortrows(groupRows, -1);
fprintf('Best %s fast candidate: mean=%g S=%g T=%g cond=%g\n', ...
    groupName, groupRows(1, 1), groupRows(1, 2), groupRows(1, 3), groupRows(1, 4));
end

function printStrictCondSummary(condLimit, rows, labels)
mask = labels == "random" & isfinite(rows(:, 1)) & rows(:, 4) <= condLimit;
if ~any(mask)
    fprintf('Best random cond<=%g fast candidate: none\n', condLimit);
    return;
end
groupRows = rows(mask, :);
groupRows = sortrows(groupRows, -1);
fprintf('Best random cond<=%g fast candidate: mean=%g S=%g T=%g cond=%g\n', ...
    condLimit, groupRows(1, 1), groupRows(1, 2), groupRows(1, 3), groupRows(1, 4));
end

function [modS, modT] = fastEffectiveModulation(stackFFT, params, carriers, otfValues)
bands = separateLatticeBandsFrequency(stackFFT, params);
modS = mean([
    maskedLeastSquaresScale(bands.C0, bands.CsPlus, otfValues, carriers.ksPixel), ...
    maskedLeastSquaresScale(bands.C0, bands.CsMinus, otfValues, -carriers.ksPixel)]);
modT = mean([
    maskedLeastSquaresScale(bands.C0, bands.CtPlus, otfValues, carriers.ktPixel), ...
    maskedLeastSquaresScale(bands.C0, bands.CtMinus, otfValues, -carriers.ktPixel)]);
modS = abs(modS);
modT = abs(modT);
end

function [modS, modT] = slowEffectiveModulation(stackFFT, params, carriers, otfValues)
bands = separateLatticeBandsFrequency(stackFFT, params);
modS = mean([
    abs(latticeGetPeak(bands.C0, bands.CsPlus, otfValues, carriers.ksPixel, 0.15)), ...
    abs(latticeGetPeak(bands.C0, bands.CsMinus, otfValues, -carriers.ksPixel, 0.15))]);
modT = mean([
    abs(latticeGetPeak(bands.C0, bands.CtPlus, otfValues, carriers.ktPixel, 0.15)), ...
    abs(latticeGetPeak(bands.C0, bands.CtMinus, otfValues, -carriers.ktPixel, 0.15))]);
end

function scale = maskedLeastSquaresScale(centerBand, sideBand, otfValues, carrierPixel)
b0 = FFT2D(centerBand, true);
b1 = FFT2D(sideBand, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));
support0 = abs(otfValues) > 1e-3;
support1 = abs(latticeFourierShift(otfValues, -carrierPixel(1), -carrierPixel(2))) > 1e-3;
mask = support0 & support1;
denominator = sum(abs(b0(mask)) .^ 2, 'all');
if denominator <= eps
    scale = NaN;
else
    scale = real(sum(conj(b0(mask)) .* b1(mask), 'all') ./ denominator);
end
end

function value = wrapPhase(value)
value = angle(exp(1i * value));
end
