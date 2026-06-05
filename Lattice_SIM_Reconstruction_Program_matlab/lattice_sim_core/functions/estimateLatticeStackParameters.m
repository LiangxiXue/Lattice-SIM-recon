function estimate = estimateLatticeStackParameters(stackFFT, params)
%ESTIMATELATTICESTACKPARAMETERS Estimate Lattice parameters from raw frequency frames.

if nargin < 2
    params = defaultLatticeSIMParams();
else
    params = defaultLatticeSIMParams(params);
end

validateLatticeSIMStack(stackFFT);

carrierParams = params;
carrierParams.phaseOffsetS = 0;
carrierParams.phaseOffsetT = 0;
carrierParams.estimatedModulationS = 1;
carrierParams.estimatedModulationT = 1;
firstPassBands = separateLatticeBandsFrequency(stackFFT, carrierParams);
[coarseCarriers, carrierDiagnostics] = estimateLatticeCarrier(firstPassBands, params);

otf = buildLatticeOTF(size(stackFFT, 1), size(stackFFT, 2), params);
phaseSearch = searchLatticePhaseOffsets(stackFFT, coarseCarriers, otf.values, params);

phaseParams = params;
phaseParams.phaseOffsetS = phaseSearch.phaseOffsetS;
phaseParams.phaseOffsetT = phaseSearch.phaseOffsetT;
phaseParams.estimatedModulationS = 1;
phaseParams.estimatedModulationT = 1;
phaseOnlyBands = separateLatticeBandsFrequency(stackFFT, phaseParams);
phaseOnlyEstimate = estimateLatticeBandParameters(phaseOnlyBands, params);

phaseSearch.phaseOffsetS = wrapPhase(phaseSearch.phaseOffsetS + phaseOnlyEstimate.phaseOffsetS);
phaseSearch.phaseOffsetT = wrapPhase(phaseSearch.phaseOffsetT + phaseOnlyEstimate.phaseOffsetT);
phaseSearch.residualPhaseOffsetS = phaseOnlyEstimate.phaseOffsetS;
phaseSearch.residualPhaseOffsetT = phaseOnlyEstimate.phaseOffsetT;
phaseParams.phaseOffsetS = phaseSearch.phaseOffsetS;
phaseParams.phaseOffsetT = phaseSearch.phaseOffsetT;

modulationSearch = searchLatticeModulations(stackFFT, phaseOnlyEstimate.carriers, otf.values, ...
    phaseSearch, params);

phaseParams.estimatedModulationS = modulationSearch.modulationS;
phaseParams.estimatedModulationT = modulationSearch.modulationT;
phaseBands = separateLatticeBandsFrequency(stackFFT, phaseParams);

residualEstimate = estimateLatticeBandParameters(phaseBands, params);

estimate = residualEstimate;
estimate.coarseCarriers = coarseCarriers;
estimate.phaseSearch = phaseSearch;
estimate.modulationSearch = modulationSearch;
estimate.phaseOffsetS = wrapPhase(phaseSearch.phaseOffsetS + residualEstimate.phaseOffsetS);
estimate.phaseOffsetT = wrapPhase(phaseSearch.phaseOffsetT + residualEstimate.phaseOffsetT);
estimate.modulationS = modulationSearch.modulationS;
estimate.modulationT = modulationSearch.modulationT;
estimate.diagnostics = residualEstimate.diagnostics;
estimate.diagnostics.estimationMode = 'stack-phase-search';
estimate.diagnostics.coarseCarrierDiagnostics = carrierDiagnostics;
estimate.diagnostics.phaseOnlyEstimate = phaseOnlyEstimate;
estimate.diagnostics.phaseSearch = phaseSearch;
estimate.diagnostics.modulationSearch = modulationSearch;
estimate.diagnostics.phaseOffsetS = estimate.phaseOffsetS;
estimate.diagnostics.phaseOffsetT = estimate.phaseOffsetT;
estimate.diagnostics.modulationS = estimate.modulationS;
estimate.diagnostics.modulationT = estimate.modulationT;
end

function phaseSearch = searchLatticePhaseOffsets(stackFFT, carriers, otfValues, params)
coarseStep = params.latticePhaseSearchCoarseStepRad;
fineStep = params.latticePhaseSearchFineStepRad;
fineRadius = params.latticePhaseSearchFineRadiusRad;

coarseAxis = -pi:coarseStep:(pi - coarseStep);
[bestS, bestT, bestScore] = searchPhaseGrid(stackFFT, carriers, otfValues, ...
    params, coarseAxis, coarseAxis);

fineAxisS = bestS - fineRadius:fineStep:bestS + fineRadius;
fineAxisT = bestT - fineRadius:fineStep:bestT + fineRadius;
[bestS, bestT, bestScore] = searchPhaseGrid(stackFFT, carriers, otfValues, ...
    params, fineAxisS, fineAxisT);

phaseSearch.phaseOffsetS = wrapPhase(bestS);
phaseSearch.phaseOffsetT = wrapPhase(bestT);
phaseSearch.score = bestScore;
phaseSearch.coarseStepRad = coarseStep;
phaseSearch.fineStepRad = fineStep;
phaseSearch.fineRadiusRad = fineRadius;
end

function modulationSearch = searchLatticeModulations(stackFFT, carriers, otfValues, phaseSearch, params)
minModulation = params.latticeModulationSearchMin;
maxModulation = params.latticeModulationSearchMax;
coarseStep = params.latticeModulationSearchCoarseStep;
fineStep = params.latticeModulationSearchFineStep;
fineRadius = params.latticeModulationSearchFineRadius;

coarseAxis = minModulation:coarseStep:maxModulation;
[bestS, bestT, bestScore] = searchModulationGrid(stackFFT, carriers, otfValues, ...
    phaseSearch, params, coarseAxis, coarseAxis);

fineAxisS = boundedAxis(bestS - fineRadius, bestS + fineRadius, fineStep, ...
    minModulation, maxModulation);
fineAxisT = boundedAxis(bestT - fineRadius, bestT + fineRadius, fineStep, ...
    minModulation, maxModulation);
[bestS, bestT, bestScore] = searchModulationGrid(stackFFT, carriers, otfValues, ...
    phaseSearch, params, fineAxisS, fineAxisT);

modulationSearch.modulationS = bestS;
modulationSearch.modulationT = bestT;
modulationSearch.score = bestScore;
modulationSearch.coarseStep = coarseStep;
modulationSearch.fineStep = fineStep;
modulationSearch.fineRadius = fineRadius;
end

function [bestS, bestT, bestScore] = searchModulationGrid(stackFFT, carriers, otfValues, ...
        phaseSearch, params, axisS, axisT)
bestS = axisS(1);
bestT = axisT(1);
bestScore = -inf;

searchParams = params;
searchParams.phaseOffsetS = phaseSearch.phaseOffsetS;
searchParams.phaseOffsetT = phaseSearch.phaseOffsetT;

for sIdx = 1:numel(axisS)
    for tIdx = 1:numel(axisT)
        searchParams.estimatedModulationS = axisS(sIdx);
        searchParams.estimatedModulationT = axisT(tIdx);
        bands = separateLatticeBandsFrequency(stackFFT, searchParams);
        score = phaseConsistencyScore(bands, carriers, otfValues, params);
        if score > bestScore
            bestScore = score;
            bestS = searchParams.estimatedModulationS;
            bestT = searchParams.estimatedModulationT;
        end
    end
end
end

function [bestS, bestT, bestScore] = searchPhaseGrid(stackFFT, carriers, otfValues, params, axisS, axisT)
bestS = 0;
bestT = 0;
bestScore = -inf;

searchParams = params;
searchParams.estimatedModulationS = 1;
searchParams.estimatedModulationT = 1;

for sIdx = 1:numel(axisS)
    for tIdx = 1:numel(axisT)
        searchParams.phaseOffsetS = wrapPhase(axisS(sIdx));
        searchParams.phaseOffsetT = wrapPhase(axisT(tIdx));
        bands = separateLatticeBandsFrequency(stackFFT, searchParams);
        score = phaseConsistencyScore(bands, carriers, otfValues, params);
        if score > bestScore
            bestScore = score;
            bestS = searchParams.phaseOffsetS;
            bestT = searchParams.phaseOffsetT;
        end
    end
end
end

function score = phaseConsistencyScore(bands, carriers, otfValues, params)
scoreS = sidebandPairScore(bands.C0, bands.CsPlus, bands.CsMinus, ...
    otfValues, carriers.ksPixel, params);
scoreT = sidebandPairScore(bands.C0, bands.CtPlus, bands.CtMinus, ...
    otfValues, carriers.ktPixel, params);
score = scoreS + scoreT;
end

function score = sidebandPairScore(centerBand, plusBand, minusBand, otfValues, carrierPixel, params)
plusScore = normalizedSidebandCorrelation(centerBand, plusBand, otfValues, ...
    carrierPixel, params);
minusScore = normalizedSidebandCorrelation(centerBand, minusBand, otfValues, ...
    -carrierPixel, params);
score = plusScore + minusScore;
end

function score = normalizedSidebandCorrelation(centerBand, sideBand, ~, carrierPixel, ~)
b0 = FFT2D(centerBand, true);
b1 = FFT2D(sideBand, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));

energy0 = sum(abs(b0).^2, 'all');
if energy0 <= eps
    score = -inf;
else
    score = -sum(abs(b1 - b0).^2, 'all') ./ energy0;
end
end

function value = wrapPhase(value)
value = angle(exp(1i * value));
end

function axis = boundedAxis(firstValue, lastValue, step, lowerBound, upperBound)
firstValue = max(firstValue, lowerBound);
lastValue = min(lastValue, upperBound);
axis = firstValue:step:lastValue;
if isempty(axis) || axis(end) < lastValue
    axis = [axis, lastValue];
end
axis = unique(axis);
end
