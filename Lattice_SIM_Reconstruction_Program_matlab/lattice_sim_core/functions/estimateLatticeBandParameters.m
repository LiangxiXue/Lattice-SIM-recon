function estimate = estimateLatticeBandParameters(bands, params)
%ESTIMATELATTICEBANDPARAMETERS Estimate Lattice carriers, phase, and modulation.

if nargin < 2
    params = defaultLatticeSIMParams();
end

[coarseCarriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(bands.C0, 1), size(bands.C0, 2), params);

try
    measurementS = measureSidebandPair(bands.C0, bands.CsPlus, bands.CsMinus, ...
        otf.values, coarseCarriers.ksPixel, params, 'S');
    measurementT = measureSidebandPair(bands.C0, bands.CtPlus, bands.CtMinus, ...
        otf.values, coarseCarriers.ktPixel, params, 'T');

    refinedCarriers = coarseCarriers;
    refinedCarriers.ksPixel = measurementS.carrierPixel;
    refinedCarriers.ktPixel = measurementT.carrierPixel;
    [h, w] = size(bands.C0);
    refinedCarriers.ksRadPerPixel = [2*pi*refinedCarriers.ksPixel(1)/w, 2*pi*refinedCarriers.ksPixel(2)/h];
    refinedCarriers.ktRadPerPixel = [2*pi*refinedCarriers.ktPixel(1)/w, 2*pi*refinedCarriers.ktPixel(2)/h];
    refinedCarriers.sidebandCarrierPixels.ksPlus = measurementS.plus.carrierPixel;
    refinedCarriers.sidebandCarrierPixels.ksMinus = measurementS.minus.carrierPixel;
    refinedCarriers.sidebandCarrierPixels.ktPlus = measurementT.plus.carrierPixel;
    refinedCarriers.sidebandCarrierPixels.ktMinus = measurementT.minus.carrierPixel;

    correlationS = measurementS.projectedCorrelation;
    correlationT = measurementT.projectedCorrelation;
    modulationS = measurementS.projectedModulation;
    modulationT = measurementT.projectedModulation;
    fitPeakS = makeProjectedPeak(measurementS);
    fitPeakT = makeProjectedPeak(measurementT);
catch err
    if ~strcmp(err.identifier, 'LatticeSIM:DegenerateCorrelationReference')
        rethrow(err);
    end

    refinedCarriers = coarseCarriers;
    peakS = makeFallbackPeak(refinedCarriers.ksPixel);
    peakT = makeFallbackPeak(refinedCarriers.ktPixel);
    correlationS = estimateFallbackCorrelation(bands.C0, bands.CsPlus, bands.CsMinus);
    correlationT = estimateFallbackCorrelation(bands.C0, bands.CtPlus, bands.CtMinus);
    modulationS = min(max(abs(correlationS), eps), 1);
    modulationT = min(max(abs(correlationT), eps), 1);
    fitPeakS = peakS;
    fitPeakT = peakT;
    carrierDiagnostics = appendDiagnosticWarning(carrierDiagnostics, ...
        'No reliable C0/sideband OTF overlap; using coarse carriers and spectrum-energy modulation estimates.');
end

estimate.carriers = refinedCarriers;
estimate.coarseCarriers = coarseCarriers;
estimate.phaseOffsetS = -angle(correlationS);
estimate.phaseOffsetT = -angle(correlationT);
estimate.modulationS = modulationS;
estimate.modulationT = modulationT;
estimate.correlationS = correlationS;
estimate.correlationT = correlationT;
estimate.fitPeakS = fitPeakS;
estimate.fitPeakT = fitPeakT;
estimate.diagnostics = carrierDiagnostics;
estimate.diagnostics.refinedCarrierS = refinedCarriers.ksPixel;
estimate.diagnostics.refinedCarrierT = refinedCarriers.ktPixel;
estimate.diagnostics.phaseOffsetS = estimate.phaseOffsetS;
estimate.diagnostics.phaseOffsetT = estimate.phaseOffsetT;
estimate.diagnostics.modulationS = estimate.modulationS;
estimate.diagnostics.modulationT = estimate.modulationT;
estimate.diagnostics.correlationMagnitudeS = abs(correlationS);
estimate.diagnostics.correlationMagnitudeT = abs(correlationT);
estimate.diagnostics.fitPeakS = estimate.fitPeakS;
estimate.diagnostics.fitPeakT = estimate.fitPeakT;
if exist('measurementS', 'var') && exist('measurementT', 'var')
    estimate.diagnostics.sidebandMeasurements.S = measurementS;
    estimate.diagnostics.sidebandMeasurements.T = measurementT;
    estimate.diagnostics.pairResidualS = measurementS.pairResidual;
    estimate.diagnostics.pairResidualT = measurementT.pairResidual;
    estimate.diagnostics.sidebandProjectionWeightsS = measurementS.weights;
    estimate.diagnostics.sidebandProjectionWeightsT = measurementT.weights;
    if measurementS.pairResidual > params.sidebandPairResidualWarningPixels
        estimate.diagnostics = appendDiagnosticWarning(estimate.diagnostics, ...
            sprintf('S sideband pair is inconsistent before projection (residual %.3f pixels).', measurementS.pairResidual));
    end
    if measurementT.pairResidual > params.sidebandPairResidualWarningPixels
        estimate.diagnostics = appendDiagnosticWarning(estimate.diagnostics, ...
            sprintf('T sideband pair is inconsistent before projection (residual %.3f pixels).', measurementT.pairResidual));
    end
end
end

function measurement = measureSidebandPair(centerBand, plusBand, minusBand, otfValues, coarseCarrier, params, label)
peakPlus = latticeFitPeak(centerBand, plusBand, otfValues, ...
    coarseCarrier, params.latticeCorrelationOverlap, ...
    params.latticeCarrierSearchStepPixels, params.latticeCarrierRefinementIterations);
peakMinus = latticeFitPeak(centerBand, minusBand, otfValues, ...
    -coarseCarrier, params.latticeCorrelationOverlap, ...
    params.latticeCarrierSearchStepPixels, params.latticeCarrierRefinementIterations);

plusCarrier = [peakPlus.kx, peakPlus.ky];
minusCarrier = [peakMinus.kx, peakMinus.ky];
correlationPlus = latticeGetPeak(centerBand, plusBand, otfValues, ...
    plusCarrier, params.latticeCorrelationOverlap);
correlationMinus = latticeGetPeak(centerBand, minusBand, otfValues, ...
    minusCarrier, params.latticeCorrelationOverlap);

weights = normalizePairWeights([abs(correlationPlus), abs(correlationMinus)]);
positiveCarrierEvidence = [plusCarrier; -minusCarrier];
carrierPixel = sum(positiveCarrierEvidence .* weights(:), 1);
positiveCorrelationEvidence = [correlationPlus, conj(correlationMinus)];
projectedCorrelation = sum(positiveCorrelationEvidence .* weights);
projectedModulation = sum(abs(positiveCorrelationEvidence) .* weights);

measurement.label = label;
measurement.plus = makeSidebandMeasurement(plusCarrier, peakPlus, correlationPlus);
measurement.minus = makeSidebandMeasurement(minusCarrier, peakMinus, correlationMinus);
measurement.weights = weights;
measurement.positiveCarrierEvidence = positiveCarrierEvidence;
measurement.carrierPixel = carrierPixel;
measurement.pairResidual = norm(positiveCarrierEvidence(1, :) - positiveCarrierEvidence(2, :));
measurement.projectedCorrelation = projectedCorrelation;
measurement.projectedModulation = min(max(projectedModulation, eps), 1);
end

function sideband = makeSidebandMeasurement(carrierPixel, peak, correlation)
sideband.carrierPixel = carrierPixel;
sideband.peak = peak;
sideband.correlation = correlation;
sideband.correlationMagnitude = abs(correlation);
sideband.correlationPhase = angle(correlation);
end

function weights = normalizePairWeights(rawWeights)
weights = rawWeights;
weights(~isfinite(weights)) = 0;
if sum(weights) <= eps
    weights = [0.5, 0.5];
else
    weights = weights ./ sum(weights);
end
end

function peak = makeProjectedPeak(measurement)
peak.kx = measurement.carrierPixel(1);
peak.ky = measurement.carrierPixel(2);
peak.resPhase = angle(measurement.projectedCorrelation);
peak.resMag = abs(measurement.projectedCorrelation);
peak.correlation = measurement.projectedCorrelation;
peak.control.plus = measurement.plus.peak.control;
peak.control.minus = measurement.minus.peak.control;
end

function peak = makeFallbackPeak(carrierPixel)
peak.kx = carrierPixel(1);
peak.ky = carrierPixel(2);
peak.resPhase = 0;
peak.resMag = 0;
peak.correlation = 0;
peak.control = zeros(10, 10, 1);
end

function correlation = estimateFallbackCorrelation(centerBand, plusBand, minusBand)
centerScale = rmsMagnitude(centerBand);
sideScale = mean([rmsMagnitude(plusBand), rmsMagnitude(minusBand)]);
if centerScale <= eps || sideScale <= eps
    magnitude = eps;
else
    magnitude = min(max(sideScale / centerScale, eps), 1);
end
correlation = complex(magnitude, 0);
end

function value = rmsMagnitude(image)
value = sqrt(mean(abs(image(:)) .^ 2));
end
