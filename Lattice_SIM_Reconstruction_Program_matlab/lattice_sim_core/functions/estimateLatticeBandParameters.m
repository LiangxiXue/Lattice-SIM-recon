function estimate = estimateLatticeBandParameters(bands, params)
%ESTIMATELATTICEBANDPARAMETERS Estimate Lattice carriers, phase, and modulation.

if nargin < 2
    params = defaultLatticeSIMParams();
end

[coarseCarriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(bands.C0, 1), size(bands.C0, 2), params);

peakS = latticeFitPeak(bands.C0, bands.CsPlus, otf.values, ...
    coarseCarriers.ksPixel, params.latticeCorrelationOverlap, ...
    params.latticeCarrierSearchStepPixels, params.latticeCarrierRefinementIterations);
peakT = latticeFitPeak(bands.C0, bands.CtPlus, otf.values, ...
    coarseCarriers.ktPixel, params.latticeCorrelationOverlap, ...
    params.latticeCarrierSearchStepPixels, params.latticeCarrierRefinementIterations);

refinedCarriers = coarseCarriers;
refinedCarriers.ksPixel = [peakS.kx, peakS.ky];
refinedCarriers.ktPixel = [peakT.kx, peakT.ky];
[h, w] = size(bands.C0);
refinedCarriers.ksRadPerPixel = [2*pi*refinedCarriers.ksPixel(1)/w, 2*pi*refinedCarriers.ksPixel(2)/h];
refinedCarriers.ktRadPerPixel = [2*pi*refinedCarriers.ktPixel(1)/w, 2*pi*refinedCarriers.ktPixel(2)/h];

correlationS = latticeGetPeak(bands.C0, bands.CsPlus, otf.values, ...
    refinedCarriers.ksPixel, params.latticeCorrelationOverlap);
correlationT = latticeGetPeak(bands.C0, bands.CtPlus, otf.values, ...
    refinedCarriers.ktPixel, params.latticeCorrelationOverlap);

estimate.carriers = refinedCarriers;
estimate.coarseCarriers = coarseCarriers;
estimate.phaseOffsetS = -angle(correlationS);
estimate.phaseOffsetT = -angle(correlationT);
estimate.modulationS = min(max(abs(correlationS), eps), 1);
estimate.modulationT = min(max(abs(correlationT), eps), 1);
estimate.correlationS = correlationS;
estimate.correlationT = correlationT;
estimate.fitPeakS = peakS;
estimate.fitPeakT = peakT;
estimate.diagnostics = carrierDiagnostics;
estimate.diagnostics.refinedCarrierS = refinedCarriers.ksPixel;
estimate.diagnostics.refinedCarrierT = refinedCarriers.ktPixel;
estimate.diagnostics.phaseOffsetS = estimate.phaseOffsetS;
estimate.diagnostics.phaseOffsetT = estimate.phaseOffsetT;
estimate.diagnostics.modulationS = estimate.modulationS;
estimate.diagnostics.modulationT = estimate.modulationT;
estimate.diagnostics.correlationMagnitudeS = abs(correlationS);
estimate.diagnostics.correlationMagnitudeT = abs(correlationT);
estimate.diagnostics.fitPeakS = peakS;
estimate.diagnostics.fitPeakT = peakT;
end
