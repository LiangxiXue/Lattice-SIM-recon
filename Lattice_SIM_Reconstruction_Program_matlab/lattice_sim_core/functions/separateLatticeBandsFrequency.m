function bands = separateLatticeBandsFrequency(stackFFT, params)
%SEPARATELATTICEBANDSFREQUENCY Demodulate five frequency-domain Lattice-SIM frames.

if nargin < 2
    params = defaultLatticeSIMParams();
end

validateLatticeSIMStack(stackFFT);

W = latticeFrequencyPhaseMatrix(params);
[h, w, ~] = size(stackFFT);
reshaped = reshape(double(stackFFT), h * w, 5);
components = reshaped / transpose(W);
components = reshape(components, h, w, 5);

bands.C0 = components(:, :, 1);
bands.CsPlus = components(:, :, 2);
bands.CsMinus = components(:, :, 3);
bands.CtPlus = components(:, :, 4);
bands.CtMinus = components(:, :, 5);
bands.domain = 'frequency';
bands.phaseMatrix = W;
bands.diagnostics.functionName = 'separateLatticeBandsFrequency';
bands.diagnostics.inputDomain = 'frequency';
bands.diagnostics.outputDomain = 'frequency';
bands.diagnostics.model = 'lattice-phase-matrix';
bands.diagnostics.phaseOffsetS = params.phaseOffsetS;
bands.diagnostics.phaseOffsetT = params.phaseOffsetT;
bands.diagnostics.modulationS = params.estimatedModulationS;
bands.diagnostics.modulationT = params.estimatedModulationT;
end

function W = latticeFrequencyPhaseMatrix(params)
phasePairs = latticePhasePairs();
W = zeros(5, 5);
for idx = 1:5
    phiS = phasePairs(idx, 1) + params.phaseOffsetS;
    phiT = phasePairs(idx, 2) + params.phaseOffsetT;
    W(idx, :) = [1, ...
        params.estimatedModulationS .* exp(1i * phiS), ...
        params.estimatedModulationS .* exp(-1i * phiS), ...
        params.estimatedModulationT .* exp(1i * phiT), ...
        params.estimatedModulationT .* exp(-1i * phiT)];
end
end
