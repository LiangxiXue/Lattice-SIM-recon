function bands = separateLatticeBandsFrequency(stackFFT, params)
%SEPARATELATTICEBANDSFREQUENCY Demodulate five frequency-domain Lattice-SIM frames.

if nargin < 2
    params = defaultLatticeSIMParams();
end

validateLatticeSIMStack(stackFFT);

W = latticePhaseMatrix(params);
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
