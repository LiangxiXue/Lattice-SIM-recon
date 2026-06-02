function bands = separateLatticeBands(stack, params)
%SEPARATELATTICEBANDS Demodulate five Lattice-SIM frames into five bands.

if nargin < 2
    params = defaultLatticeSIMParams();
end

validateLatticeSIMStack(stack);

W = latticePhaseMatrix(params);
[h, w, ~] = size(stack);
reshaped = reshape(double(stack), h * w, 5);
components = reshaped / transpose(W);
components = reshape(components, h, w, 5);

bands.C0 = components(:, :, 1);
bands.CsPlus = components(:, :, 2);
bands.CsMinus = components(:, :, 3);
bands.CtPlus = components(:, :, 4);
bands.CtMinus = components(:, :, 5);
bands.phaseMatrix = W;
if isfield(params, 'normalizationInfo')
    bands.normalization = params.normalizationInfo;
end
end
