function W = latticePhaseMatrix(params)
%LATTICEPHASEMATRIX Return the five-frame Lattice-SIM demodulation matrix.

if nargin < 1
    params = struct();
end

if isfield(params, 'phaseMatrix') && ~isempty(params.phaseMatrix)
    W = params.phaseMatrix;
    return;
end

if isfield(params, 'phasePairs') && ~isempty(params.phasePairs)
    phasePairs = params.phasePairs;
else
    phasePairs = latticePhasePairs();
end

phaseOffsetS = getParam(params, 'phaseOffsetS', 0);
phaseOffsetT = getParam(params, 'phaseOffsetT', 0);
modulationS = getParam(params, 'estimatedModulationS', 1);
modulationT = getParam(params, 'estimatedModulationT', 1);

W = zeros(5, 5);
for idx = 1:5
    phiS = phasePairs(idx, 1) + phaseOffsetS;
    phiT = phasePairs(idx, 2) + phaseOffsetT;
    W(idx, :) = [1, modulationS .* exp(1i * phiS), ...
        modulationS .* exp(-1i * phiS), ...
        modulationT .* exp(1i * phiT), ...
        modulationT .* exp(-1i * phiT)];
end
end

function value = getParam(params, name, defaultValue)
if isfield(params, name)
    value = params.(name);
else
    value = defaultValue;
end
end
