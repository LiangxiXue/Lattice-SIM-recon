function W = makeLatticePhaseMatrix(params)
%MAKELATTICEPHASEMATRIX Build the five-frame Lattice-SIM demodulation matrix.

if nargin < 1
    params = defaultLatticeSIMParams();
else
    params = defaultLatticeSIMParams(params);
end

if ~isempty(params.phaseMatrix)
    W = params.phaseMatrix;
    return;
end

phasePairs = params.phasePairs;
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
