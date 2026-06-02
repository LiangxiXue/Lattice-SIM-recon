function test_custom_lattice_phase_matrix()
%TEST_CUSTOM_LATTICE_PHASE_MATRIX Verify configurable acquisition phase model.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

stackFFT = zeros(8, 8, 5);
params = defaultLatticeSIMParams();
params.phasePairs = [
    0,       0
    pi / 2, 0
    pi,      0
    0,       pi / 2
    pi / 2, pi
];

bands = separateLatticeBandsFrequency(stackFFT, params);

expected = expectedPhaseMatrix(params.phasePairs, params);
assert(max(abs(bands.phaseMatrix(:) - expected(:))) < 1e-12);

params = rmfield(params, 'phasePairs');
params.phaseMatrix = expected;
bands = separateLatticeBandsFrequency(stackFFT, params);
assert(max(abs(bands.phaseMatrix(:) - expected(:))) < 1e-12);
end

function W = expectedPhaseMatrix(phasePairs, params)
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
