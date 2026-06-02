function test_lattice_estimation_params()
%TEST_LATTICE_ESTIMATION_PARAMS Verify defaults for two-pass Lattice estimation.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();

assert(strcmp(char(params.separationFunction), 'separateLatticeBandsFrequency'));
assert(strcmp(char(params.separationInputDomain), 'frequency'));
assert(params.phaseOffsetS == 0);
assert(params.phaseOffsetT == 0);
assert(params.estimatedModulationS == 1);
assert(params.estimatedModulationT == 1);
assert(params.enableLatticeParameterEstimation == true);
assert(params.latticeCarrierRefinementIterations == 3);
assert(params.latticeCarrierSearchStepPixels == 2.5);
assert(params.latticeCorrelationOverlap == 0.15);
assert(strcmp(char(params.carrierSearchMode), 'unconstrained'));
end
