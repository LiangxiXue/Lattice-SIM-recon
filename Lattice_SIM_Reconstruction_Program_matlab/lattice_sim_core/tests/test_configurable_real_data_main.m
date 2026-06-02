function test_configurable_real_data_main()
%TEST_CONFIGURABLE_REAL_DATA_MAIN Verify the configurable real-data entrypoint.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(coreDir);
addpath(fullfile(coreDir, 'functions'));

config = run_lattice_sim_configurable_test('dryRun', true, 'showFigures', false);

assert(strcmp(config.input.mode, 'realFiles'));
assert(numel(config.input.framePaths) == 5);
for idx = 1:5
    assert(exist(config.input.framePaths{idx}, 'file') == 2);
end

assert(config.params.pixelSizeNm > 0);
assert(config.params.emissionWavelengthNm > 0);
assert(config.params.NA > 0);
assert(isfield(config.params, 'phasePairs'));
assert(isequal(size(config.params.phasePairs), [5, 2]));
assert(isfield(config.params, 'phaseMatrix'));
assert(isequal(size(config.params.phaseMatrix), [5, 5]));
assert(strcmp(char(config.params.preprocessingMode), 'hifi-rl-fft'));
assert(strcmp(char(config.params.carrierSearchMode), 'unconstrained'));
assert(contains(config.outputDir, 'configurable_real_data_output'));
end
