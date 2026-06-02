function test_configurable_real_data_main()
%TEST_CONFIGURABLE_REAL_DATA_MAIN Verify the configurable real-data entrypoint.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(coreDir);
addpath(fullfile(coreDir, 'functions'));

scriptPath = fullfile(coreDir, 'main', 'run_lattice_sim_configurable_test.m');
sourceText = fileread(scriptPath);
sourceLines = regexp(strtrim(sourceText), '\r?\n', 'split');
assert(~startsWith(strtrim(sourceLines{1}), 'function'), ...
    'Expected run_lattice_sim_configurable_test.m to be a script, not a function.');

dryRun = true;
showFigures = false;
run(scriptPath);

assert(exist('config', 'var') == 1);
assert(exist('result', 'var') == 0);
assert(exist('saved', 'var') == 0);

assert(strcmp(config.input.mode, 'realFiles'));
assert(numel(config.input.framePaths) == 5);
assert(config.input.cropEnabled == true);
assert(isequal(config.input.cropSizePixels, [1024, 1024]));
assert(isempty(config.input.cropCenterPixels));
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
assert(config.params.enableLatticeParameterEstimation == false);
assert(contains(config.outputDir, 'configurable_real_data_output'));

clear config output result saved;
dryRun = true;
showFigures = false;
paramOverrides = struct('wiener', 0.08);
run(scriptPath);

assert(config.params.wiener == 0.08);
assert(config.params.normalizeFrames == true);
assert(strcmp(char(config.params.preprocessingMode), 'hifi-rl-fft'));
end
