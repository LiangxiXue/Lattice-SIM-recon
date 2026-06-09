function test_configurable_real_data_main()
%TEST_CONFIGURABLE_REAL_DATA_MAIN Verify the configurable real-data entrypoint.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(coreDir);
addpath(fullfile(coreDir, 'functions'));

scriptPath = fullfile(coreDir, 'main', 'main.m');
sourceText = fileread(scriptPath);
sourceLines = regexp(strtrim(sourceText), '\r?\n', 'split');
assert(~startsWith(strtrim(sourceLines{1}), 'function'), ...
    'Expected main.m to be a script, not a function.');

assert(~contains(sourceText, 'dryRun'), ...
    'The real-data entrypoint should always execute reconstruction and should not expose dryRun.');
assert(~contains(sourceText, 'clear all'), ...
    'The real-data entrypoint should not clear preset workspace configuration variables.');
assert(~contains(sourceText, 'modulationS ='), ...
    'The real-data entrypoint should not expose a manual modulationS input block.');
assert(~contains(sourceText, 'modulationT ='), ...
    'The real-data entrypoint should not expose a manual modulationT input block.');
assert(contains(sourceText, 'params.enableLatticeParameterEstimation = true'), ...
    'Lattice parameter estimation must stay enabled for this entrypoint.');
assert(contains(sourceText, 'params.carrierSearchMode = "unconstrained"'), ...
    'Carrier search must stay unconstrained for this entrypoint.');
assert(contains(sourceText, 'params.fusionMode = "hifi-two-step"'), ...
    'The real-data entrypoint should select the HiFi-style two-step fusion mode.');
assert(contains(sourceText, 'params.wienerW1 = 0.15'), ...
    'The real-data entrypoint should expose the W1 Wiener constant.');
assert(contains(sourceText, 'params.wienerW2 = 0.06'), ...
    'The real-data entrypoint should expose the W2 Wiener constant.');
assert(contains(sourceText, 'params.hifiDenominatorScaleW1 = 1.2'), ...
    'The real-data entrypoint should expose the W1 denominator scale.');
assert(contains(sourceText, 'params.hifiDenominatorScaleW2 = 0.8'), ...
    'The real-data entrypoint should expose the W2 denominator scale.');
assert(contains(sourceText, 'uigetfile'), ...
    'The real-data entrypoint should use uigetfile to select a five-frame stack.');
assert(contains(sourceText, 'exist(stackPath, ''file'') ~= 2'), ...
    'The real-data entrypoint should re-prompt when a preset stackPath points to a missing file.');
assert(~contains(sourceText, 'framePaths'), ...
    'The real-data entrypoint should not configure five separate frame paths.');
assert(~contains(sourceText, 'cropEnabled'), ...
    'The real-data entrypoint should not crop a complete five-frame stack.');

showFigures = false;
dataDir = fullfile(fileparts(coreDir), 'Lattice_SIM_test_Data');
singleFramePaths = {
    fullfile(dataDir, '1.tiff')
    fullfile(dataDir, '2.tiff')
    fullfile(dataDir, '3.tiff')
    fullfile(dataDir, '4.tiff')
    fullfile(dataDir, '5.tiff')
};
[stack, ~] = readCroppedLatticeSIMFrames(singleFramePaths, [1024, 1024], []);
stackPath = fullfile(tempdir, ['lattice_sim_stack_', char(java.util.UUID.randomUUID()), '.tif']);
for idx = 1:size(stack, 3)
    if idx == 1
        imwrite(uint8(stack(:, :, idx)), stackPath, 'tif');
    else
        imwrite(uint8(stack(:, :, idx)), stackPath, 'tif', 'WriteMode', 'append');
    end
end
outputDir = fullfile(tempdir, ['lattice_sim_output_', char(java.util.UUID.randomUUID())]);
paramOverrides = struct( ...
    'enableLatticeParameterEstimation', false, ...
    'carrierSearchMode', "axis-aligned");
run(scriptPath);

assert(exist('config', 'var') == 1);
assert(exist('result', 'var') == 1);
assert(exist('saved', 'var') == 1);

assert(strcmp(config.input.mode, 'tiffStack'));
assert(strcmp(config.input.stackPath, stackPath));
assert(exist(config.input.stackPath, 'file') == 2);

assert(config.params.pixelSizeNm > 0);
assert(config.params.emissionWavelengthNm > 0);
assert(config.params.NA > 0);
assert(isfield(config.params, 'phasePairs'));
assert(isequal(size(config.params.phasePairs), [5, 2]));
assert(isfield(config.params, 'phaseMatrix'));
assert(isequal(size(config.params.phaseMatrix), [5, 5]));
assert(strcmp(char(config.params.preprocessingMode), 'hifi-rl-fft'));
assert(strcmp(char(config.params.carrierSearchMode), 'unconstrained'));
assert(config.params.enableLatticeParameterEstimation == true);
assert(strcmp(char(config.params.fusionMode), 'hifi-two-step'));
assert(config.params.wienerW1 == 0.15);
assert(config.params.wienerW2 == 0.06);
assert(config.params.hifiDenominatorScaleW1 == 1.2);
assert(config.params.hifiDenominatorScaleW2 == 0.8);
assert(config.params.normalizeFrames == true);
assert(strcmp(config.outputDir, outputDir));
end
