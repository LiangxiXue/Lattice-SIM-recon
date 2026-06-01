%RUN_LATTICE_SIM_RECONSTRUCTION_EXAMPLE Minimal function-based Lattice-SIM example.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'functions'));

objectImage = makeLatticeSIMPhantom([128, 128]);

simParams.imageSize = size(objectImage);
simParams.ksPixel = [14, 0];
simParams.ktPixel = [0, 16];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 65;
simParams.emissionWavelengthNm = 532;
simParams.NA = 1.2;
simParams.noiseLevel = 0.02;
simParams.phaseErrorStd = 0.01;
simParams.randomSeed = 7;
simParams.useOTF = true;

[syntheticInput, truth] = simulateLatticeSIMExperiment(objectImage, simParams);

params = defaultLatticeSIMParams();
params.pixelSizeNm = 65;
params.emissionWavelengthNm = 532;
params.NA = 1.2;
params.normalizeFrames = false;

result = reconstructLatticeSIM(syntheticInput, params);

outputDir = fullfile(rootDir, 'tests', 'results', 'example_output');
saveLatticeSIMResult(result, outputDir);
save(fullfile(outputDir, 'diagnostics', 'simulation_truth.mat'), 'truth');

fprintf('Saved example outputs to %s\n', outputDir);
