function test_lattice_modulation_compensation()
%TEST_LATTICE_MODULATION_COMPENSATION Guard against double sideband scaling.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

objectImage = makeLatticeSIMPhantom([48, 48]);
simParams.imageSize = size(objectImage);
simParams.ksPixel = [6, 0];
simParams.ktPixel = [0, 7];
simParams.modulationS = 0.08;
simParams.modulationT = 0.07;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 97.5;
simParams.emissionWavelengthNm = 561;
simParams.NA = 1.42;
simParams.noiseLevel = 0.0;
simParams.phaseErrorStd = 0.0;
simParams.randomSeed = 11;
simParams.useOTF = true;

rawStack = simulateLatticeSIMExperiment(objectImage, simParams);

params = defaultLatticeSIMParams();
params.pixelSizeNm = simParams.pixelSizeNm;
params.emissionWavelengthNm = simParams.emissionWavelengthNm;
params.NA = simParams.NA;
params.modulationS = 0.7;
params.modulationT = 0.7;
params.modulationProtectionEnabled = true;
params.modulationMinReliable = 0.35;
params.modulationFallback = 0.7;
params.enableLatticeParameterEstimation = true;

result = reconstructLatticeSIM(rawStack, params);

assert(result.diagnostics.latticeEstimation.modulationS < 0.5);
assert(result.diagnostics.latticeEstimation.modulationT < 0.5);
assert(result.diagnostics.latticeEstimation.rawEstimatedModulationS < params.modulationMinReliable);
assert(result.diagnostics.latticeEstimation.rawEstimatedModulationT < params.modulationMinReliable);
assert(result.diagnostics.latticeEstimation.protectedModulationS == params.modulationFallback);
assert(result.diagnostics.latticeEstimation.protectedModulationT == params.modulationFallback);
assert(result.diagnostics.separation.modulationS == params.modulationFallback);
assert(result.diagnostics.separation.modulationT == params.modulationFallback);
assert(result.params.modulationS == params.modulationFallback);
assert(result.params.modulationT == params.modulationFallback);
assert(strcmp(result.diagnostics.combine.modulationCompensationMode, 'separation-matrix'));
end
