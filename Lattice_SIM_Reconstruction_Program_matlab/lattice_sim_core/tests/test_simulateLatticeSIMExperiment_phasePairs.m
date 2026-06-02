function test_simulateLatticeSIMExperiment_phasePairs()
%TEST_SIMULATELATTICESIMEXPERIMENT_PHASEPAIRS Verify simulated acquisition phases are configurable.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

objectImage = makeLatticeSIMPhantom([32, 32]);
phasePairs = [
    0,       0
    pi / 4, 0
    pi / 2, 0
    0,       pi / 4
    pi / 4, pi / 2
];

simParams.imageSize = size(objectImage);
simParams.ksPixel = [5, 2];
simParams.ktPixel = [-2, 6];
simParams.phasePairs = phasePairs;
simParams.randomSeed = 11;
simParams.phaseErrorStd = 0;
simParams.useOTF = false;

[~, truth] = simulateLatticeSIMExperiment(objectImage, simParams);

assert(isequal(truth.phasePairs, phasePairs));
assert(isequal(truth.ksPixel, simParams.ksPixel));
assert(isequal(truth.ktPixel, simParams.ktPixel));
end
