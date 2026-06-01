function test_simulateLatticeSIMExperiment()
%TEST_SIMULATELATTICESIMEXPERIMENT Verify OpenSIM-style Lattice forward simulation.

objectImage = makeLatticeSIMPhantom([48, 52]);

simParams.imageSize = [48, 52];
simParams.ksPixel = [7, 0];
simParams.ktPixel = [0, 8];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 65;
simParams.emissionWavelengthNm = 532;
simParams.NA = 1.2;
simParams.noiseLevel = 0.03;
simParams.phaseErrorStd = 0.02;
simParams.randomSeed = 11;
simParams.useOTF = true;

[rawStack, truth] = simulateLatticeSIMExperiment(objectImage, simParams);

assert(isequal(size(rawStack), [48, 52, 5]));
assert(isequal(size(truth.object), [48, 52]));
assert(isequal(size(truth.illumination), [48, 52, 5]));
assert(isequal(size(truth.noiseFreeStack), [48, 52, 5]));
assert(isequal(size(truth.otf), [48, 52]));
assert(isequal(truth.phasePairs, [
    0,       0
    0,       2*pi/3
    0,       4*pi/3
    2*pi/3, 0
    4*pi/3, 2*pi/3
]));
assert(norm(truth.phaseErrors(:)) > 0);
assert(norm(rawStack(:) - truth.noiseFreeStack(:)) > 0);
assert(all(isfinite(rawStack(:))));
assert(min(rawStack(:)) >= 0);

[rawStackAgain, truthAgain] = simulateLatticeSIMExperiment(objectImage, simParams);
assert(max(abs(rawStack(:) - rawStackAgain(:))) < 1e-12);
assert(max(abs(truth.phaseErrors(:) - truthAgain.phaseErrors(:))) < 1e-12);

fprintf('test_simulateLatticeSIMExperiment passed.\n');
end
