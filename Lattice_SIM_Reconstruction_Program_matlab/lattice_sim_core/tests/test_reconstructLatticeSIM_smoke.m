function test_reconstructLatticeSIM_smoke()
%TEST_RECONSTRUCTLATTICESIM_SMOKE Verify complete public API output.

imageSize = [32, 36];
[rawStack, ~] = makeSyntheticLatticeSIMStack(imageSize, ...
    'ksPixel', [6, 0], 'ktPixel', [0, 7], 'carrierAmplitude', 0.2);

params = defaultLatticeSIMParams();
params.pixelSizeNm = 65;
params.emissionWavelengthNm = 532;
params.NA = 1.2;
params.normalizeFrames = false;

result = reconstructLatticeSIM(rawStack, params);

assert(isfield(result, 'WF'));
assert(isfield(result, 'SIM'));
assert(isfield(result, 'params'));
assert(isfield(result, 'diagnostics'));
assert(isequal(size(result.WF), imageSize * 2));
assert(isequal(size(result.SIM), imageSize * 2));
assert(all(isfinite(result.WF(:))));
assert(all(isfinite(result.SIM(:))));
expectedBands = separateLatticeBands(rawStack, params);
expectedWF = imresize(abs(expectedBands.C0), 2);
assert(max(abs(result.WF(:) - expectedWF(:))) < 1e-10);
assert(~isfield(result.diagnostics, 'centerBandWF'));
assert(isfield(result.diagnostics, 'carrierS'));
assert(isfield(result.diagnostics, 'carrierT'));
assert(nnz(result.diagnostics.combine.supportMask) > nnz(result.diagnostics.otf.supportMask));

missingParams = rmfield(params, 'NA');
expectError(@() reconstructLatticeSIM(rawStack, missingParams), 'LatticeSIM:MissingPhysicalParameter');

objectImage = makeLatticeSIMPhantom(imageSize);
simParams.imageSize = imageSize;
simParams.ksPixel = [6, 0];
simParams.ktPixel = [0, 7];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = params.pixelSizeNm;
simParams.emissionWavelengthNm = params.emissionWavelengthNm;
simParams.NA = params.NA;
simParams.noiseLevel = 0.01;
simParams.phaseErrorStd = 0.0;
simParams.randomSeed = 3;
simParams.useOTF = true;
[experimentStack, experimentTruth] = simulateLatticeSIMExperiment(objectImage, simParams);

experimentResult = reconstructLatticeSIM(experimentStack, params);
assert(isequal(size(experimentResult.WF), imageSize * 2));
assert(isequal(size(experimentResult.SIM), imageSize * 2));
assert(all(isfinite(experimentResult.SIM(:))));
assert(isequal(experimentTruth.ksPixel, [6, 0]));

fprintf('test_reconstructLatticeSIM_smoke passed.\n');
end
