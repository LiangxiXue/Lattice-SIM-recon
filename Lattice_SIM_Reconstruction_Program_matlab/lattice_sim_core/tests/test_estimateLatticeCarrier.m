function test_estimateLatticeCarrier()
%TEST_ESTIMATELATTICECARRIER Verify Lattice carrier estimation after demodulation.

imageSize = [64, 64];
[rawStack, ~] = makeSyntheticLatticeSIMStack(imageSize, ...
    'ksPixel', [8, 0], 'ktPixel', [0, 10], 'carrierAmplitude', 0.25);
bands = separateLatticeBands(rawStack);

params = defaultLatticeSIMParams();
params.carrierMinRadiusPixels = 3;
params.carrierPeakWindow = 1;

[carriers, diagnostics] = estimateLatticeCarrier(bands, params);

assert(norm(carriers.ksPixel - [8, 0]) <= 1.25);
assert(norm(carriers.ktPixel - [0, 10]) <= 1.25);
assert(abs(diagnostics.orthogonalityErrorDeg) < 5);
assert(isfield(diagnostics, 'carrierSearchMaps'));

flatBands = bands;
flatBands.CsPlus = ones(imageSize);
flatBands.CsMinus = ones(imageSize);
flatBands.CtPlus = ones(imageSize);
flatBands.CtMinus = ones(imageSize);
expectError(@() estimateLatticeCarrier(flatBands, params), 'LatticeSIM:CarrierEstimationFailed');

fprintf('test_estimateLatticeCarrier passed.\n');
end
