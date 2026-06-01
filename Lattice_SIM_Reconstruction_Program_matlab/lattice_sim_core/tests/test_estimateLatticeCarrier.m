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

axisBands = makeAxisCarrierBands(imageSize);
params.carrierSearchMode = "axis-aligned";
axisCarriers = estimateLatticeCarrier(axisBands, params);
assert(norm(axisCarriers.ksPixel - [8, 0]) <= 1.25);
assert(norm(axisCarriers.ktPixel - [0, 10]) <= 1.25);

flatBands = bands;
flatBands.CsPlus = ones(imageSize);
flatBands.CsMinus = ones(imageSize);
flatBands.CtPlus = ones(imageSize);
flatBands.CtMinus = ones(imageSize);
expectError(@() estimateLatticeCarrier(flatBands, params), 'LatticeSIM:CarrierEstimationFailed');

fprintf('test_estimateLatticeCarrier passed.\n');
end

function bands = makeAxisCarrierBands(imageSize)
[h, w] = deal(imageSize(1), imageSize(2));
[x, y] = meshgrid(0:w-1, 0:h-1);
axisS = exp(1i * 2*pi * 8 * x / w);
axisT = exp(1i * 2*pi * 10 * y / h);
diagonalS = exp(1i * 2*pi * (8 * x / w + 8 * y / h));
diagonalT = exp(1i * 2*pi * (8 * x / w + 10 * y / h));

bands.C0 = ones(imageSize);
bands.CsPlus = axisS + 2 * diagonalS;
bands.CsMinus = conj(bands.CsPlus);
bands.CtPlus = axisT + 2 * diagonalT;
bands.CtMinus = conj(bands.CtPlus);
end
