function test_estimateLatticeBandParameters()
%TEST_ESTIMATELATTICEBANDPARAMETERS Verify full Lattice carrier and scale estimation.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;
params.carrierMinRadiusPixels = 2;
params.carrierWeakPeakRatio = 1.2;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCarrierRefinementIterations = 3;

image = makeLatticeSIMPhantom([64, 64]);
ks = [8, 0];
kt = [0, 10];
phiSOffset = 0.4;
phiTOffset = -0.6;
mS = 0.55;
mT = 0.35;

bands.C0 = FFT2D(image, false);
bands.CsPlus = FFT2D(mS * exp(1i * phiSOffset) .* ...
    latticeFourierShift(image, ks(1), ks(2)), false);
bands.CsMinus = FFT2D(mS * exp(-1i * phiSOffset) .* ...
    latticeFourierShift(image, -ks(1), -ks(2)), false);
bands.CtPlus = FFT2D(mT * exp(1i * phiTOffset) .* ...
    latticeFourierShift(image, kt(1), kt(2)), false);
bands.CtMinus = FFT2D(mT * exp(-1i * phiTOffset) .* ...
    latticeFourierShift(image, -kt(1), -kt(2)), false);
bands.domain = 'frequency';

estimate = estimateLatticeBandParameters(bands, params);

assert(norm(estimate.carriers.ksPixel - ks) < 0.75);
assert(norm(estimate.carriers.ktPixel - kt) < 0.75);
assert(abs(estimate.phaseOffsetS + phiSOffset) < 0.15);
assert(abs(estimate.phaseOffsetT + phiTOffset) < 0.15);
assert(abs(estimate.modulationS - mS) < 0.20);
assert(abs(estimate.modulationT - mT) < 0.20);
assert(isfield(estimate.diagnostics, 'warnings'));

params.pixelSizeNm = 19.5;
ks = [20, 0];
kt = [0, 22];
bands.CsPlus = FFT2D(mS * exp(1i * phiSOffset) .* ...
    latticeFourierShift(image, ks(1), ks(2)), false);
bands.CsMinus = FFT2D(mS * exp(-1i * phiSOffset) .* ...
    latticeFourierShift(image, -ks(1), -ks(2)), false);
bands.CtPlus = FFT2D(mT * exp(1i * phiTOffset) .* ...
    latticeFourierShift(image, kt(1), kt(2)), false);
bands.CtMinus = FFT2D(mT * exp(-1i * phiTOffset) .* ...
    latticeFourierShift(image, -kt(1), -kt(2)), false);

estimate = estimateLatticeBandParameters(bands, params);

assert(norm(estimate.carriers.ksPixel - ks) < 0.75);
assert(norm(estimate.carriers.ktPixel - kt) < 0.75);
assert(estimate.modulationS > 0);
assert(estimate.modulationS <= 1);
assert(estimate.modulationT > 0);
assert(estimate.modulationT <= 1);
assert(any(contains(estimate.diagnostics.warnings, 'overlap')));

params.pixelSizeNm = 97.5;
params.carrierWeakPeakRatio = 1.05;
ksPlus = [8, 0];
ksMinus = [-12, 0];
bands.CsPlus = FFT2D(0.15 * mS * exp(1i * phiSOffset) .* ...
    latticeFourierShift(image, ksPlus(1), ksPlus(2)), false);
bands.CsMinus = FFT2D(mS * exp(-1i * phiSOffset) .* ...
    latticeFourierShift(image, ksMinus(1), ksMinus(2)), false);
bands.CtPlus = FFT2D(mT * exp(1i * phiTOffset) .* ...
    latticeFourierShift(image, kt(1), kt(2)), false);
bands.CtMinus = FFT2D(mT * exp(-1i * phiTOffset) .* ...
    latticeFourierShift(image, -kt(1), -kt(2)), false);

estimate = estimateLatticeBandParameters(bands, params);

assert(estimate.carriers.ksPixel(1) > 10.5);
assert(estimate.carriers.ksPixel(1) < 12.75);
assert(isfield(estimate.diagnostics, 'sidebandMeasurements'));
assert(isfield(estimate.diagnostics, 'pairResidualS'));
assert(estimate.diagnostics.pairResidualS > 1);
assert(any(contains(estimate.diagnostics.warnings, 'S sideband pair is inconsistent')));
end
