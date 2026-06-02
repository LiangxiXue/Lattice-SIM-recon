function test_lattice_hifi_style_peak_helpers()
%TEST_LATTICE_HIFI_STYLE_PEAK_HELPERS Verify HiFi-style Lattice peak helpers.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

test_latticeFourierShift_matches_hifi_phase_ramp();
test_latticeGetPeak_recovers_complex_scale();
test_latticeFitPeak_refines_carrier_from_correlation();
end

function test_latticeFourierShift_matches_hifi_phase_ramp()
image = ones(8, 8);
shifted = latticeFourierShift(image, 2, -1);
[x, y] = meshgrid(0:7, 0:7);
expected = exp(2*pi*1i*((-1) * y / 8 + 2 * x / 8));

assert(max(abs(shifted(:) - expected(:))) < 1e-12);
end

function test_latticeGetPeak_recovers_complex_scale()
image = makeLatticeSIMPhantom([32, 32]);
band0 = FFT2D(image, false);
trueScale = 0.42 * exp(1i * 0.7);
band1 = trueScale .* band0;
otfValues = ones(size(image));

peak = latticeGetPeak(band0, band1, otfValues, [0, 0], 0.15);

assert(abs(abs(peak) - 0.42) < 1e-3);
assert(abs(angle(peak) - 0.7) < 1e-3);
end

function test_latticeFitPeak_refines_carrier_from_correlation()
image = makeLatticeSIMPhantom([32, 32]);
trueCarrier = [4, -3];
trueScale = 0.55 * exp(1i * -0.35);
band0 = FFT2D(image, false);
band1Space = trueScale .* latticeFourierShift(image, trueCarrier(1), trueCarrier(2));
band1 = FFT2D(band1Space, false);
otfValues = ones(size(image));

peak = latticeFitPeak(band0, band1, otfValues, [3.2, -2.4], 0.15, 2.5, 3);

assert(norm([peak.kx, peak.ky] - trueCarrier) < 0.35);
assert(abs(abs(peak.correlation) - 0.55) < 0.05);
assert(abs(angle(peak.correlation) - angle(trueScale)) < 0.15);
end
