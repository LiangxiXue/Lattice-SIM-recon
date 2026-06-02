function test_separateLatticeBandsFrequency()
%TEST_SEPARATELATTICEBANDSFREQUENCY Verify Lattice frequency-domain demodulation.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

imageSize = [16, 16];
[x, y] = meshgrid(0:imageSize(2)-1, 0:imageSize(1)-1);
C0 = 10 + zeros(imageSize);
CsPlus = 2 * exp(1i * 2*pi*(3*x/imageSize(2)));
CsMinus = conj(CsPlus);
CtPlus = 1.5 * exp(1i * 2*pi*(4*y/imageSize(1)));
CtMinus = conj(CtPlus);

expected.C0 = fftshift(fft2(C0));
expected.CsPlus = fftshift(fft2(CsPlus));
expected.CsMinus = fftshift(fft2(CsMinus));
expected.CtPlus = fftshift(fft2(CtPlus));
expected.CtMinus = fftshift(fft2(CtMinus));

phasePairs = [
    0,       0
    2*pi/3, 0
    4*pi/3, 0
    0,       2*pi/3
    2*pi/3, 4*pi/3
];
stackFFT = zeros([imageSize, 5]);
for idx = 1:5
    phiS = phasePairs(idx, 1);
    phiT = phasePairs(idx, 2);
    stackFFT(:, :, idx) = expected.C0 ...
        + expected.CsPlus .* exp(1i * phiS) ...
        + expected.CsMinus .* exp(-1i * phiS) ...
        + expected.CtPlus .* exp(1i * phiT) ...
        + expected.CtMinus .* exp(-1i * phiT);
end

params = defaultLatticeSIMParams();
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;

bands = separateLatticeBandsFrequency(stackFFT, params);

assert(max(abs(bands.C0(:) - expected.C0(:))) < 1e-8);
assert(max(abs(bands.CsPlus(:) - expected.CsPlus(:))) < 1e-8);
assert(max(abs(bands.CsMinus(:) - expected.CsMinus(:))) < 1e-8);
assert(max(abs(bands.CtPlus(:) - expected.CtPlus(:))) < 1e-8);
assert(max(abs(bands.CtMinus(:) - expected.CtMinus(:))) < 1e-8);
assert(strcmp(bands.domain, 'frequency'));
assert(strcmp(bands.diagnostics.functionName, 'separateLatticeBandsFrequency'));
assert(strcmp(bands.diagnostics.model, 'lattice-phase-matrix'));
end
