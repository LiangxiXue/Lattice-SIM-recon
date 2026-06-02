function test_hifi_style_frequency_pipeline()
%TEST_HIFI_STYLE_FREQUENCY_PIPELINE Verify reconstructLatticeSIM follows HiFi-SIM preprocessing.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

objectImage = makeLatticeSIMPhantom([32, 32]);

simParams.imageSize = size(objectImage);
simParams.ksPixel = [4, 0];
simParams.ktPixel = [0, 5];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 97.5;
simParams.emissionWavelengthNm = 561;
simParams.NA = 1.42;
simParams.noiseLevel = 0.0;
simParams.phaseErrorStd = 0.0;
simParams.randomSeed = 7;
simParams.useOTF = true;

rawStack = simulateLatticeSIMExperiment(objectImage, simParams);

params = defaultLatticeSIMParams();
params.pixelSizeNm = simParams.pixelSizeNm;
params.emissionWavelengthNm = simParams.emissionWavelengthNm;
params.NA = simParams.NA;
params.modulationS = simParams.modulationS;
params.modulationT = simParams.modulationT;

result = reconstructLatticeSIM(rawStack, params);

assert(strcmp(result.diagnostics.preprocessing.method, 'hifi-rl-fft'), ...
    'Expected reconstructLatticeSIM to use HiFi-SIM RL + FFT preprocessing.');
assert(result.diagnostics.preprocessing.deconvolutionIterations == 5, ...
    'Expected five Richardson-Lucy iterations.');
assert(strcmp(result.diagnostics.separation.functionName, 'separateLatticeBandsFrequency'), ...
    'Expected frequency-domain separation to use separateLatticeBandsFrequency.');
assert(strcmp(result.diagnostics.separation.inputDomain, 'frequency'), ...
    'Expected separation input to be frequency-domain frames.');
assert(strcmp(result.diagnostics.separation.model, 'lattice-phase-matrix'), ...
    'Expected Lattice phase matrix, not HiFi harmonic comp matrix.');
assert(isfield(result.diagnostics, 'latticeEstimation'), ...
    'Expected first-pass Lattice parameter estimation diagnostics.');
assert(isfield(result.diagnostics.latticeEstimation, 'phaseOffsetS'));
assert(isfield(result.diagnostics.latticeEstimation, 'phaseOffsetT'));
assert(isfield(result.diagnostics.latticeEstimation, 'modulationS'));
assert(isfield(result.diagnostics.latticeEstimation, 'modulationT'));
assert(strcmp(result.diagnostics.combine.bandDomain, 'frequency'), ...
    'Expected combineLatticeSpectrum to consume frequency-domain bands directly.');
assert(abs(result.diagnostics.combine.outputOtfPixelSizeNm - params.pixelSizeNm / 2) < eps, ...
    'Expected 2x output OTF to use half-size pixels.');
assert(isSpatiallyAligned(result.SIM, result.WF), ...
    'Expected SIM output to stay aligned with the C0 wide-field image.');
assertSidebandOtfSupportsMatchShiftedOrders(result.diagnostics.combine, ...
    result.diagnostics.latticeEstimation.carriers);
assert(all(isfinite(result.SIM(:))), 'SIM output must be finite.');
end

function assertSidebandOtfSupportsMatchShiftedOrders(combineDiagnostics, carriers)
expectedCenters = [
    -carriers.ksPixel
     carriers.ksPixel
    -carriers.ktPixel
     carriers.ktPixel
];
for idx = 1:4
    actualCenter = maskCentroidPixels(combineDiagnostics.bandTaperMasks{idx + 1});
    expectedCenter = expectedCenters(idx, :);
    assert(norm(actualCenter - expectedCenter) < 2, ...
        'Expected shifted sideband OTF support center to match the shifted order.');
end
end

function center = maskCentroidPixels(mask)
[h, w] = size(mask);
[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
weights = double(mask);
weightSum = sum(weights(:));
center = [sum(x(:) .* weights(:)) / weightSum, ...
    sum(y(:) .* weights(:)) / weightSum];
end

function aligned = isSpatiallyAligned(imageA, imageB)
imageA = normalizeForCorrelation(imageA);
imageB = normalizeForCorrelation(imageB);
correlation = real(ifft2(fft2(imageA) .* conj(fft2(imageB))));
[~, maxIdx] = max(correlation(:));
[row, col] = ind2sub(size(correlation), maxIdx);
rowShift = min(row - 1, size(correlation, 1) - row + 1);
colShift = min(col - 1, size(correlation, 2) - col + 1);
aligned = rowShift <= 1 && colShift <= 1;
end

function image = normalizeForCorrelation(image)
image = double(image);
image = image - mean(image(:));
scale = std(image(:));
if scale > 0
    image = image ./ scale;
end
end
