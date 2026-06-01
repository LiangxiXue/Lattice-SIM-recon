function test_combineLatticeSpectrum_modulation()
%TEST_COMBINELATTICESPECTRUM_MODULATION Verify sideband amplitude compensation.

imageSize = [128, 128];
objectImage = makeLatticeSIMPhantom(imageSize);

simParams.imageSize = imageSize;
simParams.ksPixel = [30, 0];
simParams.ktPixel = [0, 30];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 97.5;
simParams.emissionWavelengthNm = 561;
simParams.NA = 1.42;
simParams.noiseLevel = 0;
simParams.phaseErrorStd = 0;
simParams.randomSeed = 11;
simParams.useOTF = true;

[rawStack, truth] = simulateLatticeSIMExperiment(objectImage, simParams);

params = defaultLatticeSIMParams();
params.pixelSizeNm = simParams.pixelSizeNm;
params.emissionWavelengthNm = simParams.emissionWavelengthNm;
params.NA = simParams.NA;
params.modulationS = simParams.modulationS;
params.modulationT = simParams.modulationT;
params.normalizeFrames = false;

bands = separateLatticeBands(rawStack, params);
carriers = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(imageSize(1), imageSize(2), params);
[~, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params);

relativeError = rawSpectrumRelativeError(diagnostics.simSpectrum, centeredFft2(truth.object), ...
    diagnostics.wienerDenominator);
assert(relativeError < 0.05, ...
    'Expected raw combined spectrum relative error < 0.05, got %.6g.', relativeError);

fprintf('test_combineLatticeSpectrum_modulation passed.\n');
end

function errorValue = rawSpectrumRelativeError(actualSpectrum, expectedSpectrum, denominator)
supportMask = denominator > 0;
expectedMagnitude = abs(expectedSpectrum);
supportMask = supportMask & expectedMagnitude > max(expectedMagnitude(:)) * 1e-8;

expected = expectedSpectrum(supportMask);
actual = actualSpectrum(supportMask);
diff = actual - expected;
errorValue = norm(diff(:)) / max(norm(expected(:)), eps);
end

function spectrum = centeredFft2(image)
spectrum = fftshift(fft2(ifftshift(image)));
end
