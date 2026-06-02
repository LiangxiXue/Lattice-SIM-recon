%COMPARE_LATTICE_DIRECT_AND_CURRENT_SPECTRUM Save two Lattice-SIM spectrum comparison images.

scriptDir = fileparts(mfilename('fullpath'));
experimentsDir = fileparts(scriptDir);
coreDir = fileparts(experimentsDir);
functionsDir = fullfile(coreDir, 'functions');
addpath(functionsDir);

objectImage = makeLatticeSIMPhantom([128, 128]);

simParams.imageSize = size(objectImage);
simParams.ksPixel = [14, 0];
simParams.ktPixel = [0, 16];
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 65;
simParams.emissionWavelengthNm = 532;
simParams.NA = 1.2;
simParams.noiseLevel = 0.02;
simParams.phaseErrorStd = 0.01;
simParams.randomSeed = 7;
simParams.useOTF = true;

rawStack = simulateLatticeSIMExperiment(objectImage, simParams);

params = defaultLatticeSIMParams();
params.pixelSizeNm = 65;
params.emissionWavelengthNm = 532;
params.NA = 1.2;
params.normalizeFrames = false;

[stack, normalizationInfo] = normalizeSIMFrames(rawStack, params);
params.normalizationInfo = normalizationInfo;
bands = separateLatticeBands(stack, params);
[carriers, ~] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(stack, 1) * 2, size(stack, 2) * 2, params);

directCombinedSpectrum = buildDirectCombinedSpectrum(bands, carriers, otf, params);
result = reconstructLatticeSIM(rawStack, params);
currentCombinedSpectrum = result.diagnostics.combine.simSpectrum;

imwrite(spectrumForDisplay(directCombinedSpectrum), ...
    fullfile(scriptDir, 'lattice_direct_combined_spectrum.png'));
imwrite(spectrumForDisplay(currentCombinedSpectrum), ...
    fullfile(scriptDir, 'lattice_current_combined_spectrum.png'));

fprintf('Saved %s\n', fullfile(scriptDir, 'lattice_direct_combined_spectrum.png'));
fprintf('Saved %s\n', fullfile(scriptDir, 'lattice_current_combined_spectrum.png'));

function directCombinedSpectrum = buildDirectCombinedSpectrum(bands, carriers, otf, params)
[h, w] = size(bands.C0);
outputSize = [2*h, 2*w];

otf0 = otf.values;
directCombinedSpectrum = placeSpectrumAtCenterLocal(fft2cLocal(bands.C0), outputSize) ...
    .* otf0 .* otfAttenuationMaskLocal(otf0, params);

components = {
    bands.CsPlus,  carriers.ksPixel, params.modulationS
    bands.CsMinus, -carriers.ksPixel, params.modulationS
    bands.CtPlus,  carriers.ktPixel, params.modulationT
    bands.CtMinus, -carriers.ktPixel, params.modulationT
};

for idx = 1:size(components, 1)
    component = components{idx, 1};
    carrierPixel = components{idx, 2};
    sidebandAmplitude = components{idx, 3} / 2;
    freq = placeSpectrumAtCenterLocal(fft2cLocal(component), outputSize);
    freq = shiftSpectrumOnCanvasLocal(freq, -carrierPixel);
    shiftedOtf = shiftOtfByCarrierLocal(otf.values, -carrierPixel);
    attenuation = otfAttenuationMaskLocal(shiftedOtf, params);
    directCombinedSpectrum = directCombinedSpectrum + ...
        freq .* shiftedOtf .* attenuation ./ sidebandAmplitude;
end
end

function output = placeSpectrumAtCenterLocal(spectrum, outputSize)
[h, w] = size(spectrum);
output = zeros(outputSize);
rowStart = floor((outputSize(1) - h) / 2) + 1;
colStart = floor((outputSize(2) - w) / 2) + 1;
output(rowStart:rowStart+h-1, colStart:colStart+w-1) = spectrum;
end

function shifted = shiftSpectrumOnCanvasLocal(spectrum, shiftPixel)
[h, w] = size(spectrum);
[x, y] = meshgrid(0:w-1, 0:h-1);
image = ifft2cLocal(spectrum);
phaseRamp = exp(2i*pi * (shiftPixel(1) * x / w + shiftPixel(2) * y / h));
shifted = fft2cLocal(image .* phaseRamp);
end

function shiftedOtf = shiftOtfByCarrierLocal(otfValues, carrierPixel)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);
end

function attenuation = otfAttenuationMaskLocal(otfValues, params)
attenuation = ones(size(otfValues));
if ~isfield(params, 'otfAttenuationEnabled') || ~params.otfAttenuationEnabled
    return;
end

[~, maxIdx] = max(abs(otfValues(:)));
[centerRow, centerCol] = ind2sub(size(otfValues), maxIdx);
[h, w] = size(otfValues);
[x, y] = meshgrid(1:w, 1:h);
radius = hypot((x - centerCol) * 2 / w, (y - centerRow) * 2 / h);
attenuation = 1 - params.otfAttenuationStrength .* ...
    exp(-(radius .^ 2) ./ ((0.5 * params.otfAttenuationFwhm) ^ 2));
attenuation = min(max(attenuation, 0), 1);
end

function image = spectrumForDisplay(spectrum)
image = log(1 + abs(spectrum));
image = image - min(image(:));
maxValue = max(image(:));
if maxValue > 0
    image = image ./ maxValue;
end
end

function spectrum = fft2cLocal(image)
spectrum = fftshift(fft2(ifftshift(image)));
end

function image = ifft2cLocal(spectrum)
image = fftshift(ifft2(ifftshift(spectrum)));
end
