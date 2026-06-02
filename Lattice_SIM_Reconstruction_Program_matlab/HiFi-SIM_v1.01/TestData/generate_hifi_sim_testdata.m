function [rawStack, truth, outputs] = generate_hifi_sim_testdata(outputDir, objectImage, simParams)
%GENERATE_HIFI_SIM_TESTDATA Generate a 2D HiFi-SIM 3-direction, 3-phase test stack.
%
% With no inputs, this function writes a synthetic 9-frame raw stack to:
%   HiFi-SIM_v1.01/TestData/generated_hifi_sim_testdata/
%
% The frame order matches HiFiSIM.m:
%   frames 1-3: direction 1, phases 0, 2*pi/3, 4*pi/3
%   frames 4-6: direction 2, phases 0, 2*pi/3, 4*pi/3
%   frames 7-9: direction 3, phases 0, 2*pi/3, 4*pi/3

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), 'generated_hifi_sim_testdata');
end
if nargin < 2 || isempty(objectImage)
    objectImage = defaultObjectImage();
end
if nargin < 3
    simParams = struct();
end

objectImage = double(objectImage);
if ndims(objectImage) ~= 2 || isempty(objectImage) || any(~isfinite(objectImage(:)))
    error('HiFiSIM:InvalidSimulationInput', 'objectImage must be a finite 2-D numeric image.');
end

params = defaultHiFiSimulationParams(size(objectImage), simParams);
validateHiFiSimulationParams(params, size(objectImage));

if ~isempty(params.randomSeed)
    previousRng = rng();
    rng(params.randomSeed);
    cleanupRng = onCleanup(@() rng(previousRng));
end

[h, w] = size(objectImage);
[x, y] = meshgrid(0:w-1, 0:h-1);
otf = buildSimulationOTF(h, w, params);

nrFrames = params.nrDirs * params.nrPhases;
rawStack = zeros(h, w, nrFrames);
noiseFreeStack = zeros(h, w, nrFrames);
illumination = zeros(h, w, nrFrames);
carrierPixels = zeros(params.nrDirs, 2);
directionIndex = zeros(nrFrames, 1);
phaseIndex = zeros(nrFrames, 1);

frameIdx = 0;
for dirIdx = 1:params.nrDirs
    carrierPixels(dirIdx, :) = params.carrierMagnitudePixels .* ...
        [cosd(params.directionAnglesDeg(dirIdx)), sind(params.directionAnglesDeg(dirIdx))];
    phaseBase = 2*pi*(carrierPixels(dirIdx, 1) * x / w + carrierPixels(dirIdx, 2) * y / h);

    for phaseIdx = 1:params.nrPhases
        frameIdx = frameIdx + 1;
        phase = params.phaseSteps(phaseIdx) + params.phaseErrorStd .* randn();
        directionIndex(frameIdx) = dirIdx;
        phaseIndex(frameIdx) = phaseIdx;

        illumination(:, :, frameIdx) = params.meanIllumination + ...
            params.modulation .* cos(phaseBase + phase);
        illuminated = objectImage .* illumination(:, :, frameIdx);
        if params.useOTF
            imageFrame = real(ifft2(fft2(illuminated) .* fftshift(otf)));
        else
            imageFrame = illuminated;
        end

        noiseFreeStack(:, :, frameIdx) = max(imageFrame, 0);
        sigma = params.noiseLevel * std(noiseFreeStack(:, :, frameIdx), 0, 'all');
        if sigma > 0
            rawFrame = noiseFreeStack(:, :, frameIdx) + sigma .* randn(h, w);
        else
            rawFrame = noiseFreeStack(:, :, frameIdx);
        end
        rawStack(:, :, frameIdx) = max(rawFrame, 0);
    end
end

truth.object = objectImage;
truth.illumination = illumination;
truth.noiseFreeStack = noiseFreeStack;
truth.otf = otf;
truth.carrierPixels = carrierPixels;
truth.directionAnglesDeg = params.directionAnglesDeg;
truth.phaseSteps = params.phaseSteps;
truth.directionIndex = directionIndex;
truth.phaseIndex = phaseIndex;
truth.params = params;

outputs.outputDir = outputDir;
outputs.rawStackTiff = fullfile(outputDir, 'simulated_hifi_2d_3dir_3phase_9frames.tif');
outputs.objectTiff = fullfile(outputDir, 'synthetic_object.tif');
outputs.widefieldTiff = fullfile(outputDir, 'widefield_average.tif');
outputs.truthMat = fullfile(outputDir, 'diagnostics', 'hifi_simulation_truth.mat');

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end
diagnosticsDir = fullfile(outputDir, 'diagnostics');
if exist(diagnosticsDir, 'dir') ~= 7
    mkdir(diagnosticsDir);
end

imwrite(toUint16Image(objectImage), outputs.objectTiff, 'tif');
writeStackTiff(rawStack, outputs.rawStackTiff);
imwrite(toUint16Image(mean(rawStack, 3)), outputs.widefieldTiff, 'tif');
save(outputs.truthMat, 'truth', 'rawStack');

fprintf('HiFi-SIM simulated raw stack: %s\n', outputs.rawStackTiff);
fprintf('HiFi-SIM simulation truth: %s\n', outputs.truthMat);
end

function objectImage = defaultObjectImage()
thisDir = fileparts(mfilename('fullpath'));
testpatPath = fullfile(thisDir, '..', '..', 'lattice_sim_core', 'testpat.tiff');
cropSize = 256;
if exist(testpatPath, 'file') == 2
    fullImage = double(imread(testpatPath));
    rowStart = floor((size(fullImage, 1) - cropSize) / 2) + 1;
    colStart = floor((size(fullImage, 2) - cropSize) / 2) + 1;
    objectImage = fullImage(rowStart:rowStart+cropSize-1, colStart:colStart+cropSize-1);
else
    objectImage = peaks(cropSize);
end
objectImage = objectImage - min(objectImage(:));
maxValue = max(objectImage(:));
if maxValue > 0
    objectImage = objectImage ./ maxValue;
end
end

function params = defaultHiFiSimulationParams(imageSize, userParams)
params.imageSize = imageSize;
params.nrDirs = 3;
params.nrPhases = 3;
params.directionAnglesDeg = [0, 60, 120];
params.phaseSteps = [0, 2*pi/3, 4*pi/3];
params.pixelSizeNm = 78.6;
params.emissionWavelengthNm = 525;
params.NA = 1.42;
params.carrierCutoffFraction = 0.55;
params.carrierMagnitudePixels = [];
params.modulation = 0.75;
params.meanIllumination = 1.0;
params.noiseLevel = 0.02;
params.phaseErrorStd = 0.0;
params.randomSeed = 13;
params.useOTF = true;

names = fieldnames(userParams);
for idx = 1:numel(names)
    params.(names{idx}) = userParams.(names{idx});
end

if isempty(params.carrierMagnitudePixels)
    cutoffCyclesPerNm = 2 * params.NA / params.emissionWavelengthNm;
    cyclesPerPixelBin = 1 / (imageSize(2) * params.pixelSizeNm);
    params.carrierMagnitudePixels = params.carrierCutoffFraction * ...
        cutoffCyclesPerNm / cyclesPerPixelBin;
end
end

function validateHiFiSimulationParams(params, objectSize)
if ~isequal(params.imageSize, objectSize)
    error('HiFiSIM:InvalidSimulationParameter', 'params.imageSize must match objectImage size.');
end
if params.nrDirs ~= 3 || params.nrPhases ~= 3
    error('HiFiSIM:InvalidSimulationParameter', 'This generator targets 3 directions x 3 phases.');
end
if numel(params.directionAnglesDeg) ~= params.nrDirs
    error('HiFiSIM:InvalidSimulationParameter', 'directionAnglesDeg must contain 3 values.');
end
if numel(params.phaseSteps) ~= params.nrPhases
    error('HiFiSIM:InvalidSimulationParameter', 'phaseSteps must contain 3 values.');
end
positiveFields = {'pixelSizeNm', 'emissionWavelengthNm', 'NA', ...
    'carrierMagnitudePixels', 'meanIllumination'};
for idx = 1:numel(positiveFields)
    name = positiveFields{idx};
    if ~isnumeric(params.(name)) || ~isscalar(params.(name)) || params.(name) <= 0
        error('HiFiSIM:InvalidSimulationParameter', 'params.%s must be a positive scalar.', name);
    end
end
if params.modulation < 0 || params.meanIllumination <= params.modulation
    error('HiFiSIM:InvalidSimulationParameter', ...
        'meanIllumination must exceed nonnegative modulation.');
end
if params.noiseLevel < 0 || params.phaseErrorStd < 0
    error('HiFiSIM:InvalidSimulationParameter', 'Noise and phase-error levels must be nonnegative.');
end
end

function otf = buildSimulationOTF(h, w, params)
fxAxis = ((1:w) - floor(w/2) - 1) ./ (w * params.pixelSizeNm);
fyAxis = ((1:h) - floor(h/2) - 1) ./ (h * params.pixelSizeNm);
[fx, fy] = meshgrid(fxAxis, fyAxis);
cutoff = 2 * params.NA / params.emissionWavelengthNm;
rho = hypot(fx, fy) / cutoff;
otf = zeros(h, w);
inside = rho <= 1;
otf(inside) = (2 / pi) * (acos(rho(inside)) - rho(inside) .* sqrt(1 - rho(inside).^2));
end

function writeStackTiff(stack, outputPath)
if exist(outputPath, 'file') == 2
    delete(outputPath);
end
for idx = 1:size(stack, 3)
    frame = toUint16Image(stack(:, :, idx));
    if idx == 1
        imwrite(frame, outputPath, 'tif');
    else
        imwrite(frame, outputPath, 'tif', 'WriteMode', 'append');
    end
end
end

function image = toUint16Image(image)
image = double(image);
image = image - min(image(:));
maxValue = max(image(:));
if maxValue > 0
    image = image ./ maxValue;
end
image = uint16(round(65535 .* image));
end
