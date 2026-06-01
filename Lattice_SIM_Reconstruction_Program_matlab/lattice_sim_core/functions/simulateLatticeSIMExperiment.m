function [rawStack, truth] = simulateLatticeSIMExperiment(objectImage, simParams)
%SIMULATELATTICESIMEXPERIMENT Generate OpenSIM-style five-frame Lattice-SIM data.

if nargin < 1 || isempty(objectImage)
    objectImage = makeLatticeSIMPhantom([256, 256]);
end
if nargin < 2
    simParams = struct();
end

objectImage = double(objectImage);
if ndims(objectImage) ~= 2 || isempty(objectImage) || any(~isfinite(objectImage(:)))
    error('LatticeSIM:InvalidSimulationParameter', 'objectImage must be a finite 2-D numeric image.');
end

simParams = defaultSimulationParams(size(objectImage), simParams);
validateSimulationParams(simParams, size(objectImage));

if ~isempty(simParams.randomSeed)
    previousRng = rng();
    rng(simParams.randomSeed);
    rngCleanup = onCleanup(@() rng(previousRng));
end

[h, w] = size(objectImage);
[x, y] = meshgrid(0:w-1, 0:h-1);

phasePairs = latticePhasePairs();

phaseErrors = simParams.phaseErrorStd .* randn(5, 2);
illumination = zeros(h, w, 5);
noiseFreeStack = zeros(h, w, 5);
rawStack = zeros(h, w, 5);

otf = buildSimulationOTF(h, w, simParams);

ks = simParams.ksPixel;
kt = simParams.ktPixel;
phaseSBase = 2*pi*(ks(1) * x / w + ks(2) * y / h);
phaseTBase = 2*pi*(kt(1) * x / w + kt(2) * y / h);

for frameIdx = 1:5
    phiS = phasePairs(frameIdx, 1) + phaseErrors(frameIdx, 1);
    phiT = phasePairs(frameIdx, 2) + phaseErrors(frameIdx, 2);
    illumination(:, :, frameIdx) = simParams.meanIllumination ...
        + simParams.modulationS * cos(phaseSBase + phiS) ...
        + simParams.modulationT * cos(phaseTBase + phiT);
    illuminated = objectImage .* illumination(:, :, frameIdx);
    if simParams.useOTF
        imageFrame = real(ifft2(fft2(illuminated) .* fftshift(otf)));
    else
        imageFrame = illuminated;
    end
    noiseFreeStack(:, :, frameIdx) = max(imageFrame, 0);
    sigma = simParams.noiseLevel * std(noiseFreeStack(:, :, frameIdx), 0, 'all');
    if sigma > 0
        noisyFrame = noiseFreeStack(:, :, frameIdx) + sigma .* randn(h, w);
    else
        noisyFrame = noiseFreeStack(:, :, frameIdx);
    end
    rawStack(:, :, frameIdx) = max(noisyFrame, 0);
end

truth.object = objectImage;
truth.illumination = illumination;
truth.noiseFreeStack = noiseFreeStack;
truth.otf = otf;
truth.ksPixel = simParams.ksPixel;
truth.ktPixel = simParams.ktPixel;
truth.phasePairs = phasePairs;
truth.phaseErrors = phaseErrors;
truth.params = simParams;
end

function params = defaultSimulationParams(imageSize, userParams)
params.imageSize = imageSize;
params.ksPixel = [8, 0];
params.ktPixel = [0, 10];
params.modulationS = 0.45;
params.modulationT = 0.40;
params.meanIllumination = 1.0;
params.pixelSizeNm = 65;
params.emissionWavelengthNm = 532;
params.NA = 1.2;
params.noiseLevel = 0.02;
params.phaseErrorStd = 0.0;
params.randomSeed = 1;
params.useOTF = true;

names = fieldnames(userParams);
for idx = 1:numel(names)
    params.(names{idx}) = userParams.(names{idx});
end
end

function validateSimulationParams(params, objectSize)
if ~isequal(params.imageSize, objectSize)
    error('LatticeSIM:InvalidSimulationParameter', 'params.imageSize must match objectImage size.');
end
positiveFields = {'modulationS', 'modulationT', 'meanIllumination', ...
    'pixelSizeNm', 'emissionWavelengthNm', 'NA'};
for idx = 1:numel(positiveFields)
    name = positiveFields{idx};
    if ~isnumeric(params.(name)) || ~isscalar(params.(name)) || params.(name) <= 0
        error('LatticeSIM:InvalidSimulationParameter', 'params.%s must be a positive scalar.', name);
    end
end
if params.meanIllumination <= params.modulationS + params.modulationT
    error('LatticeSIM:InvalidSimulationParameter', ...
        'meanIllumination must exceed modulationS + modulationT to keep illumination nonnegative.');
end
if params.noiseLevel < 0 || params.phaseErrorStd < 0
    error('LatticeSIM:InvalidSimulationParameter', 'Noise and phase-error levels must be nonnegative.');
end
if numel(params.ksPixel) ~= 2 || numel(params.ktPixel) ~= 2
    error('LatticeSIM:InvalidSimulationParameter', 'ksPixel and ktPixel must be two-element vectors.');
end
end

function otf = buildSimulationOTF(h, w, params)
[fx, fy] = simulationFrequencyGrid(h, w, params.pixelSizeNm);
cutoff = 2 * params.NA / params.emissionWavelengthNm;
rho = hypot(fx, fy) / cutoff;
otf = zeros(h, w);
inside = rho <= 1;
otf(inside) = (2 / pi) * (acos(rho(inside)) - rho(inside) .* sqrt(1 - rho(inside).^2));
end

function [fx, fy] = simulationFrequencyGrid(h, w, pixelSizeNm)
fxAxis = ((1:w) - floor(w/2) - 1) ./ (w * pixelSizeNm);
fyAxis = ((1:h) - floor(h/2) - 1) ./ (h * pixelSizeNm);
[fx, fy] = meshgrid(fxAxis, fyAxis);
end
