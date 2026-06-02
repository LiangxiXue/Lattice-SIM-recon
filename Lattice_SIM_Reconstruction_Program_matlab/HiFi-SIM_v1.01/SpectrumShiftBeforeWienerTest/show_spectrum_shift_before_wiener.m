function debugResult = show_spectrum_shift_before_wiener(rawTiffPath)
%SHOW_SPECTRUM_SHIFT_BEFORE_WIENER Display HiFi-SIM spectrum-combined image before Wiener filtering.
%
% Usage:
%   show_spectrum_shift_before_wiener
%   show_spectrum_shift_before_wiener('/path/to/raw_sim_stack.tif')
%
% This script follows the HiFi-SIM GUI path through frequency shifting and
% direct spectrum combination, then stops before Wk1/Wk2 Wiener filtering.

scriptDir = fileparts(mfilename('fullpath'));
hifiDir = fileparts(scriptDir);
mainFunDir = fullfile(hifiDir, 'Main_fun');
addpath(mainFunDir);

if nargin < 1 || isempty(rawTiffPath)
    rawTiffPath = fullfile(hifiDir, 'TestData', '2D-SIM(3 angles_3 phases)_9 frames_Group1.tif');
end

outputDir = fullfile(scriptDir, 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

param = default2DSIMParams();
[param.Iraw, param.Size1, param.Size2, param.imgSize, param.Format] = loadRawStack(rawTiffPath, param);
param.filename = rawTiffPath;
param.cyclesPerMicron = 1 / (param.imgSize * param.micronsPerPixel);
param.attStrength = 0;
param.OtfProvider = SimOtfProvider(param, param.NA, param.lambda, 1);

fprintf('Loaded raw stack: %s\n', rawTiffPath);
fprintf('Image size: %d x %d, frames: %d\n', param.Size2, param.Size1, size(param.Iraw, 3));

IIrawFFT = preprocessToFFT(param);
[param, K0] = estimateParametersLikeHiFiSIM(param, IIrawFFT);
[fftDirectlyCombined, shiftedBands] = combineShiftedSpectrumLikeHiFiSIM(param, IIrawFFT);

beforeWienerImage = real(ifft2(fftshift(fftDirectlyCombined)));
beforeWienerImage(beforeWienerImage < 0) = 0;

frequencyMagnitude = log(1 + abs(fftDirectlyCombined));
frequencyMagnitude = frequencyMagnitude ./ max(frequencyMagnitude(:));

displayImage = normalizeForDisplay(beforeWienerImage);
displaySpectrum = normalizeForDisplay(frequencyMagnitude);

imwrite(displayImage, fullfile(outputDir, 'before_wiener_spatial.png'));
imwrite(displaySpectrum, fullfile(outputDir, 'before_wiener_spectrum.png'));
save(fullfile(outputDir, 'before_wiener_debug.mat'), ...
    'rawTiffPath', 'param', 'K0', 'IIrawFFT', 'shiftedBands', ...
    'fftDirectlyCombined', 'beforeWienerImage', 'frequencyMagnitude');

figure('Name', 'HiFi-SIM before Wiener: spatial image');
imagesc(beforeWienerImage);
axis image off;
colormap hot;
colorbar;
title('After spectrum shift/combination, before Wiener');

figure('Name', 'HiFi-SIM before Wiener: spectrum magnitude');
imagesc(frequencyMagnitude);
axis image off;
colormap hot;
colorbar;
title('log(1 + abs(fftDirectlyCombined))');

debugResult = struct();
debugResult.rawTiffPath = rawTiffPath;
debugResult.outputDir = outputDir;
debugResult.param = param;
debugResult.K0 = K0;
debugResult.fftDirectlyCombined = fftDirectlyCombined;
debugResult.beforeWienerImage = beforeWienerImage;
debugResult.frequencyMagnitude = frequencyMagnitude;

fprintf('Saved spatial image: %s\n', fullfile(outputDir, 'before_wiener_spatial.png'));
fprintf('Saved spectrum image: %s\n', fullfile(outputDir, 'before_wiener_spectrum.png'));
fprintf('Saved debug MAT: %s\n', fullfile(outputDir, 'before_wiener_debug.mat'));
end

function param = default2DSIMParams()
param = struct();
param.nrDirs = 3;
param.nrPhases = 3;
param.nrBands = 2;
param.phaOff = 0;
param.fac = ones(1, 2);
param.micronsPerPixel = 78.6e-3;
param.lambda = 525;
param.NA = 1.42;
param.attStrength = 0.90;
end

function [Iraw, size1, size2, imgSize, format] = loadRawStack(rawTiffPath, param)
info = imfinfo(rawTiffPath);
expectedFrames = param.nrDirs * param.nrPhases;
if numel(info) ~= expectedFrames
    error('HiFiSIMDebug:FrameCountMismatch', ...
        'Expected %d frames, but %s has %d frames.', expectedFrames, rawTiffPath, numel(info));
end

imgSize = max(info(1).Width, info(1).Height);
Iraw0 = zeros(info(1).Height, info(1).Width, expectedFrames);
Iraw = zeros(imgSize, imgSize, expectedFrames);
for frameIdx = 1:expectedFrames
    Iraw0(:, :, frameIdx) = double(imread(rawTiffPath, frameIdx));
    if info(1).Width == info(1).Height || info(1).Width > info(1).Height
        Iraw(1:info(1).Height, :, frameIdx) = Iraw0(1:info(1).Height, :, frameIdx);
    else
        Iraw(:, 1:info(1).Width, frameIdx) = Iraw0(:, 1:info(1).Width, frameIdx);
    end
end

size1 = info(1).Height;
size2 = info(1).Width;
format = info(1).Format;
end

function IIrawFFT = preprocessToFFT(param)
psf = abs(otf2psf(param.OtfProvider.otf));
temp = importImages(param.Iraw);
IIraw = deconvlucy(temp, psf, 5);

frameCount = param.nrDirs * param.nrPhases;
IIrawFFT = zeros(size(IIraw, 1), size(IIraw, 2), frameCount);
for frameIdx = 1:frameCount
    IIrawFFT(:, :, frameIdx) = FFT2D(IIraw(:, :, frameIdx), false);
end
end

function [param, K0] = estimateParametersLikeHiFiSIM(param, IIrawFFT)
NPixel = size(param.Iraw, 1);
cnt = [NPixel / 2 + 1, NPixel / 2 + 1];
param.cutoff = 1000 / (0.5 * param.lambda / param.NA);
[x, y] = meshgrid(1:NPixel, 1:NPixel);
rad = sqrt((y - cnt(1)).^2 + (x - cnt(2)).^2);
mask = double(rad <= 1.0 * (param.cutoff / param.cyclesPerMicron + 1));
notchFilter0 = getotfAtt(NPixel, param.OtfProvider.cyclesPerMicron, 0.5 * param.cutoff, 0, 0);
notchFilter = notchFilter0 .* mask;
mask2 = double(rad <= 1.10 * (param.cutoff / param.cyclesPerMicron + 1));
notchFilter2 = notchFilter0 .* mask2;

crossCorrelation = zeros(size(mask2, 1), size(mask2, 2), param.nrDirs);
k0 = zeros(1, param.nrDirs);
separateCache = cell(1, param.nrDirs);

for dirIdx = 1:param.nrDirs
    frameRange = (dirIdx - 1) * param.nrPhases + 1:dirIdx * param.nrPhases;
    separateII = separateBands(IIrawFFT(:, :, frameRange), 0, param.nrBands, ones(1, param.nrBands));
    separateCache{dirIdx} = separateII;

    c0 = separateII(:, :, 1);
    c2 = separateII(:, :, 2);
    c0 = c0 ./ max(abs(c0(:)));
    c2 = c2 ./ max(abs(c2(:)));
    c0 = c0 .* notchFilter;
    c2 = c2 .* notchFilter;
    c0 = FFT2D(c0, false);
    c2 = FFT2D(c2, false);
    c2 = c2 .* conj(c0);
    c2 = c2 ./ max(c2(:));
    vec = fftshift(FFT2D(c2, true));
    crossCorrelation(:, :, dirIdx) = vec;

    temp = vec .* notchFilter2;
    temp = log(1 + abs(temp));
    temp = temp ./ max(temp(:));
    [yPos, xPos] = find(temp == max(temp(:)));
    k0(dirIdx) = sqrt((xPos(1) - cnt(1))^2 + (yPos(1) - cnt(2))^2);
end

flag = 0;
if param.nrDirs > 2 && max(k0) - min(k0) > 8
    flag = 1;
    Kobject = min(k0);
    maskOuter = rad >= (Kobject + 1);
    maskInner = rad <= (Kobject - 1);
end

K0 = zeros(1, param.nrDirs);
for dirIdx = 1:param.nrDirs
    vec = crossCorrelation(:, :, dirIdx);
    if flag == 1
        vec(maskOuter) = 0;
        vec(maskInner) = 0;
    end
    temp = vec .* notchFilter2;
    temp = log(1 + abs(temp));
    temp = temp ./ max(temp(:));
    [yPos, xPos] = find(temp == max(temp(:)));

    peak.xPos = xPos(1);
    peak.yPos = yPos(1);
    cntrl = zeros(10, 30);
    overlap = 0.15;
    step = 2.5;
    bn1 = (param.nrBands - 1) * 2;
    kx = peak.xPos - cnt(2);
    ky = peak.yPos - cnt(1);

    separateII = separateCache{dirIdx};
    c0 = separateII(:, :, 1) ./ max(abs(reshape(separateII(:, :, 1), [], 1)));
    c2 = separateII(:, :, 2) ./ max(abs(reshape(separateII(:, :, 2), [], 1)));
    [peak, ~] = fitPeak(c0, c2, 1, bn1, param.OtfProvider, -kx, -ky, overlap, step, cntrl);

    p1 = getPeak(separateII(:, :, 1), separateII(:, :, 2), 1, 2, ...
        param.OtfProvider, peak.kx, peak.ky, overlap);
    param.Dir(dirIdx).px = -peak.kx;
    param.Dir(dirIdx).py = -peak.ky;
    param.Dir(dirIdx).phaOff = -phase(p1);
    modulation = abs(p1);
    if modulation > 1.0
        modulation = 1.0;
    end
    if modulation < 0.35
        modulation = 0.7;
    end
    param.Dir(dirIdx).modul = modulation;
    K0(dirIdx) = sqrt(param.Dir(dirIdx).px^2 + param.Dir(dirIdx).py^2);
end
end

function [fftDirectlyCombined, shiftedBands] = combineShiftedSpectrumLikeHiFiSIM(param, IIrawFFT)
h = size(param.Iraw, 1);
w = size(param.Iraw, 2);
fftDirectlyCombined = zeros(h * 2, w * 2);
shiftedBands = cell(1, param.nrDirs);

param.attStrength = 0.90;
param.a = 1.0;
param.attFWHM = 1.0;
param.OtfProvider = SimOtfProvider(param, param.NA, param.lambda, param.a);

for dirIdx = 1:param.nrDirs
    par = param.Dir(dirIdx);
    fac = ones(1, param.nrBands);
    fac(2:param.nrBands) = param.Dir(dirIdx).modul(1:param.nrBands - 1);
    frameRange = (dirIdx - 1) * param.nrPhases + 1:dirIdx * param.nrPhases;
    separate = separateBands(IIrawFFT(:, :, frameRange), par.phaOff, param.nrBands, fac);

    shifted = zeros(2 * h, 2 * w, param.nrPhases);
    shifted(:, :, 1) = placeFreq(separate(:, :, 1));
    for bandIdx = 2:param.nrBands
        pos = bandIdx * 2 - 2;
        neg = bandIdx * 2 - 1;
        shifted(:, :, pos) = placeFreq(separate(:, :, pos));
        shifted(:, :, neg) = placeFreq(separate(:, :, neg));
        shifted(:, :, pos) = NfourierShift(shifted(:, :, pos), ...
            -(bandIdx - 1) * par.px, -(bandIdx - 1) * par.py);
        shifted(:, :, neg) = NfourierShift(shifted(:, :, neg), ...
            (bandIdx - 1) * par.px, (bandIdx - 1) * par.py);
    end

    shifted(:, :, 1) = applyOtf(shifted(:, :, 1), param.OtfProvider, 1, 0, 0, 1, 0);
    for bandIdx = 2:param.nrBands
        pos = bandIdx * 2 - 2;
        neg = bandIdx * 2 - 1;
        shifted(:, :, pos) = applyOtf(shifted(:, :, pos), param.OtfProvider, ...
            bandIdx, -(bandIdx - 1) * par.px, -(bandIdx - 1) * par.py, 1, 0);
        shifted(:, :, neg) = applyOtf(shifted(:, :, neg), param.OtfProvider, ...
            bandIdx, (bandIdx - 1) * par.px, (bandIdx - 1) * par.py, 1, 0);
    end

    for bandImageIdx = 1:param.nrBands * 2 - 1
        fftDirectlyCombined = fftDirectlyCombined + shifted(:, :, bandImageIdx);
    end
    shiftedBands{dirIdx} = shifted;
end
end

function image = normalizeForDisplay(image)
image = real(image);
image = image - min(image(:));
maxValue = max(image(:));
if maxValue > 0
    image = image ./ maxValue;
end
end
