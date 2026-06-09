function test_lattice_shifted_component_visualization()
%TEST_LATTICE_SHIFTED_COMPONENT_VISUALIZATION Export shifted component spectra for current data.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

stack = readCurrentFiveFrameStack(fullfile(coreDir, 'roi.tif'));
params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 610;
params.NA = 1.42;
params.normalizeFrames = true;
params.phaseOffsetS = -3.0063;
params.phaseOffsetT = -0.2506;
[normalizedStack, ~] = normalizeSIMFrames(stack, params);
stackFFT = zeros(size(normalizedStack));
for idx = 1:5
    stackFFT(:, :, idx) = FFT2D(normalizedStack(:, :, idx), false);
end
bands = separateLatticeBandsFrequency(stackFFT, params);
[carriers, ~] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(stackFFT, 1), size(stackFFT, 2), params);
[~, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params);

assert(isfield(diagnostics, 'shiftedBandSpectra'), ...
    'Fusion diagnostics should expose the five shifted component spectra before Wiener fusion.');
assert(isfield(diagnostics, 'shiftedBandNames'), ...
    'Fusion diagnostics should name the shifted component spectra.');
assert(isfield(diagnostics, 'shiftedBandShiftPixels'), ...
    'Fusion diagnostics should record the pixel shift applied to each component.');
assert(isfield(diagnostics, 'shiftedBandPeakPixels'), ...
    'Fusion diagnostics should record where each shifted component peaks on the output canvas.');
assert(isfield(diagnostics, 'fusionNumeratorContributions'), ...
    'Fusion diagnostics should expose the five OTF-weighted spectra that enter fftDirectlyCombined.');
assert(numel(diagnostics.shiftedBandSpectra) == 5, ...
    'Expected C0, CsPlus, CsMinus, CtPlus, and CtMinus shifted spectra.');
assert(numel(diagnostics.fusionNumeratorContributions) == 5, ...
    'Expected one OTF-weighted fusion numerator contribution for each shifted spectrum.');
assert(isequal(size(diagnostics.shiftedBandShiftPixels), [5, 2]), ...
    'Expected one [x, y] shift vector per shifted component.');

outputDir = fullfile(coreDir, 'output', 'diagnostics', 'shifted_components');
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

displayTiles = cell(1, 5);
weightedTiles = cell(1, 5);
weightedSum = zeros(size(diagnostics.simSpectrum));
for idx = 1:5
    spectrum = diagnostics.shiftedBandSpectra{idx};
    weightedSpectrum = diagnostics.fusionNumeratorContributions{idx};
    assert(isequal(size(spectrum), size(diagnostics.simSpectrum)), ...
        'Each shifted component spectrum should live on the same 2x fusion canvas as the final SIM spectrum.');
    assert(isequal(size(weightedSpectrum), size(diagnostics.fftDirectlyCombined)), ...
        'Each OTF-weighted contribution should live on the same canvas as fftDirectlyCombined.');
    assert(any(abs(spectrum(:)) > 0), ...
        'Shifted component spectrum should not be empty.');
    assert(any(abs(weightedSpectrum(:)) > 0), ...
        'OTF-weighted contribution should not be empty.');
    displayTiles{idx} = spectrumToUint16LogImage(spectrum);
    weightedTiles{idx} = spectrumToUint16LogImage(weightedSpectrum);
    weightedSum = weightedSum + weightedSpectrum;
    imwrite(displayTiles{idx}, fullfile(outputDir, ...
        sprintf('%02d_%s_shifted_spectrum.tif', idx, diagnostics.shiftedBandNames{idx})), 'tif');
    imwrite(weightedTiles{idx}, fullfile(outputDir, ...
        sprintf('%02d_%s_fusion_numerator_contribution.tif', idx, diagnostics.shiftedBandNames{idx})), 'tif');
end
assert(max(abs(weightedSum(:) - diagnostics.fftDirectlyCombined(:))) < 1e-10 * ...
    max(1, max(abs(diagnostics.fftDirectlyCombined(:)))), ...
    'The five OTF-weighted contributions should sum exactly to fftDirectlyCombined.');

montageImage = makeFiveTileMontage(displayTiles);
imwrite(montageImage, fullfile(outputDir, 'shifted_component_spectra_montage.tif'), 'tif');
weightedMontageImage = makeFiveTileMontage(weightedTiles);
imwrite(weightedMontageImage, fullfile(outputDir, 'fusion_numerator_contributions_montage.tif'), 'tif');

shiftedBandNames = diagnostics.shiftedBandNames;
shiftedBandShiftPixels = diagnostics.shiftedBandShiftPixels;
shiftedBandPeakPixels = diagnostics.shiftedBandPeakPixels;
shiftedBandCentroidPixels = diagnostics.shiftedBandCentroidPixels;
save(fullfile(outputDir, 'shifted_component_positions.mat'), ...
    'shiftedBandNames', 'shiftedBandShiftPixels', ...
    'shiftedBandPeakPixels', 'shiftedBandCentroidPixels');
writePositionSummary(outputDir, shiftedBandNames, shiftedBandShiftPixels, ...
    shiftedBandPeakPixels, shiftedBandCentroidPixels);
end

function stack = readCurrentFiveFrameStack(stackPath)
stack = [];
for idx = 1:5
    frame = double(imread(stackPath, idx));
    if idx == 1
        stack = zeros(size(frame, 1), size(frame, 2), 5);
    end
    stack(:, :, idx) = frame;
end
end

function image = spectrumToUint16LogImage(spectrum)
displayImage = log10(1 + abs(spectrum));
displayImage = displayImage - min(displayImage(:));
maxValue = max(displayImage(:));
if maxValue > 0
    displayImage = displayImage ./ maxValue;
end
image = uint16(round(displayImage * double(intmax('uint16'))));
end

function montageImage = makeFiveTileMontage(tiles)
[h, w] = size(tiles{1});
montageImage = zeros(2*h, 3*w, 'uint16');
for idx = 1:5
    row = floor((idx - 1) / 3);
    col = mod(idx - 1, 3);
    rows = row*h + (1:h);
    cols = col*w + (1:w);
    montageImage(rows, cols) = tiles{idx};
end
end

function writePositionSummary(outputDir, names, shifts, peaks, centroids)
fid = fopen(fullfile(outputDir, 'shifted_component_positions.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'name\tshift_x_px\tshift_y_px\tpeak_x_px\tpeak_y_px\tcentroid_x_px\tcentroid_y_px\n');
for idx = 1:numel(names)
    fprintf(fid, '%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n', ...
        names{idx}, shifts(idx, 1), shifts(idx, 2), peaks(idx, 1), peaks(idx, 2), ...
        centroids(idx, 1), centroids(idx, 2));
end
end
