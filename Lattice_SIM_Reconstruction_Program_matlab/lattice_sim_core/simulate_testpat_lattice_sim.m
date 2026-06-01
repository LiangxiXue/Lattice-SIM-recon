%SIMULATE_TESTPAT_LATTICE_SIM Simulate and reconstruct Lattice-SIM from testpat.tiff.
%
% Outputs are written to:
%   testpat_lattice_simulation_output/
%     cropped_testpat.tif
%     simulated_raw_stack.tif
%     Wide-field-C0.tif
%     Wide-field-D3D4D5.tif
%     Lattice-SIM.tif
%     diagnostics/
%       rawStack.mat
%       result.mat
%       simulation_truth.mat

clear;
clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'functions'));

inputPath = fullfile(rootDir, 'testpat.tiff');
if exist(inputPath, 'file') ~= 2
    error('LatticeSIM:MissingInputImage', 'Expected test image not found: %s', inputPath);
end

objectImageFull = double(imread(inputPath));
cropSize = 256;
if size(objectImageFull, 1) < cropSize || size(objectImageFull, 2) < cropSize
    error('LatticeSIM:InvalidInputImage', 'testpat.tiff must be at least %d x %d pixels.', cropSize, cropSize);
end

rowStart = floor((size(objectImageFull, 1) - cropSize) / 2) + 1;
colStart = floor((size(objectImageFull, 2) - cropSize) / 2) + 1;
objectImage = objectImageFull(rowStart:rowStart+cropSize-1, colStart:colStart+cropSize-1);

simParams.imageSize = size(objectImage);
simParams.expansionFactor = 1.6;
simParams.modulationS = 0.45;
simParams.modulationT = 0.40;
simParams.meanIllumination = 1.0;
simParams.pixelSizeNm = 97.5;
simParams.emissionWavelengthNm = 561;
simParams.NA = 1.42;
simParams.noiseLevel = 0.05;
simParams.phaseErrorStd = 0.01;
simParams.randomSeed = 7;
simParams.useOTF = true;
[simParams.ksPixel, simParams.ktPixel] = expansionFactorToCarrierPixels( ...
    simParams.expansionFactor, simParams.imageSize, simParams.pixelSizeNm, ...
    simParams.emissionWavelengthNm, simParams.NA);

[rawStack, truth] = simulateLatticeSIMExperiment(objectImage, simParams);

reconParams = defaultLatticeSIMParams();
reconParams.pixelSizeNm = simParams.pixelSizeNm;
reconParams.emissionWavelengthNm = simParams.emissionWavelengthNm;
reconParams.NA = simParams.NA;
reconParams.modulationS = simParams.modulationS;
reconParams.modulationT = simParams.modulationT;
reconParams.normalizeFrames = true;

result = reconstructLatticeSIM(rawStack, reconParams);
normalizedRawStack = normalizeSIMFrames(rawStack, reconParams);
comparisonBands = separateLatticeBands(normalizedRawStack, reconParams);
widefieldC0 = imresize(abs(comparisonBands.C0), 2);
widefieldD3D4D5 = imresize(mean(normalizedRawStack(:, :, 3:5), 3), 2);

outputDir = fullfile(rootDir, 'testpat_lattice_simulation_output');
diagnosticsDir = fullfile(outputDir, 'diagnostics');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(diagnosticsDir, 'dir')
    mkdir(diagnosticsDir);
end

imwrite(toUnitImage(objectImage), fullfile(outputDir, 'cropped_testpat.tif'), 'tif');
writeStackTiff(rawStack, fullfile(outputDir, 'simulated_raw_stack.tif'));
deleteIfExists(fullfile(outputDir, 'Wide-field.tif'));
imwrite(toUnitImage(widefieldC0), fullfile(outputDir, 'Wide-field-C0.tif'), 'tif');
imwrite(toUnitImage(widefieldD3D4D5), fullfile(outputDir, 'Wide-field-D3D4D5.tif'), 'tif');
imwrite(toUnitImage(result.SIM), fullfile(outputDir, 'Lattice-SIM.tif'), 'tif');

save(fullfile(diagnosticsDir, 'rawStack.mat'), 'rawStack');
save(fullfile(diagnosticsDir, 'result.mat'), 'result');
save(fullfile(diagnosticsDir, 'simulation_truth.mat'), 'truth');

fprintf('Cropped test image: %s\n', fullfile(outputDir, 'cropped_testpat.tif'));
fprintf('Simulated raw stack: %s\n', fullfile(outputDir, 'simulated_raw_stack.tif'));
fprintf('C0 wide-field comparison: %s\n', fullfile(outputDir, 'Wide-field-C0.tif'));
fprintf('D3-D4-D5 wide-field comparison: %s\n', fullfile(outputDir, 'Wide-field-D3D4D5.tif'));
fprintf('Lattice-SIM result: %s\n', fullfile(outputDir, 'Lattice-SIM.tif'));

function writeStackTiff(stack, outputPath)
if exist(outputPath, 'file') == 2
    delete(outputPath);
end
for idx = 1:size(stack, 3)
    frame = toUnitImage(stack(:, :, idx));
    if idx == 1
        imwrite(frame, outputPath, 'tif');
    else
        imwrite(frame, outputPath, 'tif', 'WriteMode', 'append');
    end
end
end

function image = toUnitImage(image)
image = double(image);
image = image - min(image(:));
maxValue = max(image(:));
if maxValue > 0
    image = image ./ maxValue;
end
image = im2uint16(image);
end

function deleteIfExists(path)
if exist(path, 'file') == 2
    delete(path);
end
end
