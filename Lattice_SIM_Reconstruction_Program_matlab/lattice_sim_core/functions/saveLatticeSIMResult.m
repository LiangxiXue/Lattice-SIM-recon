function saved = saveLatticeSIMResult(result, outputDir)
%SAVELATTICESIMRESULT Save reconstruction outputs without modifying the result.

if ~isstruct(result) || ~isfield(result, 'WF') || ~isfield(result, 'SIM')
    error('LatticeSIM:InvalidResult', 'Result must contain WF and SIM fields.');
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diagnosticsDir = fullfile(outputDir, 'diagnostics');
if ~exist(diagnosticsDir, 'dir')
    mkdir(diagnosticsDir);
end

wfPath = fullfile(outputDir, 'Wide-field.tif');
simPath = fullfile(outputDir, 'Lattice-SIM.tif');
imwrite(toUint16Image(result.WF), wfPath, 'tif');
imwrite(toUint16Image(result.SIM), simPath, 'tif');

resultMatPath = fullfile(diagnosticsDir, 'result.mat');
save(resultMatPath, 'result', '-v7.3');

saved.widefieldPath = wfPath;
saved.simPath = simPath;
saved.diagnosticsDir = diagnosticsDir;
saved.resultMatPath = resultMatPath;
end

function image = toUint16Image(image)
image = double(image);
image = image - min(image(:));
maxValue = max(image(:));
if maxValue > 0
    image = image ./ maxValue;
end
image = uint16(round(image .* double(intmax('uint16'))));
end
