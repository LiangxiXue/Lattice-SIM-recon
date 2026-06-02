function [stack, metadata] = readCroppedLatticeSIMFrames(framePaths, cropSizePixels, cropCenterPixels)
%READCROPPEDLATTICESIMFRAMES Read a center ROI from five separate TIFF files.

paths = cellstr(framePaths(:));
if numel(paths) ~= 5
    error('LatticeSIM:InvalidFrameCount', 'Expected exactly five TIFF files, got %d.', numel(paths));
end
if nargin < 2 || isempty(cropSizePixels)
    error('LatticeSIM:InvalidCropSize', 'cropSizePixels must be [height, width].');
end
if nargin < 3
    cropCenterPixels = [];
end

cropSizePixels = double(cropSizePixels(:)');
if numel(cropSizePixels) ~= 2 || any(~isfinite(cropSizePixels)) || ...
        any(cropSizePixels <= 0) || any(fix(cropSizePixels) ~= cropSizePixels)
    error('LatticeSIM:InvalidCropSize', 'cropSizePixels must be two positive integer values.');
end

for idx = 1:5
    if exist(paths{idx}, 'file') ~= 2
        error('LatticeSIM:FileNotFound', 'Input TIFF file does not exist: %s', paths{idx});
    end
    info = imfinfo(paths{idx});
    if idx == 1
        originalFrameSize = [info(1).Height, info(1).Width];
        originalClass = class(imread(paths{idx}, 'PixelRegion', {[1, 1], [1, 1]}));
        [rowRange, colRange, cropCenterPixels] = cropRanges(originalFrameSize, cropSizePixels, cropCenterPixels);
        stack = zeros(cropSizePixels(1), cropSizePixels(2), 5);
    elseif ~isequal([info(1).Height, info(1).Width], originalFrameSize)
        error('LatticeSIM:InconsistentFrameSize', 'All TIFF files must have identical dimensions.');
    end

    frame = imread(paths{idx}, 'PixelRegion', ...
        {[rowRange(1), rowRange(end)], [colRange(1), colRange(end)]});
    stack(:, :, idx) = double(frame);
end

metadata.sourceType = 'croppedTiffFiles';
metadata.originalClass = originalClass;
metadata.originalFrameSize = originalFrameSize;
metadata.frameSize = cropSizePixels;
metadata.cropSizePixels = cropSizePixels;
metadata.cropCenterPixels = cropCenterPixels;
metadata.rowRange = rowRange;
metadata.colRange = colRange;
metadata.numFrames = 5;
end

function [rowRange, colRange, centerPixels] = cropRanges(frameSize, cropSizePixels, centerPixels)
if isempty(centerPixels)
    rowStart = floor((frameSize(1) - cropSizePixels(1)) / 2) + 1;
    colStart = floor((frameSize(2) - cropSizePixels(2)) / 2) + 1;
    centerPixels = [rowStart + (cropSizePixels(1) - 1) / 2, ...
        colStart + (cropSizePixels(2) - 1) / 2];
else
    centerPixels = double(centerPixels(:)');
    if numel(centerPixels) ~= 2 || any(~isfinite(centerPixels))
        error('LatticeSIM:InvalidCropCenter', 'cropCenterPixels must be empty or [row, col].');
    end
    rowStart = round(centerPixels(1) - (cropSizePixels(1) - 1) / 2);
    colStart = round(centerPixels(2) - (cropSizePixels(2) - 1) / 2);
end

rowEnd = rowStart + cropSizePixels(1) - 1;
colEnd = colStart + cropSizePixels(2) - 1;
if rowStart < 1 || colStart < 1 || rowEnd > frameSize(1) || colEnd > frameSize(2)
    error('LatticeSIM:CropOutsideFrame', ...
        'Requested crop [%d x %d] does not fit inside frame [%d x %d].', ...
        cropSizePixels(1), cropSizePixels(2), frameSize(1), frameSize(2));
end

rowRange = rowStart:rowEnd;
colRange = colStart:colEnd;
end
