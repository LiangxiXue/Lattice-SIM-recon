function [stack, metadata] = readLatticeSIMInput(input)
%READLATTICESIMINPUT Read supported Lattice-SIM input forms as H x W x 5 double.

if isnumeric(input)
    metadata.sourceType = 'array';
    metadata.originalClass = class(input);
    stack = double(input);
    validateLatticeSIMStack(stack);
    metadata.frameSize = size(stack(:, :, 1));
    metadata.numFrames = size(stack, 3);
    return;
end

if ischar(input) || (isstring(input) && isscalar(input))
    path = char(input);
    if exist(path, 'file') ~= 2
        error('LatticeSIM:FileNotFound', 'Input TIFF stack does not exist: %s', path);
    end
    info = imfinfo(path);
    if numel(info) ~= 5
        error('LatticeSIM:InvalidFrameCount', 'Expected exactly five TIFF pages, got %d.', numel(info));
    end
    first = imread(path, 1);
    stack = zeros(size(first, 1), size(first, 2), 5);
    stack(:, :, 1) = double(first);
    for idx = 2:5
        frame = imread(path, idx);
        if ~isequal(size(frame), size(first))
            error('LatticeSIM:InconsistentFrameSize', 'All TIFF pages must have identical dimensions.');
        end
        stack(:, :, idx) = double(frame);
    end
    metadata.sourceType = 'tiffStack';
    metadata.originalClass = class(first);
    metadata.frameSize = size(first);
    metadata.numFrames = 5;
    return;
end

if iscell(input) || isstring(input)
    paths = cellstr(input(:));
    if numel(paths) ~= 5
        error('LatticeSIM:InvalidFrameCount', 'Expected exactly five TIFF files, got %d.', numel(paths));
    end
    for idx = 1:5
        if exist(paths{idx}, 'file') ~= 2
            error('LatticeSIM:FileNotFound', 'Input TIFF file does not exist: %s', paths{idx});
        end
        frame = imread(paths{idx});
        if idx == 1
            firstSize = size(frame);
            stack = zeros(firstSize(1), firstSize(2), 5);
            originalClass = class(frame);
        elseif ~isequal(size(frame), firstSize)
            error('LatticeSIM:InconsistentFrameSize', 'All TIFF files must have identical dimensions.');
        end
        stack(:, :, idx) = double(frame);
    end
    metadata.sourceType = 'tiffFiles';
    metadata.originalClass = originalClass;
    metadata.frameSize = firstSize;
    metadata.numFrames = 5;
    return;
end

error('LatticeSIM:InvalidInput', 'Input must be a TIFF path, five TIFF paths, or an H x W x 5 numeric array.');
end
