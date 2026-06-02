function test_readCroppedLatticeSIMFrames()
%TEST_READCROPPEDLATTICESIMFRAMES Verify center ROI loading for five frames.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

dataDir = tempname;
mkdir(dataDir);
cleanup = onCleanup(@() removeDirectory(dataDir));

framePaths = cell(5, 1);
baseFrame = reshape(uint16(1:80), 8, 10);
for idx = 1:5
    framePaths{idx} = fullfile(dataDir, sprintf('%d.tiff', idx));
    imwrite(baseFrame + uint16(idx * 1000), framePaths{idx}, 'tiff');
end

[stack, metadata] = readCroppedLatticeSIMFrames(framePaths, [4, 6], []);

assert(isequal(size(stack), [4, 6, 5]));
assert(isequal(stack(:, :, 1), double(baseFrame(3:6, 3:8) + uint16(1000))));
assert(strcmp(metadata.sourceType, 'croppedTiffFiles'));
assert(isequal(metadata.originalFrameSize, [8, 10]));
assert(isequal(metadata.cropSizePixels, [4, 6]));
assert(isequal(metadata.rowRange, 3:6));
assert(isequal(metadata.colRange, 3:8));
assert(metadata.numFrames == 5);
end

function removeDirectory(path)
if exist(path, 'dir')
    rmdir(path, 's');
end
end
