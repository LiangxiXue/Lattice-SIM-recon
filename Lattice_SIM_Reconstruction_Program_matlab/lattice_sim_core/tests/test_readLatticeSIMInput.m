function test_readLatticeSIMInput()
%TEST_READLATTICESIMINPUT Verify supported input forms and validation errors.

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));

stack = reshape(uint16(1:8*7*5), [8, 7, 5]);

[loaded, meta] = readLatticeSIMInput(stack);
assert(isa(loaded, 'double'));
assert(isequal(size(loaded), [8, 7, 5]));
assert(strcmp(meta.sourceType, 'array'));
assert(strcmp(meta.originalClass, 'uint16'));

stackPath = fullfile(tmpDir, 'stack.tif');
writeTempTiffStack(stack, stackPath);
[loadedStack, metaStack] = readLatticeSIMInput(stackPath);
assert(isequal(size(loadedStack), [8, 7, 5]));
assert(strcmp(metaStack.sourceType, 'tiffStack'));

paths = cell(1, 5);
for idx = 1:5
    paths{idx} = fullfile(tmpDir, sprintf('frame_%d.tif', idx));
    imwrite(stack(:, :, idx), paths{idx}, 'tif');
end
[loadedFiles, metaFiles] = readLatticeSIMInput(paths);
assert(isequal(size(loadedFiles), [8, 7, 5]));
assert(strcmp(metaFiles.sourceType, 'tiffFiles'));

expectError(@() readLatticeSIMInput(stack(:, :, 1:4)), 'LatticeSIM:InvalidFrameCount');

badPaths = paths;
imwrite(uint16(ones(9, 7)), badPaths{5}, 'tif');
expectError(@() readLatticeSIMInput(badPaths), 'LatticeSIM:InconsistentFrameSize');

fprintf('test_readLatticeSIMInput passed.\n');
end
