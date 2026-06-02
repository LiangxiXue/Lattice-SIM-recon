function test_saveLatticeSIMResult()
%TEST_SAVELATTICESIMRESULT Verify reconstruction outputs are saved robustly.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() removeDirectory(outputDir));

result.WF = reshape(1:16, 4, 4);
result.SIM = reshape(16:-1:1, 4, 4);
result.params = struct();
result.diagnostics = struct();

saved = saveLatticeSIMResult(result, outputDir);

assert(exist(saved.widefieldPath, 'file') == 2);
assert(exist(saved.simPath, 'file') == 2);
assert(exist(saved.resultMatPath, 'file') == 2);

fid = fopen(saved.resultMatPath, 'r');
assert(fid > 0);
fileHeader = fread(fid, 64, '*char')';
fclose(fid);

assert(contains(fileHeader, 'MATLAB 7.3 MAT-file'), ...
    'Expected result.mat to use MAT-file v7.3/HDF5 format.');
end

function removeDirectory(path)
if exist(path, 'dir')
    rmdir(path, 's');
end
end
