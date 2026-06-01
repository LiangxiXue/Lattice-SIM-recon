function test_saveLatticeSIMResult()
%TEST_SAVELATTICESIMRESULT Verify saving stays outside the core algorithm.

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));

result.WF = reshape(1:16, [4, 4]);
result.SIM = reshape(16:-1:1, [4, 4]);
result.params = defaultLatticeSIMParams();
result.diagnostics.warnings = {};

saved = saveLatticeSIMResult(result, tmpDir);

assert(exist(fullfile(tmpDir, 'Wide-field.tif'), 'file') == 2);
assert(exist(fullfile(tmpDir, 'Lattice-SIM.tif'), 'file') == 2);
assert(exist(fullfile(tmpDir, 'diagnostics'), 'dir') == 7);
assert(exist(saved.resultMatPath, 'file') == 2);

fprintf('test_saveLatticeSIMResult passed.\n');
end
