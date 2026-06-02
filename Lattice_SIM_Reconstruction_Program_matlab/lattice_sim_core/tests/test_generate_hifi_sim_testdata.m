function test_generate_hifi_sim_testdata()
%TEST_GENERATE_HIFI_SIM_TESTDATA Verify the HiFi-SIM 2D test-data generator contract.

matlabRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
generatorDir = fullfile(matlabRoot, 'HiFi-SIM_v1.01', 'TestData');
addpath(generatorDir);

outputDir = fullfile(tempdir(), ['hifi_sim_testdata_', char(java.util.UUID.randomUUID())]);
cleanup = onCleanup(@() removeTestDir(outputDir));

objectImage = peaks(48);
objectImage = objectImage - min(objectImage(:));

[rawStack, truth, outputs] = generate_hifi_sim_testdata(outputDir, objectImage);

assert(isequal(size(rawStack), [48, 48, 9]));
assert(numel(imfinfo(outputs.rawStackTiff)) == 9);
assert(exist(outputs.truthMat, 'file') == 2);
assert(isequal(truth.directionIndex(:)', [1 1 1 2 2 2 3 3 3]));
assert(isequal(truth.phaseIndex(:)', [1 2 3 1 2 3 1 2 3]));
assert(all(isfinite(rawStack(:))));
assert(min(rawStack(:)) >= 0);
end

function removeTestDir(path)
if exist(path, 'dir') == 7
    rmdir(path, 's');
end
end
