function run_all_lattice_sim_tests()
%RUN_ALL_LATTICE_SIM_TESTS Run the lightweight Lattice-SIM test suite.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'functions'));
addpath(fullfile(rootDir, 'tests'));
addpath(fullfile(rootDir, 'tests', 'helpers'));

resultsDir = fullfile(rootDir, 'tests', 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
diaryPath = fullfile(resultsDir, 'run_all_lattice_sim_tests_latest.txt');
if exist(diaryPath, 'file')
    delete(diaryPath);
end
diary(diaryPath);
diaryCleanup = onCleanup(@() diary('off'));

test_readLatticeSIMInput();
test_normalizeSIMFrames();
test_separateLatticeBands();
test_expansionFactorToCarrierPixels();
test_simulateLatticeSIMExperiment();
test_estimateLatticeCarrier();
test_reconstructLatticeSIM_smoke();
test_saveLatticeSIMResult();

fprintf('All Lattice-SIM tests passed.\n');
end
