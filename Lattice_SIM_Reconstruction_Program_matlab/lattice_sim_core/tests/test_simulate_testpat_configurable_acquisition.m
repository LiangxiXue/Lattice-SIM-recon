function test_simulate_testpat_configurable_acquisition()
%TEST_SIMULATE_TESTPAT_CONFIGURABLE_ACQUISITION Verify testpat script config plumbing.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
scriptPath = fullfile(coreDir, 'simulation_tests', 'simulate_testpat_lattice_sim.m');
scriptText = fileread(scriptPath);

configStart = strfind(scriptText, '%% Editable simulation and reconstruction parameters.');
firstRead = strfind(scriptText, 'objectImageFull = double(imread(inputPath));');
assert(~isempty(configStart));
assert(~isempty(firstRead));
assert(configStart(1) < firstRead(1));

editableBlock = scriptText(configStart(1):firstRead(1));
assert(contains(editableBlock, 'cropSize = 256;'));
assert(contains(editableBlock, 'frequencyMode = "expansion-factor-angle";'));
assert(contains(editableBlock, 'expansionFactor = 1.6;'));
assert(contains(editableBlock, 'carrierAngleTarget = "s";'));
assert(contains(editableBlock, 'carrierAngleDeg = 18;'));
assert(contains(editableBlock, '[ksPixel, ktPixel] = carrierPixelsFromExpansionAngle('));
assert(contains(editableBlock, 'acquisitionPhasePairs = ['));
assert(contains(editableBlock, 'reconstructionPhaseMatrix = [];'));
assert(contains(editableBlock, 'simParams.phasePairs = acquisitionPhasePairs;'));
assert(contains(editableBlock, 'reconParams.phasePairs = acquisitionPhasePairs;'));
assert(contains(editableBlock, 'reconParams.phaseMatrix = makeLatticePhaseMatrix(reconParams);'));
end
