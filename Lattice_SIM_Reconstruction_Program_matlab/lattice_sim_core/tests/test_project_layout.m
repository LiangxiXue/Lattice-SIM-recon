function test_project_layout()
%TEST_PROJECT_LAYOUT Verify the Lattice-SIM core folder organization.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);

assert(exist(fullfile(coreDir, 'main'), 'dir') == 7);
assert(exist(fullfile(coreDir, 'functions'), 'dir') == 7);
assert(exist(fullfile(coreDir, 'tests'), 'dir') == 7);
assert(exist(fullfile(coreDir, 'simulation_tests'), 'dir') == 7);
assert(exist(fullfile(coreDir, 'debug_experiments'), 'dir') == 7);

assert(exist(fullfile(coreDir, 'main', 'run_lattice_sim_configurable_test.m'), 'file') == 2);
assert(exist(fullfile(coreDir, 'simulation_tests', 'simulate_testpat_lattice_sim.m'), 'file') == 2);
assert(exist(fullfile(coreDir, 'simulation_tests', 'testpat.tiff'), 'file') == 2);
assert(exist(fullfile(coreDir, 'debug_experiments', 'component_separation_debug'), 'dir') == 7);
assert(exist(fullfile(coreDir, 'debug_experiments', 'direct_spectrum_comparison'), 'dir') == 7);

assert(exist(fullfile(coreDir, 'run_lattice_sim_configurable_test.m'), 'file') == 0);
assert(exist(fullfile(coreDir, 'simulate_testpat_lattice_sim.m'), 'file') == 0);
assert(exist(fullfile(coreDir, 'testpat.tiff'), 'file') == 0);
assert(exist(fullfile(coreDir, 'component_separation_debug'), 'dir') == 0);
assert(exist(fullfile(coreDir, 'direct_spectrum_comparison'), 'dir') == 0);
end
