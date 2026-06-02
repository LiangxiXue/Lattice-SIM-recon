function output = run_lattice_sim_configurable_test(varargin)
%RUN_LATTICE_SIM_CONFIGURABLE_TEST Reconstruct the bundled five-frame test data.
%
% Edit the configuration block at the top of this file when testing
% microscope settings, acquisition phase settings, and reconstruction modes.

%% Runtime switches. These can also be overridden by name-value arguments.
parser = inputParser;
parser.addParameter('dryRun', false, @(v) islogical(v) || isnumeric(v));
parser.addParameter('showFigures', false, @(v) islogical(v) || isnumeric(v));
parser.parse(varargin{:});
dryRun = logical(parser.Results.dryRun);
showFigures = logical(parser.Results.showFigures);

%% Input and output configuration.
coreDir = fileparts(mfilename('fullpath'));
programDir = fileparts(coreDir);
addpath(fullfile(coreDir, 'functions'));
dataDir = fullfile(programDir, 'Lattice_SIM_test_Data');
outputDir = fullfile(coreDir, 'configurable_real_data_output');
framePaths = {
    fullfile(dataDir, '1.tiff')
    fullfile(dataDir, '2.tiff')
    fullfile(dataDir, '3.tiff')
    fullfile(dataDir, '4.tiff')
    fullfile(dataDir, '5.tiff')
};

%% Microscope and acquisition parameters.
pixelSizeNm = 19.5;
emissionWavelengthNm = 532;
objectiveNA = 1.42;
modulationS = 0.7;
modulationT = 0.7;

%% Captured five-frame phase model.
% phasePairs columns are [s-direction phase, t-direction phase], in radians.
% Leave phaseMatrix empty to derive the 5 x 5 demodulation matrix from
% phasePairs. Set phaseMatrix to a custom 5 x 5 matrix when the acquisition
% calibration provides the full model directly.
phasePairs = [
    0,       0
    2*pi/3, 0
    4*pi/3, 0
    0,       2*pi/3
    2*pi/3, 4*pi/3
];
phaseMatrix = [];

%% Reconstruction algorithm parameters and modes.
params = defaultLatticeSIMParams();
params.pixelSizeNm = pixelSizeNm;
params.emissionWavelengthNm = emissionWavelengthNm;
params.NA = objectiveNA;
params.modulationS = modulationS;
params.modulationT = modulationT;
params.phasePairs = phasePairs;
params.phaseMatrix = phaseMatrix;
params.normalizeFrames = true;
params.preprocessingMode = "hifi-rl-fft";
params.deconvolutionIterations = 5;
params.separationFunction = "separateLatticeBandsFrequency";
params.separationInputDomain = "frequency";
params.enableLatticeParameterEstimation = true;
params.latticeCarrierRefinementIterations = 3;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCorrelationOverlap = 0.15;
params.carrierSearchMode = "unconstrained";
params.carrierAxisToleranceDeg = 15;
params.wiener = 0.04;
params.apodizationEnabled = true;
params.apodizationMode = "radial-gaussian";
params.apodizationStrength = 0.4;
params.apodizationRadius = 0.5;
params.supportThreshold = 1e-4;
params.reliabilityThreshold = 1e-3;
params.otfClipThreshold = 1e-4;
params.otfTaperLow = 1e-4;
params.otfTaperHigh = 5e-3;
params.otfAttenuationEnabled = true;
params.otfAttenuationStrength = 0.15;
params.otfAttenuationFwhm = 0.25;
params.notchScale = 0.5;
params.outputScaleMode = "none";
params = defaultLatticeSIMParams(params);

displayParams = params;
if isempty(displayParams.phaseMatrix)
    displayParams.phaseMatrix = makeLatticePhaseMatrix(displayParams);
end

config.input.mode = 'realFiles';
config.input.dataDir = dataDir;
config.input.framePaths = framePaths;
config.outputDir = outputDir;
config.params = displayParams;
config.showFigures = showFigures;

if dryRun
    output = config;
    return;
end

result = reconstructLatticeSIM(framePaths, params);
saved = saveLatticeSIMResult(result, outputDir);

if showFigures
    figure('Name', 'Wide-field');
    imagesc(result.WF);
    axis image off;
    colormap gray;
    colorbar;
    title('Wide-field');

    figure('Name', 'Lattice-SIM');
    imagesc(result.SIM);
    axis image off;
    colormap gray;
    colorbar;
    title('Lattice-SIM');
end

output.config = config;
output.result = result;
output.saved = saved;

fprintf('Wide-field result: %s\n', saved.widefieldPath);
fprintf('Lattice-SIM result: %s\n', saved.simPath);
fprintf('Diagnostics: %s\n', saved.diagnosticsDir);
end
