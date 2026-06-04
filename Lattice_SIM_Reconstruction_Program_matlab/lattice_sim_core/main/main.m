%RUN_LATTICE_SIM_CONFIGURABLE_TEST Reconstruct the bundled five-frame test data.
%
% This is a script. Edit the configuration block below, or define variables
% such as showFigures, stackPath, microscopeParams, and phasePairs in the
% workspace before running this file. The script leaves config, result, saved,
% and output in the workspace.

close all;
clear all;
clc;
%% Runtime switches.
if ~exist('showFigures', 'var') || isempty(showFigures)
    showFigures = false;
end
showFigures = logical(showFigures);

%% Input and output configuration.
mainDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(mainDir);
programDir = fileparts(coreDir);
addpath(fullfile(coreDir, 'functions'));

if exist('stackPath', 'var') && ~isempty(stackPath)
    stackPath = char(stackPath);
    if exist(stackPath, 'file') ~= 2
        warning('LatticeSIM:MissingPresetStackPath', ...
            'Preset stackPath does not exist and will be reselected: %s', stackPath);
        stackPath = '';
    end
end

if ~exist('stackPath', 'var') || isempty(stackPath)
    [stackFileName, stackDir] = uigetfile( ...
        {'*.tif;*.tiff', 'TIFF stack (*.tif, *.tiff)'}, ...
        'Select five-frame Lattice-SIM TIFF stack', ...
        programDir);
    if isequal(stackFileName, 0)
        error('LatticeSIM:InputSelectionCancelled', ...
            'No five-frame Lattice-SIM TIFF stack was selected.');
    end
    stackPath = fullfile(stackDir, stackFileName);
end
stackPath = char(stackPath);

if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = fullfile(coreDir, 'output');
end

%% Microscope and captured five-frame phase model.
if ~exist('microscopeParams', 'var') || isempty(microscopeParams)
    microscopeParams = struct();
end
microscopeDefaults = struct( ...
    'pixelSizeNm', 97.5, ...
    'emissionWavelengthNm', 561, ...
    'objectiveNA', 1.42);
microscopeNames = fieldnames(microscopeDefaults);
for microscopeIdx = 1:numel(microscopeNames)
    name = microscopeNames{microscopeIdx};
    if ~isfield(microscopeParams, name) || isempty(microscopeParams.(name))
        microscopeParams.(name) = microscopeDefaults.(name);
    end
end

% phasePairs columns are [s-direction phase, t-direction phase], in radians.
% Leave phaseMatrix empty to derive the 5 x 5 demodulation matrix from
% phasePairs. Set phaseMatrix to a custom 5 x 5 matrix when the acquisition
% calibration provides the full model directly.
if ~exist('phasePairs', 'var') || isempty(phasePairs)
    phasePairs = [
        0,       0
        2*pi/3, 0
        4*pi/3, 0
        0,       2*pi/3
        2*pi/3, 4*pi/3
    ];
end
if ~exist('phaseMatrix', 'var')
    phaseMatrix = [];
end

%% Reconstruction algorithm parameters and modes.
params = defaultLatticeSIMParams();
params.pixelSizeNm = microscopeParams.pixelSizeNm;
params.emissionWavelengthNm = microscopeParams.emissionWavelengthNm;
params.NA = microscopeParams.objectiveNA;
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

if exist('paramOverrides', 'var') && ~isempty(paramOverrides)
    overrideNames = fieldnames(paramOverrides);
    for overrideIdx = 1:numel(overrideNames)
        params.(overrideNames{overrideIdx}) = paramOverrides.(overrideNames{overrideIdx});
    end
end
params.enableLatticeParameterEstimation = true;
params.carrierSearchMode = "unconstrained";
params = defaultLatticeSIMParams(params);
params.enableLatticeParameterEstimation = true;
params.carrierSearchMode = "unconstrained";

displayParams = params;
if isempty(displayParams.phaseMatrix)
    displayParams.phaseMatrix = makeLatticePhaseMatrix(displayParams);
end

config.input.mode = 'tiffStack';
config.input.stackPath = stackPath;
config.outputDir = outputDir;
config.params = displayParams;
config.showFigures = showFigures;

result = reconstructLatticeSIM(stackPath, params);
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
