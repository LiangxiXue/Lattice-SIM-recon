%RUN_LATTICE_SIM_CONFIGURABLE_TEST Reconstruct the bundled five-frame test data.
%
% This is a script. Edit the parameter block below, or define variables such
% as dryRun, pixelSizeNm, and phasePairs in the workspace before running this
% file. The script leaves config, result, saved, and output in the workspace.

%% Runtime switches.
if ~exist('dryRun', 'var') || isempty(dryRun)
    dryRun = false;
end
if ~exist('showFigures', 'var') || isempty(showFigures)
    showFigures = false;
end
dryRun = logical(dryRun);
showFigures = logical(showFigures);

%% Input and output configuration.
mainDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(mainDir);
programDir = fileparts(coreDir);
addpath(fullfile(coreDir, 'functions'));

if ~exist('dataDir', 'var') || isempty(dataDir)
    dataDir = fullfile(programDir, 'Lattice_SIM_test_Data');
end
if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = fullfile(coreDir, 'configurable_real_data_output');
end
if ~exist('framePaths', 'var') || isempty(framePaths)
    framePaths = {
        fullfile(dataDir, '1.tiff')
        fullfile(dataDir, '2.tiff')
        fullfile(dataDir, '3.tiff')
        fullfile(dataDir, '4.tiff')
        fullfile(dataDir, '5.tiff')
    };
end
if ~exist('cropEnabled', 'var') || isempty(cropEnabled)
    cropEnabled = true;
end
if ~exist('cropSizePixels', 'var') || isempty(cropSizePixels)
    cropSizePixels = [1024, 1024];
end
if ~exist('cropCenterPixels', 'var')
    cropCenterPixels = [];
end
cropEnabled = logical(cropEnabled);

%% Microscope and acquisition parameters.
if ~exist('pixelSizeNm', 'var') || isempty(pixelSizeNm)
    pixelSizeNm = 19.5;
end
if ~exist('emissionWavelengthNm', 'var') || isempty(emissionWavelengthNm)
    emissionWavelengthNm = 532;
end
if ~exist('objectiveNA', 'var') || isempty(objectiveNA)
    objectiveNA = 1.42;
end
if ~exist('modulationS', 'var') || isempty(modulationS)
    modulationS = 0.7;
end
if ~exist('modulationT', 'var') || isempty(modulationT)
    modulationT = 0.7;
end

%% Captured five-frame phase model.
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
params.enableLatticeParameterEstimation = false;
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
params = defaultLatticeSIMParams(params);

displayParams = params;
if isempty(displayParams.phaseMatrix)
    displayParams.phaseMatrix = makeLatticePhaseMatrix(displayParams);
end

config.input.mode = 'realFiles';
config.input.dataDir = dataDir;
config.input.framePaths = framePaths;
config.input.cropEnabled = cropEnabled;
config.input.cropSizePixels = cropSizePixels;
config.input.cropCenterPixels = cropCenterPixels;
config.outputDir = outputDir;
config.params = displayParams;
config.showFigures = showFigures;

if dryRun
    if exist('result', 'var')
        clear result;
    end
    if exist('saved', 'var')
        clear saved;
    end
    output = config;
else
    if cropEnabled
        [reconstructionInput, cropMetadata] = readCroppedLatticeSIMFrames( ...
            framePaths, cropSizePixels, cropCenterPixels);
        config.input.cropMetadata = cropMetadata;
    else
        reconstructionInput = framePaths;
    end

    result = reconstructLatticeSIM(reconstructionInput, params);
    if cropEnabled
        result.diagnostics.input.crop = cropMetadata;
    end
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
