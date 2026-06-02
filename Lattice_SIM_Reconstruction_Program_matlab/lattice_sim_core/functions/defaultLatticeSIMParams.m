function params = defaultLatticeSIMParams(userParams)
%DEFAULTLATTICESIMPARAMS Return default parameters for 2D five-frame Lattice-SIM.

if nargin < 1 || isempty(userParams)
    userParams = struct();
end

params.pixelSizeNm = [];
params.emissionWavelengthNm = [];
params.NA = [];
params.modulationS = 0.7;
params.modulationT = 0.7;
params.wiener = 0.04;
params.apodizationEnabled = true;
params.apodizationStrength = 0.4;
params.apodizationMode = "radial-gaussian";
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
params.normalizeFrames = false;
params.preprocessingMode = "hifi-rl-fft";
params.deconvolutionIterations = 5;
params.separationFunction = "separateLatticeBandsFrequency";
params.separationInputDomain = "frequency";
params.phaOff = 0;
params.nrBands = 3;
params.bandFactors = ones(1, 3);
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;
params.enableLatticeParameterEstimation = true;
params.latticeCarrierRefinementIterations = 3;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCorrelationOverlap = 0.15;
params.carrierMinRadiusPixels = [];
params.carrierPeakWindow = 1;
params.carrierWeakPeakRatio = 4;
params.carrierSearchMode = "axis-aligned";
params.carrierAxisToleranceDeg = 15;

if ~isstruct(userParams)
    error('LatticeSIM:InvalidParams', 'Parameters must be a struct.');
end

names = fieldnames(userParams);
for idx = 1:numel(names)
    params.(names{idx}) = userParams.(names{idx});
end

validateLatticeSIMParams(params, false);
end
