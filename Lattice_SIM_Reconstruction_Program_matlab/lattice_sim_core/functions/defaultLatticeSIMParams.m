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
params.fusionMode = "single-step";
params.wiener = 0.04;
params.wienerW1 = 1.2;
params.wienerW2 = 0.1;
params.hifiDenominatorScaleW1 = 1.0;
params.hifiDenominatorScaleW2 = 1.0;
params.hifiCenterDenominatorWeight = 0.5;
params.hifiW1SidebandAttenuationScale = 1.0;
params.hifiW2CenterAttenuationScale = 1.05;
params.hifiW2SidebandAttenuationScale = 1.15;
params.apodizationEnabled = true;
params.apodizationStrength = 0.4;
params.apodizationMode = "radial-gaussian";
params.apodizationRadius = 0.5;
params.supportThreshold = 1e-4;
params.reliabilityThreshold = 1e-3;
params.otfClipThreshold = 1e-4;
params.otfTaperLow = 1e-4;
params.otfTaperHigh = 5e-3;
params.hifiOtfA = 0.85;
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
params.phasePairs = [
    0,       0
    2*pi/3, 0
    4*pi/3, 0
    0,       2*pi/3
    2*pi/3, 4*pi/3
];
params.phaseMatrix = [];
params.phaOff = 0;
params.nrBands = 3;
params.bandFactors = ones(1, 3);
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;
params.modulationProtectionEnabled = true;
params.modulationMinReliable = 0.35;
params.modulationFallback = 0.7;
params.enableLatticeParameterEstimation = true;
params.latticeCarrierRefinementIterations = 3;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCorrelationOverlap = 0.15;
params.latticePhaseSearchCoarseStepRad = pi / 12;
params.latticePhaseSearchFineStepRad = pi / 72;
params.latticePhaseSearchFineRadiusRad = pi / 12;
params.latticeModulationSearchMin = 0.05;
params.latticeModulationSearchMax = 1.0;
params.latticeModulationSearchCoarseStep = 0.05;
params.latticeModulationSearchFineStep = 0.01;
params.latticeModulationSearchFineRadius = 0.05;
params.sidebandPairResidualWarningPixels = 1.0;
params.carrierMinRadiusPixels = [];
params.carrierPeakWindow = 1;
params.carrierWeakPeakRatio = 4;
params.carrierSearchMode = "unconstrained";
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
