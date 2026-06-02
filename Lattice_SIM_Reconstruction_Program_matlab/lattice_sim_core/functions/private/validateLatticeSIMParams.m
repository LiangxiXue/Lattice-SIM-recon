function validateLatticeSIMParams(params, requirePhysical)
%VALIDATELATTICESIMPARAMS Validate user-facing reconstruction parameters.

if nargin < 2
    requirePhysical = false;
end

if params.modulationS <= 0 || params.modulationT <= 0
    error('LatticeSIM:InvalidModulation', 'Modulation factors must be positive.');
end
if params.wiener < 0
    error('LatticeSIM:InvalidWiener', 'Wiener parameter must be non-negative.');
end
if ~any(strcmp(char(params.preprocessingMode), {'hifi-rl-fft'}))
    error('LatticeSIM:InvalidPreprocessingMode', ...
        'Preprocessing mode must be "hifi-rl-fft".');
end
if params.deconvolutionIterations < 0 || fix(params.deconvolutionIterations) ~= params.deconvolutionIterations
    error('LatticeSIM:InvalidDeconvolutionIterations', ...
        'Deconvolution iterations must be a non-negative integer.');
end
if ~any(strcmp(char(params.separationFunction), {'separateLatticeBandsFrequency'}))
    error('LatticeSIM:InvalidSeparationFunction', ...
        'Separation function must be "separateLatticeBandsFrequency".');
end
if ~any(strcmp(char(params.separationInputDomain), {'frequency'}))
    error('LatticeSIM:InvalidSeparationInputDomain', ...
        'Separation input domain must be "frequency".');
end
if ~isfield(params, 'phasePairs') || ~isnumeric(params.phasePairs) || ...
        ~isequal(size(params.phasePairs), [5, 2]) || any(~isfinite(params.phasePairs(:)))
    error('LatticeSIM:InvalidPhasePairs', 'phasePairs must be a finite 5 x 2 numeric matrix.');
end
if isfield(params, 'phaseMatrix') && ~isempty(params.phaseMatrix) && ...
        (~isnumeric(params.phaseMatrix) || ~isequal(size(params.phaseMatrix), [5, 5]) || ...
        any(~isfinite(params.phaseMatrix(:))))
    error('LatticeSIM:InvalidPhaseMatrix', 'phaseMatrix must be empty or a finite 5 x 5 numeric matrix.');
end
if params.nrBands ~= 3
    error('LatticeSIM:InvalidBandCount', 'Lattice-SIM core expects nrBands = 3 for five-frame separation.');
end
if numel(params.bandFactors) ~= params.nrBands
    error('LatticeSIM:InvalidBandFactors', 'bandFactors must contain one value per band.');
end
if ~isnumeric(params.phaseOffsetS) || ~isscalar(params.phaseOffsetS)
    error('LatticeSIM:InvalidPhaseOffset', 'phaseOffsetS must be a scalar.');
end
if ~isnumeric(params.phaseOffsetT) || ~isscalar(params.phaseOffsetT)
    error('LatticeSIM:InvalidPhaseOffset', 'phaseOffsetT must be a scalar.');
end
if params.estimatedModulationS <= 0 || params.estimatedModulationT <= 0
    error('LatticeSIM:InvalidEstimatedModulation', ...
        'Estimated modulation values must be positive.');
end
if ~islogical(params.enableLatticeParameterEstimation) && ...
        ~(isnumeric(params.enableLatticeParameterEstimation) && isscalar(params.enableLatticeParameterEstimation))
    error('LatticeSIM:InvalidParameterEstimationFlag', ...
        'enableLatticeParameterEstimation must be a scalar logical or numeric value.');
end
if params.latticeCarrierRefinementIterations < 0 || ...
        fix(params.latticeCarrierRefinementIterations) ~= params.latticeCarrierRefinementIterations
    error('LatticeSIM:InvalidCarrierRefinementIterations', ...
        'latticeCarrierRefinementIterations must be a non-negative integer.');
end
if params.latticeCarrierSearchStepPixels <= 0
    error('LatticeSIM:InvalidCarrierSearchStep', ...
        'latticeCarrierSearchStepPixels must be positive.');
end
if params.latticeCorrelationOverlap <= 0 || params.latticeCorrelationOverlap >= 1
    error('LatticeSIM:InvalidCorrelationOverlap', ...
        'latticeCorrelationOverlap must be in the open interval (0, 1).');
end
if params.supportThreshold < 0 || params.supportThreshold > 1
    error('LatticeSIM:InvalidSupportThreshold', ...
        'Support threshold must be in the range [0, 1].');
end
if params.reliabilityThreshold < 0 || params.reliabilityThreshold > 1
    error('LatticeSIM:InvalidReliabilityThreshold', ...
        'Reliability threshold must be in the range [0, 1].');
end
if params.otfClipThreshold < 0 || params.otfClipThreshold > 1
    error('LatticeSIM:InvalidOtfClipThreshold', ...
        'OTF clip threshold must be in the range [0, 1].');
end
if params.otfTaperLow <= 0 || params.otfTaperLow > 1
    error('LatticeSIM:InvalidOtfTaperLow', ...
        'OTF taper low threshold must be in the range (0, 1].');
end
if params.otfTaperHigh <= 0 || params.otfTaperHigh > 1
    error('LatticeSIM:InvalidOtfTaperHigh', ...
        'OTF taper high threshold must be in the range (0, 1].');
end
if params.otfTaperLow >= params.otfTaperHigh
    error('LatticeSIM:InvalidOtfTaperRange', ...
        'otfTaperLow must be smaller than otfTaperHigh.');
end
if params.otfAttenuationStrength < 0 || params.otfAttenuationStrength >= 1
    error('LatticeSIM:InvalidOtfAttenuationStrength', ...
        'OTF attenuation strength must be in the range [0, 1).');
end
if params.otfAttenuationFwhm <= 0
    error('LatticeSIM:InvalidOtfAttenuationFwhm', ...
        'OTF attenuation FWHM must be positive.');
end
if params.apodizationRadius <= 0
    error('LatticeSIM:InvalidApodizationRadius', ...
        'Apodization radius must be positive.');
end
if ~any(strcmp(char(params.apodizationMode), {'support-distance', 'radial-gaussian'}))
    error('LatticeSIM:InvalidApodizationMode', ...
        'Apodization mode must be "support-distance" or "radial-gaussian".');
end
if ~any(strcmp(char(params.carrierSearchMode), {'axis-aligned', 'unconstrained'}))
    error('LatticeSIM:InvalidCarrierSearchMode', ...
        'Carrier search mode must be "axis-aligned" or "unconstrained".');
end
if params.carrierAxisToleranceDeg <= 0 || params.carrierAxisToleranceDeg > 90
    error('LatticeSIM:InvalidCarrierAxisTolerance', ...
        'Carrier axis tolerance must be in the range (0, 90] degrees.');
end
if requirePhysical
    required = {'pixelSizeNm', 'emissionWavelengthNm', 'NA'};
    for idx = 1:numel(required)
        name = required{idx};
        if ~isfield(params, name) || isempty(params.(name)) || ...
                ~isnumeric(params.(name)) || ~isscalar(params.(name)) || params.(name) <= 0
            error('LatticeSIM:MissingPhysicalParameter', ...
                'params.%s must be a positive scalar for OTF-based reconstruction.', name);
        end
    end
end
end
