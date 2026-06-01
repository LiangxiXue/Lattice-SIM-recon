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
