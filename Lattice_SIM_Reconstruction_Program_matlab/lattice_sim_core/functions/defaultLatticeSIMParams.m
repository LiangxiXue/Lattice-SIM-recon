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
params.apodizationEnabled = false;
params.apodizationStrength = 0.4;
params.notchScale = 0.5;
params.outputScaleMode = "none";
params.normalizeFrames = true;
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
