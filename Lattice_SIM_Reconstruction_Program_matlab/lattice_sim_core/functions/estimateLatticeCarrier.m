function [carriers, diagnostics] = estimateLatticeCarrier(bands, params)
%ESTIMATELATTICECARRIER Estimate orthogonal Lattice carrier vectors.

if nargin < 2
    params = defaultLatticeSIMParams();
end

[h, w] = size(bands.C0);
if isempty(params.carrierMinRadiusPixels)
    minRadius = max(2, round(min(h, w) * params.notchScale / 15));
else
    minRadius = params.carrierMinRadiusPixels;
end

if strcmp(char(params.carrierSearchMode), 'axis-aligned')
    sAxis = 'x';
    tAxis = 'y';
else
    sAxis = 'none';
    tAxis = 'none';
end
bandDomain = getBandDomain(bands);

[ksPlus, mapSPlus, ratioSPlus] = findCarrierPeak(bands.CsPlus, minRadius, ...
    params.carrierPeakWindow, sAxis, params.carrierAxisToleranceDeg, bandDomain);
[ksMinus, mapSMinus, ratioSMinus] = findCarrierPeak(bands.CsMinus, minRadius, ...
    params.carrierPeakWindow, sAxis, params.carrierAxisToleranceDeg, bandDomain);
[ktPlus, mapTPlus, ratioTPlus] = findCarrierPeak(bands.CtPlus, minRadius, ...
    params.carrierPeakWindow, tAxis, params.carrierAxisToleranceDeg, bandDomain);
[ktMinus, mapTMinus, ratioTMinus] = findCarrierPeak(bands.CtMinus, minRadius, ...
    params.carrierPeakWindow, tAxis, params.carrierAxisToleranceDeg, bandDomain);

ks = (ksPlus - ksMinus) / 2;
kt = (ktPlus - ktMinus) / 2;

if norm(ks) < minRadius || norm(kt) < minRadius
    error('LatticeSIM:CarrierEstimationFailed', 'Carrier peak is too close to the DC region.');
end

carriers.ksPixel = ks;
carriers.ktPixel = kt;
carriers.ksRadPerPixel = [2*pi*ks(1)/w, 2*pi*ks(2)/h];
carriers.ktRadPerPixel = [2*pi*kt(1)/w, 2*pi*kt(2)/h];
carriers.searchMaps.s = max(mapSPlus, mapSMinus);
carriers.searchMaps.t = max(mapTPlus, mapTMinus);
carriers.peakStrengthS = min(ratioSPlus, ratioSMinus);
carriers.peakStrengthT = min(ratioTPlus, ratioTMinus);

angleS = atan2d(ks(2), ks(1));
angleT = atan2d(kt(2), kt(1));
orthogonalityErrorDeg = abs(90 - abs(wrapTo180Local(angleT - angleS)));

diagnostics.carrierS = carriers.ksPixel;
diagnostics.carrierT = carriers.ktPixel;
diagnostics.carrierMagnitudeS = norm(ks);
diagnostics.carrierMagnitudeT = norm(kt);
diagnostics.carrierAngleDeg = [angleS, angleT];
diagnostics.orthogonalityErrorDeg = orthogonalityErrorDeg;
diagnostics.carrierSearchMaps = carriers.searchMaps;
diagnostics.carrierSearchMode = params.carrierSearchMode;
diagnostics.bandDomain = bandDomain;
diagnostics.warnings = {};

if carriers.peakStrengthS < params.carrierWeakPeakRatio
    diagnostics.warnings{end + 1} = 'S carrier peak is weak or ambiguous.';
end
if carriers.peakStrengthT < params.carrierWeakPeakRatio
    diagnostics.warnings{end + 1} = 'T carrier peak is weak or ambiguous.';
end
if orthogonalityErrorDeg > 10
    diagnostics.warnings{end + 1} = 'Estimated carrier directions are not close to orthogonal.';
end
end

function domain = getBandDomain(bands)
if isfield(bands, 'domain')
    domain = char(bands.domain);
else
    domain = 'space';
end
end

function [carrierPixel, searchMap, peakRatio] = findCarrierPeak(component, minRadius, peakWindow, axisName, axisToleranceDeg, domain)
if strcmp(domain, 'frequency')
    F = abs(component);
else
    F = abs(fftshift(fft2(ifftshift(component))));
end
[h, w] = size(F);
[xGrid, yGrid] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
radius = hypot(xGrid, yGrid);
searchMap = F;
searchMap(radius < minRadius) = 0;
searchMap(~axisSearchMask(xGrid, yGrid, axisName, axisToleranceDeg)) = 0;
searchMap(~isfinite(searchMap)) = 0;

[peakValue, linearIdx] = max(searchMap(:));
if peakValue <= 0
    error('LatticeSIM:CarrierEstimationFailed', 'Could not find a nonzero carrier peak.');
end

background = median(searchMap(searchMap > 0));
if isempty(background) || background <= 0
    background = eps;
end
peakRatio = peakValue / background;
if peakRatio < 1.5
    error('LatticeSIM:CarrierEstimationFailed', 'Carrier peak is too weak or ambiguous.');
end

[row, col] = ind2sub(size(searchMap), linearIdx);
rowRange = max(1, row - peakWindow):min(h, row + peakWindow);
colRange = max(1, col - peakWindow):min(w, col + peakWindow);
localMap = searchMap(rowRange, colRange);
[localX, localY] = meshgrid(colRange - floor(w/2) - 1, rowRange - floor(h/2) - 1);
weightSum = sum(localMap(:));
carrierPixel = [sum(localX(:) .* localMap(:)) / weightSum, ...
    sum(localY(:) .* localMap(:)) / weightSum];
end

function mask = axisSearchMask(xGrid, yGrid, axisName, axisToleranceDeg)
if strcmp(axisName, 'none')
    mask = true(size(xGrid));
    return;
end

angleFromXAxis = atan2d(abs(yGrid), abs(xGrid));
if strcmp(axisName, 'x')
    mask = angleFromXAxis <= axisToleranceDeg;
elseif strcmp(axisName, 'y')
    mask = abs(angleFromXAxis - 90) <= axisToleranceDeg;
else
    error('LatticeSIM:InvalidCarrierSearchAxis', 'Unsupported carrier search axis: %s', axisName);
end
end

function angleDeg = wrapTo180Local(angleDeg)
angleDeg = mod(angleDeg + 180, 360) - 180;
end
