function [normalized, info] = normalizeSIMFrames(stack, params)
%NORMALIZESIMFRAMES Normalize each raw frame to the stack-wide mean intensity.

if nargin < 2
    params = defaultLatticeSIMParams();
end

validateLatticeSIMStack(stack);

if isfield(params, 'normalizeFrames') && ~params.normalizeFrames
    normalized = double(stack);
    info.frameMeans = squeeze(mean(mean(normalized, 1), 2));
    info.targetMean = mean(info.frameMeans);
    info.scaleFactors = ones(5, 1);
    return;
end

normalized = double(stack);
frameMeans = squeeze(mean(mean(normalized, 1), 2));
if any(~isfinite(frameMeans)) || any(frameMeans <= 0)
    error('LatticeSIM:InvalidFrameMean', 'All frames must have finite positive mean intensity for normalization.');
end

targetMean = mean(frameMeans);
scaleFactors = targetMean ./ frameMeans;
for idx = 1:5
    normalized(:, :, idx) = normalized(:, :, idx) .* scaleFactors(idx);
end

info.frameMeans = frameMeans;
info.targetMean = targetMean;
info.scaleFactors = scaleFactors;
end
