function [spectrum, apo] = applyLatticeApodization(spectrum, supportMask, params)
%APPLYLATTICEAPODIZATION Optionally taper the combined spectrum support.

apo = double(supportMask);
if ~isfield(params, 'apodizationEnabled') || ~params.apodizationEnabled
    spectrum = spectrum .* apo;
    return;
end

if isfield(params, 'apodizationMode') && strcmp(char(params.apodizationMode), 'radial-gaussian')
    [h, w] = size(spectrum);
    [x, y] = meshgrid(1:w, 1:h);
    center = [floor(h/2) + 1, floor(w/2) + 1];
    radius = hypot((x - center(2)) * 2 / w, (y - center(1)) * 2 / h);
    apo = exp(-0.5 * (radius ./ params.apodizationRadius .* sqrt(2*log(2))).^2);
    apo = apo .* double(supportMask);
elseif exist('bwdist', 'file') == 2
    dist = bwdist(~supportMask);
    if max(dist(:)) > 0
        apo = (dist ./ max(dist(:))) .^ params.apodizationStrength;
    end
end
spectrum = spectrum .* apo;
end
