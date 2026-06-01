function spectrum = applyLatticeApodization(spectrum, supportMask, params)
%APPLYLATTICEAPODIZATION Optionally taper the combined spectrum support.

if ~isfield(params, 'apodizationEnabled') || ~params.apodizationEnabled
    return;
end

if exist('bwdist', 'file') == 2
    dist = bwdist(~supportMask);
    if max(dist(:)) > 0
        apo = (dist ./ max(dist(:))) .^ params.apodizationStrength;
        spectrum = spectrum .* apo;
    end
end
end
