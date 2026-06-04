function otf = buildLatticeOTF(imageHeight, imageWidth, params)
%BUILDLATTICEOTF Build a replaceable theoretical incoherent OTF model.

validateLatticeSIMParams(params, true);

[fx, fy] = makeFrequencyGrid(imageHeight, imageWidth, params.pixelSizeNm);
cutoff = 2 * params.NA / params.emissionWavelengthNm;
rho = hypot(fx, fy) / cutoff;

otfValues = zeros(imageHeight, imageWidth);
inside = rho <= 1;
idealValues = (2 / pi) * (acos(rho(inside)) - rho(inside) .* sqrt(1 - rho(inside).^2));
otfValues(inside) = idealValues .* (params.hifiOtfA .^ rho(inside));
if params.otfAttenuationEnabled
    otfValues = otfValues .* hifiAttenuationMask(fx, fy, params);
end

otf.values = otfValues;
otf.supportMask = otfValues > 0;
otf.fxCyclesPerNm = fx;
otf.fyCyclesPerNm = fy;
otf.cutoffCyclesPerNm = cutoff;
otf.pixelSizeNm = params.pixelSizeNm;
otf.hifiOtfA = params.hifiOtfA;
otf.attenuationEnabled = params.otfAttenuationEnabled;
otf.attenuationStrength = params.otfAttenuationStrength;
otf.attenuationFwhm = params.otfAttenuationFwhm;
end

function attenuation = hifiAttenuationMask(fx, fy, params)
radius = hypot(fx, fy);
fwhmCyclesPerNm = params.otfAttenuationFwhm * 1e-3;
attenuation = 1 - params.otfAttenuationStrength .* ...
    exp(-(radius .^ 2) ./ ((0.5 * fwhmCyclesPerNm) ^ 2));
attenuation = min(max(attenuation, 0), 1);
end
