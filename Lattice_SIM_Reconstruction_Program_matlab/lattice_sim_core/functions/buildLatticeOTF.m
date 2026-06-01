function otf = buildLatticeOTF(imageHeight, imageWidth, params)
%BUILDLATTICEOTF Build a replaceable theoretical incoherent OTF model.

validateLatticeSIMParams(params, true);

[fx, fy] = makeFrequencyGrid(imageHeight, imageWidth, params.pixelSizeNm);
cutoff = 2 * params.NA / params.emissionWavelengthNm;
rho = hypot(fx, fy) / cutoff;

otfValues = zeros(imageHeight, imageWidth);
inside = rho <= 1;
otfValues(inside) = (2 / pi) * (acos(rho(inside)) - rho(inside) .* sqrt(1 - rho(inside).^2));

otf.values = otfValues;
otf.supportMask = otfValues > 0;
otf.fxCyclesPerNm = fx;
otf.fyCyclesPerNm = fy;
otf.cutoffCyclesPerNm = cutoff;
otf.pixelSizeNm = params.pixelSizeNm;
end
