function [SIM, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params)
%COMBINELATTICESPECTRUM Combine Lattice bands with OTF/Wiener damping.

[h, w] = size(bands.C0);
[x, y] = meshgrid(0:w-1, 0:h-1);

freq0 = fft2c(bands.C0);
num = freq0 .* conj(otf.values);
den = abs(otf.values).^2;

components = {
    bands.CsPlus,  carriers.ksRadPerPixel,  carriers.ksPixel, params.modulationS
    bands.CsMinus, -carriers.ksRadPerPixel, -carriers.ksPixel, params.modulationS
    bands.CtPlus,  carriers.ktRadPerPixel,  carriers.ktPixel, params.modulationT
    bands.CtMinus, -carriers.ktRadPerPixel, -carriers.ktPixel, params.modulationT
};

for idx = 1:size(components, 1)
    component = components{idx, 1};
    carrierRad = components{idx, 2};
    carrierPixel = components{idx, 3};
    sidebandAmplitude = components{idx, 4} / 2;
    phaseRamp = exp(-1i * (carrierRad(1) * x + carrierRad(2) * y));
    centered = component .* phaseRamp;
    freq = fft2c(centered);
    shiftedOtf = shiftOtfByCarrier(otf.values, carrierPixel);
    centerRow = floor(h/2) + 1;
    centerCol = floor(w/2) + 1;
    phaseReference = freq(centerRow, centerCol) / (freq0(centerRow, centerCol) + eps);
    freq = freq .* exp(-1i * angle(phaseReference));
    num = num + freq .* conj(shiftedOtf) ./ sidebandAmplitude;
    den = den + abs(shiftedOtf).^2;
end

combinedSpectrum = num ./ (den + params.wiener);
supportMask = den > 0;
combinedSpectrum = applyLatticeApodization(combinedSpectrum, supportMask, params);

SIM = imresize(abs(ifft2c(combinedSpectrum)), 2);

if any(~isfinite(SIM(:)))
    error('LatticeSIM:InvalidReconstruction', 'Reconstructed SIM image contains NaN or Inf values.');
end

diagnostics.simSpectrum = combinedSpectrum;
diagnostics.wienerDenominator = den;
diagnostics.supportMask = supportMask;
end

function shiftedOtf = shiftOtfByCarrier(otfValues, carrierPixel)
[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);
end
