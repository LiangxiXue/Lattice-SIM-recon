function separate = separateBands(IrawFFT, phaOff, bands, fac)
%SEPARATEBANDS HiFi-SIM frequency-domain band separation.

phaPerBand = (bands * 2) - 1;
phases = zeros(1, phaPerBand);
for p = 1:phaPerBand
    phases(p) = (2 * pi * (p - 1)) / phaPerBand + phaOff;
end
separate = separateBands_final(IrawFFT, phases, bands, fac);
end

function separate = separateBands_final(IrawFFT, phases, bands, fac)
if fac == 0
    fac = zeros(1, bands);
    fac(1) = 1;
    for idx = 2:bands
        fac(idx) = 0.5;
    end
else
    for idx = 2:bands
        fac(idx) = fac(idx) * 0.5;
    end
end

comp = zeros(1, bands * 2 - 1);
comp(1) = 0;
for idx = 2:bands
    comp((idx - 1) * 2) = idx - 1;
    comp((idx - 1) * 2 + 1) = -(idx - 1);
end

compfac = zeros(1, bands * 2 - 1);
compfac(1) = fac(1);
for idx = 2:bands
    compfac((idx - 1) * 2) = fac(idx);
    compfac((idx - 1) * 2 + 1) = fac(idx);
end

W = exp(1i * phases' * comp);
for idx = 1:bands * 2 - 1
    W(idx, :) = W(idx, :) .* compfac;
end

frameCount = size(phases, 2);
siz = size(IrawFFT(:, :, 1));
S = reshape(IrawFFT, [prod(siz), frameCount]) * pinv(W)';
separate = reshape(S, [siz, bands * 2 - 1]);
end
