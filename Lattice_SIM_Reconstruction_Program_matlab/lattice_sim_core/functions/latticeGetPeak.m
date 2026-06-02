function ret = latticeGetPeak(band0, band1, otfValues, carrierPixel, weightLimit)
%LATTICEGETPEAK HiFi-style complex C0/sideband correlation for Lattice bands.

if nargin < 5 || isempty(weightLimit)
    weightLimit = 0.15;
end

dist = 0.15;
[b0, b1] = latticeCommonRegion(band0, band1, otfValues, carrierPixel, dist, weightLimit, true);

b0 = FFT2D(b0, true);
b1 = FFT2D(b1, true);
b1 = latticeFourierShift(b1, -carrierPixel(1), -carrierPixel(2));
b1 = b1 .* conj(b0);

denominator = sum(abs(b0).^2, 'all');
if denominator <= eps
    error('LatticeSIM:DegenerateCorrelationReference', ...
        'C0 overlap region has zero energy.');
end
ret = sum(b1, 'all') ./ denominator;
end
