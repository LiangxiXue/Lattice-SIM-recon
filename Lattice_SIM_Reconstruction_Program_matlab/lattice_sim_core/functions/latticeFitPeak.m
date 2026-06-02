function peak = latticeFitPeak(band0, band1, otfValues, carrierPixel, weightLimit, search, iterations)
%LATTICEFITPEAK Refine a Lattice carrier with HiFi-style correlation search.

if nargin < 5 || isempty(weightLimit)
    weightLimit = 0.15;
end
if nargin < 6 || isempty(search)
    search = 2.5;
end
if nargin < 7 || isempty(iterations)
    iterations = 3;
end

kx = carrierPixel(1);
ky = carrierPixel(2);
bestCorrelation = 0;
control = zeros(10, 10, max(1, iterations));

for iter = 1:iterations
    [b0, b1] = latticeCommonRegion(band0, band1, otfValues, [kx, ky], 0.15, weightLimit, true);
    b0 = FFT2D(b0, true);
    b1 = FFT2D(b1, true);

    denominator = sum(abs(b0).^2, 'all');
    if denominator <= eps
        error('LatticeSIM:DegenerateCorrelationReference', ...
            'C0 overlap region has zero energy.');
    end

    corr = zeros(10, 10);
    cmax = 0;
    cmin = inf;
    newKx = kx;
    newKy = ky;

    for yi = 1:10
        for xi = 1:10
            xpos = kx + (((xi - 1) - 4.5) / 4.5) * search;
            ypos = ky + (((yi - 1) - 4.5) / 4.5) * search;
            b1s = latticeFourierShift(b1, -xpos, -ypos);
            value = sum(b1s .* conj(b0), 'all') ./ denominator;
            corr(xi, yi) = value;
            if abs(value) > cmax
                cmax = abs(value);
                newKx = xpos;
                newKy = ypos;
                bestCorrelation = value;
            end
            if abs(value) < cmin
                cmin = abs(value);
            end
        end
    end

    if cmax > cmin
        control(:, :, iter) = (abs(corr).' - cmin) ./ (cmax - cmin);
    end

    kx = newKx;
    ky = newKy;
    search = search / 3;
end

peak.kx = kx;
peak.ky = ky;
peak.resPhase = angle(bestCorrelation);
peak.resMag = abs(bestCorrelation);
peak.correlation = bestCorrelation;
peak.control = control;
end
