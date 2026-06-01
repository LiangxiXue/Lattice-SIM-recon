function [rawStack, bands] = makeSyntheticLatticeSIMStack(imageSize, varargin)
%MAKESYNTHETICLATTICESIMSTACK Build deterministic five-frame Lattice-SIM data.

parser = inputParser;
parser.addRequired('imageSize', @(v) isnumeric(v) && numel(v) == 2);
parser.addParameter('ksPixel', [6, 0], @(v) isnumeric(v) && numel(v) == 2);
parser.addParameter('ktPixel', [0, 7], @(v) isnumeric(v) && numel(v) == 2);
parser.addParameter('carrierAmplitude', 0.2, @(v) isnumeric(v) && isscalar(v));
parser.parse(imageSize, varargin{:});

h = imageSize(1);
w = imageSize(2);
[x, y] = meshgrid(0:w-1, 0:h-1);

base = 100 + 7 * cos(2*pi*x/w) + 5 * sin(2*pi*y/h);
bands.C0 = base;

amp = parser.Results.carrierAmplitude;
ks = parser.Results.ksPixel;
kt = parser.Results.ktPixel;

if amp == 0
    bands.CsPlus = 2 + 0.2 * cos(2*pi*x/w) + 1i * 0.1 * sin(2*pi*y/h);
    bands.CsMinus = conj(bands.CsPlus);
    bands.CtPlus = 4 + 0.3 * cos(2*pi*(x+y)/(w+h)) + 1i * 0.2 * sin(2*pi*(x-y)/(w+h));
    bands.CtMinus = conj(bands.CtPlus);
else
    phaseS = 2*pi*(ks(1) * x / w + ks(2) * y / h);
    phaseT = 2*pi*(kt(1) * x / w + kt(2) * y / h);
    bands.CsPlus = amp * base .* exp(1i * phaseS);
    bands.CsMinus = amp * base .* exp(-1i * phaseS);
    bands.CtPlus = amp * base .* exp(1i * phaseT);
    bands.CtMinus = amp * base .* exp(-1i * phaseT);
end

W = latticePhaseMatrix();
componentStack = cat(3, bands.C0, bands.CsPlus, bands.CsMinus, ...
    bands.CtPlus, bands.CtMinus);
rawStack = zeros(h, w, 5);
for frameIdx = 1:5
    frame = zeros(h, w);
    for bandIdx = 1:5
        frame = frame + W(frameIdx, bandIdx) .* componentStack(:, :, bandIdx);
    end
    rawStack(:, :, frameIdx) = real(frame);
end

end
