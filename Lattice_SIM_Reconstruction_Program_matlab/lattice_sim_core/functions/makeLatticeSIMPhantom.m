function objectImage = makeLatticeSIMPhantom(imageSize)
%MAKELATTICESIMPHANTOM Create a deterministic synthetic specimen image.

if nargin < 1
    imageSize = [256, 256];
end
if ~isnumeric(imageSize) || numel(imageSize) ~= 2
    error('LatticeSIM:InvalidSimulationParameter', 'imageSize must be a two-element numeric vector.');
end

h = imageSize(1);
w = imageSize(2);
[x, y] = meshgrid(linspace(-1, 1, w), linspace(-1, 1, h));

objectImage = 0.10 + 0.08 * cos(2*pi*2.2*x) + 0.06 * sin(2*pi*1.7*y);

spots = [
    -0.45, -0.35, 0.08, 1.10
     0.28, -0.42, 0.06, 0.85
    -0.15,  0.20, 0.10, 0.75
     0.46,  0.24, 0.07, 0.95
     0.04,  0.55, 0.05, 0.70
];

for idx = 1:size(spots, 1)
    x0 = spots(idx, 1);
    y0 = spots(idx, 2);
    sigma = spots(idx, 3);
    amplitude = spots(idx, 4);
    objectImage = objectImage + amplitude * exp(-((x - x0).^2 + (y - y0).^2) / (2*sigma^2));
end

objectImage = objectImage - min(objectImage(:));
maxValue = max(objectImage(:));
if maxValue > 0
    objectImage = objectImage ./ maxValue;
end
objectImage = 100 * objectImage;
end
