function result = reconstructLatticeSIM(input, params)
%RECONSTRUCTLATTICESIM Reconstruct 2D five-frame Lattice-SIM data.

if nargin < 2
    params = struct();
end
params = defaultLatticeSIMParams(params);
validateLatticeSIMParams(params, true);

[rawStack, inputMetadata] = readLatticeSIMInput(input);
[stack, preprocessingInfo] = preprocessLatticeSIMFrames(rawStack, params);
params.preprocessingInfo = preprocessingInfo;

firstPassBands = separateLatticeBandsFrequency(stack, params);
if params.enableLatticeParameterEstimation
    latticeEstimate = estimateLatticeBandParameters(firstPassBands, params);
    params.phaseOffsetS = latticeEstimate.phaseOffsetS;
    params.phaseOffsetT = latticeEstimate.phaseOffsetT;
    params.estimatedModulationS = latticeEstimate.modulationS;
    params.estimatedModulationT = latticeEstimate.modulationT;
    params.modulationS = latticeEstimate.modulationS;
    params.modulationT = latticeEstimate.modulationT;
    bands = separateLatticeBandsFrequency(stack, params);
    carriers = latticeEstimate.carriers;
    carrierDiagnostics = latticeEstimate.diagnostics;
else
    bands = firstPassBands;
    [carriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
    latticeEstimate = struct();
end
otf = buildLatticeOTF(size(stack, 1), size(stack, 2), params);
[SIM, combineDiagnostics] = combineLatticeSpectrum(bands, carriers, otf, params);

result.WF = makeWidefieldFromRawFrames(rawStack);
result.SIM = SIM;
result.params = params;
result.diagnostics = carrierDiagnostics;
result.diagnostics.input = inputMetadata;
result.diagnostics.preprocessing = preprocessingInfo;
result.diagnostics.separation = bands.diagnostics;
result.diagnostics.latticeEstimation = latticeEstimate;
result.diagnostics.otf = rmfield(otf, {'values'});
result.diagnostics.combine = combineDiagnostics;
end

function [stack, info] = preprocessLatticeSIMFrames(rawStack, params)
validateLatticeSIMStack(rawStack);

rawStack = double(rawStack);
[h, w, frameCount] = size(rawStack);
otf = buildLatticeOTF(h, w, params);
psf = abs(otf2psf(ifftshift(otf.values)));

windowedStack = importImages(rawStack);
deconvolvedStack = deconvlucy(windowedStack, psf, params.deconvolutionIterations);

stack = zeros(h, w, frameCount);
for idx = 1:frameCount
    stack(:, :, idx) = FFT2D(deconvolvedStack(:, :, idx), false);
end

info.method = char(params.preprocessingMode);
info.inputDomain = 'space';
info.outputDomain = 'frequency';
info.deconvolutionIterations = params.deconvolutionIterations;
info.psfSize = size(psf);
info.frameMeans = squeeze(mean(mean(rawStack, 1), 2));
end

function image = makeWidefieldFromRawFrames(rawStack)
image = mean(double(rawStack(:, :, 3:5)), 3);
image = imresize(image, 2);
end
