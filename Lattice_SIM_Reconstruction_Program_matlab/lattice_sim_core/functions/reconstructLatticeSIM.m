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

if ~params.enableLatticeParameterEstimation
    params.estimatedModulationS = params.modulationS;
    params.estimatedModulationT = params.modulationT;
end
firstPassBands = separateLatticeBandsFrequency(stack, params);
if params.enableLatticeParameterEstimation
    latticeEstimate = estimateLatticeBandParameters(firstPassBands, params);
    rawModulationS = latticeEstimate.modulationS;
    rawModulationT = latticeEstimate.modulationT;
    [protectedModulationS, protectionAppliedS] = protectLatticeModulation(rawModulationS, params);
    [protectedModulationT, protectionAppliedT] = protectLatticeModulation(rawModulationT, params);
    params.phaseOffsetS = latticeEstimate.phaseOffsetS;
    params.phaseOffsetT = latticeEstimate.phaseOffsetT;
    params.estimatedModulationS = protectedModulationS;
    params.estimatedModulationT = protectedModulationT;
    params.modulationS = protectedModulationS;
    params.modulationT = protectedModulationT;
    latticeEstimate.rawEstimatedModulationS = rawModulationS;
    latticeEstimate.rawEstimatedModulationT = rawModulationT;
    latticeEstimate.protectedModulationS = protectedModulationS;
    latticeEstimate.protectedModulationT = protectedModulationT;
    latticeEstimate.modulationProtectionAppliedS = protectionAppliedS;
    latticeEstimate.modulationProtectionAppliedT = protectionAppliedT;
    latticeEstimate.diagnostics.rawEstimatedModulationS = rawModulationS;
    latticeEstimate.diagnostics.rawEstimatedModulationT = rawModulationT;
    latticeEstimate.diagnostics.protectedModulationS = protectedModulationS;
    latticeEstimate.diagnostics.protectedModulationT = protectedModulationT;
    latticeEstimate.diagnostics.modulationProtectionAppliedS = protectionAppliedS;
    latticeEstimate.diagnostics.modulationProtectionAppliedT = protectionAppliedT;
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

function [modulation, protectionApplied] = protectLatticeModulation(rawModulation, params)
modulation = rawModulation;
protectionApplied = false;
if isfield(params, 'modulationProtectionEnabled') && params.modulationProtectionEnabled && ...
        rawModulation < params.modulationMinReliable
    modulation = params.modulationFallback;
    protectionApplied = true;
end
end
