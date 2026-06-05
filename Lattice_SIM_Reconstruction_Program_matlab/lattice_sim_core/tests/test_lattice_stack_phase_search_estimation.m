function test_lattice_stack_phase_search_estimation()
%TEST_LATTICE_STACK_PHASE_SEARCH_ESTIMATION Recover phase/modulation from raw frames.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;
params.carrierMinRadiusPixels = 2;
params.carrierWeakPeakRatio = 1.05;
params.latticePhaseSearchCoarseStepRad = pi / 6;
params.latticePhaseSearchFineStepRad = pi / 36;
params.latticePhaseSearchFineRadiusRad = pi / 6;

image = makeLatticeSIMPhantom([48, 48]);
ks = [7, 2];
kt = [-3, 8];
truePhaseS = 1.15;
truePhaseT = -0.85;
trueModulationS = 0.52;
trueModulationT = 0.38;

truthParams = params;
truthParams.phaseOffsetS = truePhaseS;
truthParams.phaseOffsetT = truePhaseT;
truthParams.estimatedModulationS = trueModulationS;
truthParams.estimatedModulationT = trueModulationT;

stackFFT = makeFrequencyStack(image, ks, kt, truthParams);
estimate = estimateLatticeStackParameters(stackFFT, params);

assert(norm(estimate.carriers.ksPixel - ks) < 0.75);
assert(norm(estimate.carriers.ktPixel - kt) < 0.75);
assert(angularError(estimate.phaseOffsetS, truePhaseS) < 0.12);
assert(angularError(estimate.phaseOffsetT, truePhaseT) < 0.12);
assert(abs(estimate.modulationS - trueModulationS) < 0.15);
assert(abs(estimate.modulationT - trueModulationT) < 0.15);
assert(strcmp(estimate.diagnostics.estimationMode, 'stack-phase-search'));
end

function stackFFT = makeFrequencyStack(image, ks, kt, params)
components = cat(3, ...
    FFT2D(image, false), ...
    FFT2D(latticeFourierShift(image, ks(1), ks(2)), false), ...
    FFT2D(latticeFourierShift(image, -ks(1), -ks(2)), false), ...
    FFT2D(latticeFourierShift(image, kt(1), kt(2)), false), ...
    FFT2D(latticeFourierShift(image, -kt(1), -kt(2)), false));

W = makeLatticePhaseMatrix(params);
[h, w, ~] = size(components);
componentRows = reshape(components, h * w, 5);
stackRows = componentRows * transpose(W);
stackFFT = reshape(stackRows, h, w, 5);
end

function err = angularError(actual, expected)
err = abs(angle(exp(1i * (actual - expected))));
end
