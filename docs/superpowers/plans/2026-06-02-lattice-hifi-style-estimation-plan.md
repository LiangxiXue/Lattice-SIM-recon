# Lattice HiFi-Style Estimation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Lattice-SIM reconstruction pipeline that keeps HiFi-SIM-style preprocessing and two-pass parameter estimation, while using a Lattice-specific frequency-domain separation model.

**Architecture:** Keep `deconvlucy + FFT2D` preprocessing from the current `reconstructLatticeSIM` path. Replace direct HiFi-SIM `separateBands` use with a frequency-domain Lattice phase matrix, then estimate `ks/kt`, phase offsets, and modulation from C0/sideband overlap before a second frequency-domain separation. Leave `applyOtf`, HiFi-SIM W1/W2 Wiener filters, and `apodize_gauss` out of scope for this plan.

**Tech Stack:** MATLAB R2024a, `lattice_sim_core/functions`, local MATLAB function tests launched with `/Applications/MATLAB_R2024a.app/bin/matlab -batch`.

---

## Scope Boundary

This plan implements:

- HiFi-SIM-style preprocessing already requested: `importImages`, `deconvlucy(..., psf, 5)`, `FFT2D`.
- Lattice-specific frequency-domain separation, not HiFi-SIM `comp = [0, 1, -1, 2, -2]`.
- First-pass ideal separation using `latticePhasePairs`.
- Estimation of `ks`, `kt`, `phiS_offset`, `phiT_offset`, `mS`, and `mT`.
- Second-pass separation using estimated phase offsets and modulation.
- Existing `combineLatticeSpectrum` remains the final synthesis backend.

This plan does not implement:

- HiFi-SIM `applyOtf`.
- HiFi-SIM `WienerFilterW1_*` / `WienerFilterW2_*`.
- HiFi-SIM `writeApoVector`.
- HiFi-SIM `apodize_gauss`.

## File Structure

- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/separateLatticeBandsFrequency.m`
  - Responsibility: frequency-domain Lattice-SIM demodulation using `[C0, Cs+, Cs-, Ct+, Ct-]`.
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/estimateLatticeBandParameters.m`
  - Responsibility: estimate carrier vectors, phase offsets, and modulation values from first-pass frequency-domain bands.
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/latticeCommonRegion.m`
  - Responsibility: select overlap regions between C0 and shifted sidebands for robust correlation, adapted from HiFi-SIM `commonRegion` concept.
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/latticeBandCorrelation.m`
  - Responsibility: compute a complex C0/sideband correlation. `angle` gives phase offset; `abs` gives modulation.
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/reconstructLatticeSIM.m`
  - Responsibility: call first-pass separation, parameter estimation, second-pass separation, then combine.
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/defaultLatticeSIMParams.m`
  - Responsibility: add estimation defaults.
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/validateLatticeSIMParams.m`
  - Responsibility: validate estimation defaults.
- Modify or replace: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_hifi_style_frequency_pipeline.m`
  - Responsibility: assert the final path uses Lattice-specific frequency separation, not HiFi-SIM `separateBands`.
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_separateLatticeBandsFrequency.m`
  - Responsibility: validate the frequency-domain Lattice matrix on synthetic known components.
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_estimateLatticeBandParameters.m`
  - Responsibility: validate phase/modulation estimation on known synthetic frequency bands.

---

## Task 1: Add Lattice-Specific Frequency-Domain Separation

**Files:**
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/separateLatticeBandsFrequency.m`
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_separateLatticeBandsFrequency.m`

- [ ] **Step 1: Write the failing test**

Create `test_separateLatticeBandsFrequency.m`:

```matlab
function test_separateLatticeBandsFrequency()
testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

imageSize = [16, 16];
[x, y] = meshgrid(0:imageSize(2)-1, 0:imageSize(1)-1);
C0 = 10 + zeros(imageSize);
CsPlus = 2 * exp(1i * 2*pi*(3*x/imageSize(2)));
CsMinus = conj(CsPlus);
CtPlus = 1.5 * exp(1i * 2*pi*(4*y/imageSize(1)));
CtMinus = conj(CtPlus);

expected.C0 = fftshift(fft2(C0));
expected.CsPlus = fftshift(fft2(CsPlus));
expected.CsMinus = fftshift(fft2(CsMinus));
expected.CtPlus = fftshift(fft2(CtPlus));
expected.CtMinus = fftshift(fft2(CtMinus));

phasePairs = latticePhasePairs();
stackFFT = zeros([imageSize, 5]);
for idx = 1:5
    phiS = phasePairs(idx, 1);
    phiT = phasePairs(idx, 2);
    stackFFT(:, :, idx) = expected.C0 ...
        + expected.CsPlus .* exp(1i * phiS) ...
        + expected.CsMinus .* exp(-1i * phiS) ...
        + expected.CtPlus .* exp(1i * phiT) ...
        + expected.CtMinus .* exp(-1i * phiT);
end

params = defaultLatticeSIMParams();
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;

bands = separateLatticeBandsFrequency(stackFFT, params);

assert(max(abs(bands.C0(:) - expected.C0(:))) < 1e-8);
assert(max(abs(bands.CsPlus(:) - expected.CsPlus(:))) < 1e-8);
assert(max(abs(bands.CsMinus(:) - expected.CsMinus(:))) < 1e-8);
assert(max(abs(bands.CtPlus(:) - expected.CtPlus(:))) < 1e-8);
assert(max(abs(bands.CtMinus(:) - expected.CtMinus(:))) < 1e-8);
assert(strcmp(bands.domain, 'frequency'));
assert(strcmp(bands.diagnostics.model, 'lattice-phase-matrix'));
end
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_separateLatticeBandsFrequency"
```

Expected: failure because `separateLatticeBandsFrequency` does not exist.

- [ ] **Step 3: Implement the frequency-domain Lattice separator**

Create `separateLatticeBandsFrequency.m`:

```matlab
function bands = separateLatticeBandsFrequency(stackFFT, params)
%SEPARATELATTICEBANDSFREQUENCY Demodulate five frequency-domain Lattice-SIM frames.

if nargin < 2
    params = defaultLatticeSIMParams();
end

validateLatticeSIMStack(stackFFT);

W = latticeFrequencyPhaseMatrix(params);
[h, w, ~] = size(stackFFT);
reshaped = reshape(double(stackFFT), h * w, 5);
components = reshaped / transpose(W);
components = reshape(components, h, w, 5);

bands.C0 = components(:, :, 1);
bands.CsPlus = components(:, :, 2);
bands.CsMinus = components(:, :, 3);
bands.CtPlus = components(:, :, 4);
bands.CtMinus = components(:, :, 5);
bands.domain = 'frequency';
bands.phaseMatrix = W;
bands.diagnostics.model = 'lattice-phase-matrix';
bands.diagnostics.inputDomain = 'frequency';
bands.diagnostics.outputDomain = 'frequency';
bands.diagnostics.phaseOffsetS = params.phaseOffsetS;
bands.diagnostics.phaseOffsetT = params.phaseOffsetT;
bands.diagnostics.modulationS = params.estimatedModulationS;
bands.diagnostics.modulationT = params.estimatedModulationT;
end

function W = latticeFrequencyPhaseMatrix(params)
phasePairs = latticePhasePairs();
W = zeros(5, 5);
for idx = 1:5
    phiS = phasePairs(idx, 1) + params.phaseOffsetS;
    phiT = phasePairs(idx, 2) + params.phaseOffsetT;
    W(idx, :) = [1, ...
        params.estimatedModulationS .* exp(1i * phiS), ...
        params.estimatedModulationS .* exp(-1i * phiS), ...
        params.estimatedModulationT .* exp(1i * phiT), ...
        params.estimatedModulationT .* exp(-1i * phiT)];
end
end
```

- [ ] **Step 4: Run the test and verify it passes**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_separateLatticeBandsFrequency"
```

Expected: exit code 0.

---

## Task 2: Add Parameters For Two-Pass Lattice Estimation

**Files:**
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/defaultLatticeSIMParams.m`
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/validateLatticeSIMParams.m`
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_lattice_estimation_params.m`

- [ ] **Step 1: Write the failing parameter test**

Create `test_lattice_estimation_params.m`:

```matlab
function test_lattice_estimation_params()
testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();

assert(strcmp(char(params.separationFunction), 'separateLatticeBandsFrequency'));
assert(params.phaseOffsetS == 0);
assert(params.phaseOffsetT == 0);
assert(params.estimatedModulationS == 1);
assert(params.estimatedModulationT == 1);
assert(params.enableLatticeParameterEstimation == true);
assert(params.latticeCarrierRefinementIterations == 3);
assert(params.latticeCarrierSearchStepPixels == 2.5);
assert(params.latticeCorrelationOverlap == 0.15);
end
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_lattice_estimation_params"
```

Expected: failure because defaults still name `separateBands` and do not include the Lattice estimation fields.

- [ ] **Step 3: Add defaults**

In `defaultLatticeSIMParams.m`, replace the current direct HiFi-SIM separator defaults with:

```matlab
params.separationFunction = "separateLatticeBandsFrequency";
params.separationInputDomain = "frequency";
params.phaseOffsetS = 0;
params.phaseOffsetT = 0;
params.estimatedModulationS = 1;
params.estimatedModulationT = 1;
params.enableLatticeParameterEstimation = true;
params.latticeCarrierRefinementIterations = 3;
params.latticeCarrierSearchStepPixels = 2.5;
params.latticeCorrelationOverlap = 0.15;
```

- [ ] **Step 4: Add validation**

In `validateLatticeSIMParams.m`, replace the `separateBands`-only validation with:

```matlab
if ~any(strcmp(char(params.separationFunction), {'separateLatticeBandsFrequency'}))
    error('LatticeSIM:InvalidSeparationFunction', ...
        'Separation function must be "separateLatticeBandsFrequency".');
end
if ~isnumeric(params.phaseOffsetS) || ~isscalar(params.phaseOffsetS)
    error('LatticeSIM:InvalidPhaseOffset', 'phaseOffsetS must be a scalar.');
end
if ~isnumeric(params.phaseOffsetT) || ~isscalar(params.phaseOffsetT)
    error('LatticeSIM:InvalidPhaseOffset', 'phaseOffsetT must be a scalar.');
end
if params.estimatedModulationS <= 0 || params.estimatedModulationT <= 0
    error('LatticeSIM:InvalidEstimatedModulation', ...
        'Estimated modulation values must be positive.');
end
if params.latticeCarrierRefinementIterations < 0 || ...
        fix(params.latticeCarrierRefinementIterations) ~= params.latticeCarrierRefinementIterations
    error('LatticeSIM:InvalidCarrierRefinementIterations', ...
        'latticeCarrierRefinementIterations must be a non-negative integer.');
end
if params.latticeCarrierSearchStepPixels <= 0
    error('LatticeSIM:InvalidCarrierSearchStep', ...
        'latticeCarrierSearchStepPixels must be positive.');
end
if params.latticeCorrelationOverlap <= 0 || params.latticeCorrelationOverlap >= 1
    error('LatticeSIM:InvalidCorrelationOverlap', ...
        'latticeCorrelationOverlap must be in the open interval (0, 1).');
end
```

- [ ] **Step 5: Run the test and verify it passes**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_lattice_estimation_params"
```

Expected: exit code 0.

---

## Task 3: Implement C0/Sideband Common-Region Correlation

**Files:**
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/latticeCommonRegion.m`
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/latticeBandCorrelation.m`
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_estimateLatticeBandParameters.m`

- [ ] **Step 1: Write the failing correlation test**

Create `test_estimateLatticeBandParameters.m` with the first test:

```matlab
function test_estimateLatticeBandParameters()
test_latticeBandCorrelation_recovers_complex_scale();
end

function test_latticeBandCorrelation_recovers_complex_scale()
testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

image = makeLatticeSIMPhantom([32, 32]);
band0 = fftshift(fft2(image));
trueScale = 0.42 * exp(1i * 0.7);
band1 = trueScale .* band0;

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;
otf = buildLatticeOTF(32, 32, params);

scale = latticeBandCorrelation(band0, band1, otf.values, [0, 0], params);

assert(abs(abs(scale) - 0.42) < 1e-3);
assert(abs(angle(scale) - 0.7) < 1e-3);
end
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_estimateLatticeBandParameters"
```

Expected: failure because `latticeBandCorrelation` does not exist.

- [ ] **Step 3: Implement common-region selection**

Create `private/latticeCommonRegion.m`:

```matlab
function mask = latticeCommonRegion(otfValues, carrierPixel, params)
%LATTICECOMMONREGION Select reliable overlap between C0 and shifted sideband OTF.

[h, w] = size(otfValues);
[colGrid, rowGrid] = meshgrid(1:w, 1:h);
shiftedOtf = interp2(colGrid, rowGrid, otfValues, ...
    colGrid + carrierPixel(1), rowGrid + carrierPixel(2), 'linear', 0);

threshold = max(abs(otfValues(:))) * params.otfClipThreshold;
mask = abs(otfValues) > threshold & abs(shiftedOtf) > threshold;

[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
radius = hypot(x, y);
mask(radius < max(2, min(h, w) * params.latticeCorrelationOverlap * 0.1)) = false;
end
```

- [ ] **Step 4: Implement complex correlation**

Create `private/latticeBandCorrelation.m`:

```matlab
function scale = latticeBandCorrelation(band0, band1, otfValues, carrierPixel, params)
%LATTICEBANDCORRELATION Estimate complex sideband scale relative to C0.

mask = latticeCommonRegion(otfValues, carrierPixel, params);
if nnz(mask) < 8
    error('LatticeSIM:InsufficientOverlap', ...
        'Not enough common OTF support to estimate sideband correlation.');
end

reference = band0(mask);
target = band1(mask);
denominator = sum(abs(reference).^2, 'all');
if denominator <= eps
    error('LatticeSIM:DegenerateCorrelationReference', ...
        'C0 overlap region has zero energy.');
end

scale = sum(target .* conj(reference), 'all') ./ denominator;
end
```

- [ ] **Step 5: Run the test and verify it passes**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_estimateLatticeBandParameters"
```

Expected: exit code 0.

---

## Task 4: Estimate Lattice Carrier, Phase Offset, And Modulation

**Files:**
- Create: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/estimateLatticeBandParameters.m`
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_estimateLatticeBandParameters.m`

- [ ] **Step 1: Extend the failing test**

Append to `test_estimateLatticeBandParameters.m`:

```matlab
function test_estimateLatticeBandParameters_recovers_offsets_and_modulation()
testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;

image = makeLatticeSIMPhantom([64, 64]);
[x, y] = meshgrid(0:63, 0:63);
ks = [8, 0];
kt = [0, 10];
phiSOffset = 0.4;
phiTOffset = -0.6;
mS = 0.55;
mT = 0.35;

bands.C0 = fftshift(fft2(image));
bands.CsPlus = mS * exp(1i * phiSOffset) .* fftshift(fft2(image .* exp(1i * 2*pi*(ks(1)*x/64 + ks(2)*y/64))));
bands.CsMinus = conj(flipud(fliplr(bands.CsPlus)));
bands.CtPlus = mT * exp(1i * phiTOffset) .* fftshift(fft2(image .* exp(1i * 2*pi*(kt(1)*x/64 + kt(2)*y/64))));
bands.CtMinus = conj(flipud(fliplr(bands.CtPlus)));
bands.domain = 'frequency';

estimate = estimateLatticeBandParameters(bands, params);

assert(norm(estimate.carriers.ksPixel - ks) < 2);
assert(norm(estimate.carriers.ktPixel - kt) < 2);
assert(abs(estimate.phaseOffsetS + phiSOffset) < 0.25);
assert(abs(estimate.phaseOffsetT + phiTOffset) < 0.25);
assert(abs(estimate.modulationS - mS) < 0.2);
assert(abs(estimate.modulationT - mT) < 0.2);
end
```

Also update the top function:

```matlab
function test_estimateLatticeBandParameters()
test_latticeBandCorrelation_recovers_complex_scale();
test_estimateLatticeBandParameters_recovers_offsets_and_modulation();
end
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_estimateLatticeBandParameters"
```

Expected: failure because `estimateLatticeBandParameters` does not exist.

- [ ] **Step 3: Implement estimator**

Create `estimateLatticeBandParameters.m`:

```matlab
function estimate = estimateLatticeBandParameters(bands, params)
%ESTIMATELATTICEBANDPARAMETERS Estimate Lattice carrier, phase offsets, and modulation.

[carriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(bands.C0, 1), size(bands.C0, 2), params);

pS = latticeBandCorrelation(bands.C0, bands.CsPlus, otf.values, carriers.ksPixel, params);
pT = latticeBandCorrelation(bands.C0, bands.CtPlus, otf.values, carriers.ktPixel, params);

estimate.carriers = carriers;
estimate.phaseOffsetS = -angle(pS);
estimate.phaseOffsetT = -angle(pT);
estimate.modulationS = min(max(abs(pS), eps), 1);
estimate.modulationT = min(max(abs(pT), eps), 1);
estimate.correlationS = pS;
estimate.correlationT = pT;
estimate.diagnostics = carrierDiagnostics;
estimate.diagnostics.phaseOffsetS = estimate.phaseOffsetS;
estimate.diagnostics.phaseOffsetT = estimate.phaseOffsetT;
estimate.diagnostics.modulationS = estimate.modulationS;
estimate.diagnostics.modulationT = estimate.modulationT;
end
```

- [ ] **Step 4: Run and verify it passes**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_estimateLatticeBandParameters"
```

Expected: exit code 0.

---

## Task 5: Wire Two-Pass Separation Into reconstructLatticeSIM

**Files:**
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/reconstructLatticeSIM.m`
- Modify: `Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_hifi_style_frequency_pipeline.m`

- [ ] **Step 1: Replace the pipeline test expectation**

Update `test_hifi_style_frequency_pipeline.m` assertions:

```matlab
assert(strcmp(result.diagnostics.preprocessing.method, 'hifi-rl-fft'), ...
    'Expected reconstructLatticeSIM to use HiFi-SIM RL + FFT preprocessing.');
assert(result.diagnostics.preprocessing.deconvolutionIterations == 5, ...
    'Expected five Richardson-Lucy iterations.');
assert(strcmp(result.diagnostics.separation.functionName, 'separateLatticeBandsFrequency'), ...
    'Expected Lattice-specific frequency-domain separation.');
assert(strcmp(result.diagnostics.separation.inputDomain, 'frequency'), ...
    'Expected separation input to be frequency-domain frames.');
assert(strcmp(result.diagnostics.separation.model, 'lattice-phase-matrix'), ...
    'Expected Lattice phase matrix, not HiFi-SIM harmonic comp matrix.');
assert(isfield(result.diagnostics, 'latticeEstimation'), ...
    'Expected first-pass Lattice parameter estimation diagnostics.');
assert(isfield(result.diagnostics.latticeEstimation, 'phaseOffsetS'));
assert(isfield(result.diagnostics.latticeEstimation, 'phaseOffsetT'));
assert(isfield(result.diagnostics.latticeEstimation, 'modulationS'));
assert(isfield(result.diagnostics.latticeEstimation, 'modulationT'));
assert(strcmp(result.diagnostics.combine.bandDomain, 'frequency'), ...
    'Expected combineLatticeSpectrum to consume frequency-domain bands directly.');
assert(all(isfinite(result.SIM(:))), 'SIM output must be finite.');
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_hifi_style_frequency_pipeline"
```

Expected: failure because `reconstructLatticeSIM` still calls direct HiFi-SIM `separateBands`.

- [ ] **Step 3: Implement two-pass separation**

In `reconstructLatticeSIM.m`, replace:

```matlab
separate = separateBands(stack, params.phaOff, params.nrBands, params.bandFactors);
bands = packSeparatedFrequencyBands(separate, params);
[carriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
```

with:

```matlab
firstPassBands = separateLatticeBandsFrequency(stack, params);

if params.enableLatticeParameterEstimation
    latticeEstimate = estimateLatticeBandParameters(firstPassBands, params);
    params.phaseOffsetS = latticeEstimate.phaseOffsetS;
    params.phaseOffsetT = latticeEstimate.phaseOffsetT;
    params.estimatedModulationS = latticeEstimate.modulationS;
    params.estimatedModulationT = latticeEstimate.modulationT;
    bands = separateLatticeBandsFrequency(stack, params);
    carriers = latticeEstimate.carriers;
    carrierDiagnostics = latticeEstimate.diagnostics;
else
    bands = firstPassBands;
    [carriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
    latticeEstimate = struct();
end
```

Remove the now-obsolete local `packSeparatedFrequencyBands` helper from `reconstructLatticeSIM.m`.

Add diagnostics:

```matlab
result.diagnostics.latticeEstimation = latticeEstimate;
```

- [ ] **Step 4: Run the pipeline test**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_hifi_style_frequency_pipeline"
```

Expected: exit code 0.

---

## Task 6: Verify End-To-End Simulation Still Runs

**Files:**
- No source file changes expected.

- [ ] **Step 1: Run all new targeted tests**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests'); test_lattice_estimation_params; test_separateLatticeBandsFrequency; test_estimateLatticeBandParameters; test_hifi_style_frequency_pipeline"
```

Expected: exit code 0.

- [ ] **Step 2: Run the existing testpat simulation**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "run('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/simulate_testpat_lattice_sim.m')"
```

Expected: exit code 0 and these files are reported:

```text
cropped_testpat.tif
simulated_raw_stack.tif
Wide-field-C0.tif
Wide-field-D3D4D5.tif
Lattice-SIM.tif
```

- [ ] **Step 3: Inspect diagnostics**

Run:

```bash
/Applications/MATLAB_R2024a.app/bin/matlab -batch "load('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/testpat_lattice_simulation_output/diagnostics/result.mat'); disp(result.diagnostics.preprocessing.method); disp(result.diagnostics.separation.model); disp(result.diagnostics.latticeEstimation.phaseOffsetS); disp(result.diagnostics.latticeEstimation.phaseOffsetT); disp(result.diagnostics.latticeEstimation.modulationS); disp(result.diagnostics.latticeEstimation.modulationT);"
```

Expected:

```text
hifi-rl-fft
lattice-phase-matrix
```

The four numeric estimation values should print as finite scalars.

---

## Self-Review

- Spec coverage: The plan covers HiFi-style preprocessing retention, Lattice-specific frequency-domain separation, first-pass estimation, phase/modulation estimation, second-pass separation, and final end-to-end verification.
- Explicitly excluded: HiFi-SIM `applyOtf`, W1/W2 Wiener filters, `writeApoVector`, and `apodize_gauss`.
- Placeholder scan: No `TBD`, `TODO`, or unspecified tests remain.
- Type consistency: The plan consistently uses `phaseOffsetS`, `phaseOffsetT`, `estimatedModulationS`, `estimatedModulationT`, `separateLatticeBandsFrequency`, and `estimateLatticeBandParameters`.
