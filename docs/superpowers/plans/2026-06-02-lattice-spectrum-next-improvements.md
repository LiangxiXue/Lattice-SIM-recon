# Lattice Spectrum Next Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the visible central square and axis-aligned discontinuities in the reconstructed Lattice-SIM FFT by improving spectrum synthesis normalization, overlap blending, and diagnostics.

**Architecture:** Keep the current double-size Fourier canvas and existing reconstruction entry points. Replace hard support transitions and unstable per-band normalization with explicit per-band weighted accumulation, overlap diagnostics, smooth confidence masks, and a final display crop/normalization stage that does not hide frequency-domain defects.

**Tech Stack:** MATLAB R2024a, existing `lattice_sim_core/functions`, existing MATLAB test runner `run_all_lattice_sim_tests`.

---

## Current Verified State

- Current repo path: `/Users/xueliangxi/Lattice-SIM-recon`.
- Main reconstruction code: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions`.
- Current output generator: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/simulate_testpat_lattice_sim.m`.
- Already changed and tested:
  - frequency synthesis now uses a double-size Fourier canvas;
  - explicit `finalSupportMask` exists;
  - radial Gaussian apodization is available and enabled by default;
  - final inverse transform was changed from `abs(ifft2c(...))` to `real(ifft2c(...))`;
  - `run_all_lattice_sim_tests` passes after the `real(ifft)` change.
- Observed result:
  - `real(ifft)` alone has small visual impact;
  - the remaining obvious square is therefore more likely caused by hard transitions in support/overlap weighting, per-band OTF normalization, or inconsistent blending between the original central band and shifted sidebands.

---

## Files To Modify

- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
  - Main work: split per-band contributions, normalize weights, compute smooth overlap/confidence maps, apply final support.
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/applyLatticeApodization.m`
  - Main work: support confidence-mask input and avoid creating hard square boundaries.
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/defaultLatticeSIMParams.m`
  - Main work: add tunable blend thresholds and default values.
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/validateLatticeSIMParams.m`
  - Main work: validate new parameters.
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m`
  - Main work: add coverage for smooth masks and overlap map behavior.
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_reconstructLatticeSIM_smoke.m`
  - Main work: assert new diagnostics are present and finite.
- Optional create: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_latticeSpectrumWeights.m`
  - Use only if weight helper functions are split out of `combineLatticeSpectrum.m`.

---

## Task 1: Add Quantitative Diagnostics For The Square Artifact

**Files:**
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m`

- [ ] **Step 1: Write the failing test**

Add assertions after the existing diagnostics checks in `test_combineLatticeSpectrum_modulation.m`:

```matlab
assert(isfield(diagnostics, 'bandWeightSum'));
assert(isfield(diagnostics, 'bandCoverageCount'));
assert(isfield(diagnostics, 'overlapTransitionMap'));
assert(isequal(size(diagnostics.bandWeightSum), size(diagnostics.simSpectrum)));
assert(isequal(size(diagnostics.bandCoverageCount), size(diagnostics.simSpectrum)));
assert(isequal(size(diagnostics.overlapTransitionMap), size(diagnostics.simSpectrum)));
assert(all(isfinite(diagnostics.bandWeightSum(:))));
assert(all(isfinite(diagnostics.overlapTransitionMap(:))));
assert(max(diagnostics.bandCoverageCount(:)) >= 3);
```

- [ ] **Step 2: Run the test and confirm failure**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_combineLatticeSpectrum_modulation"
```

Expected: failure because the new diagnostics do not exist yet.

- [ ] **Step 3: Implement minimal diagnostics**

In `combineLatticeSpectrum.m`, initialize per-band tracking before the first band is added:

```matlab
bandWeightSum = zeros(outputSize);
bandCoverageCount = zeros(outputSize);
```

For each band, compute a band reliability weight from the shifted OTF before adding it:

```matlab
bandWeight = abs(shiftedOtf).^2;
bandWeight(~otfMask) = 0;
bandWeightSum = bandWeightSum + bandWeight;
bandCoverageCount = bandCoverageCount + double(otfMask);
```

For the central band, use:

```matlab
bandWeight = abs(otf0).^2;
bandWeight(~otf0Mask) = 0;
bandWeightSum = bandWeightSum + bandWeight;
bandCoverageCount = bandCoverageCount + double(otf0Mask);
```

At the end, store:

```matlab
diagnostics.bandWeightSum = bandWeightSum;
diagnostics.bandCoverageCount = bandCoverageCount;
diagnostics.overlapTransitionMap = abs(del2(double(bandCoverageCount)));
```

- [ ] **Step 4: Run targeted test**

Run the same command from Step 2.

Expected: `test_combineLatticeSpectrum_modulation passed.`

- [ ] **Step 5: Commit**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m
git commit -m "test: add lattice spectrum overlap diagnostics"
```

---

## Task 2: Replace Hard OTF Masks With Smooth Per-Band Taper Masks

**Files:**
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/defaultLatticeSIMParams.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/validateLatticeSIMParams.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m`

- [ ] **Step 1: Write failing parameter and behavior tests**

In `test_combineLatticeSpectrum_modulation.m`, add:

```matlab
assert(isfield(params, 'otfTaperLow'));
assert(isfield(params, 'otfTaperHigh'));
assert(params.otfTaperLow < params.otfTaperHigh);
assert(isfield(diagnostics, 'bandTaperMasks'));
assert(numel(diagnostics.bandTaperMasks) == 5);
for maskIdx = 1:numel(diagnostics.bandTaperMasks)
    taperMask = diagnostics.bandTaperMasks{maskIdx};
    assert(all(taperMask(:) >= 0));
    assert(all(taperMask(:) <= 1));
end
```

- [ ] **Step 2: Run and confirm failure**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_combineLatticeSpectrum_modulation"
```

Expected: failure because `otfTaperLow`, `otfTaperHigh`, and `bandTaperMasks` do not exist.

- [ ] **Step 3: Add default parameters**

In `defaultLatticeSIMParams.m`, add defaults near the existing support and OTF thresholds:

```matlab
defaults.otfTaperLow = 1e-4;
defaults.otfTaperHigh = 5e-3;
```

In `validateLatticeSIMParams.m`, add:

```matlab
mustBePositiveScalar(params.otfTaperLow, 'otfTaperLow');
mustBePositiveScalar(params.otfTaperHigh, 'otfTaperHigh');
if params.otfTaperLow >= params.otfTaperHigh
    error('LatticeSIM:InvalidParameter', 'otfTaperLow must be smaller than otfTaperHigh.');
end
```

- [ ] **Step 4: Implement smooth taper helper**

In `combineLatticeSpectrum.m`, add a local helper:

```matlab
function taper = smoothOtfTaper(otfValues, params)
maxOtf = max(otfValues(:));
low = maxOtf * params.otfTaperLow;
high = maxOtf * params.otfTaperHigh;
t = (otfValues - low) ./ (high - low + eps);
t = min(max(t, 0), 1);
taper = t .* t .* (3 - 2*t);
end
```

Use this taper instead of binary-only multiplication:

```matlab
otfTaper = smoothOtfTaper(shiftedOtf, params);
otfMask = otfTaper > 0;
freq = freq .* otfTaper;
```

For the central band:

```matlab
otf0Taper = smoothOtfTaper(otf0, params);
otf0Mask = otf0Taper > 0;
freq0 = freq0 .* otf0Taper;
```

Store:

```matlab
bandTaperMasks = {otf0Taper};
bandTaperMasks{end + 1} = otfTaper;
diagnostics.bandTaperMasks = bandTaperMasks;
```

- [ ] **Step 5: Run targeted test**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_combineLatticeSpectrum_modulation"
```

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/defaultLatticeSIMParams.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/private/validateLatticeSIMParams.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m
git commit -m "fix: smooth lattice otf band tapers"
```

---

## Task 3: Change Final Support From Binary Cutoff To Smooth Confidence Mask

**Files:**
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/applyLatticeApodization.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_reconstructLatticeSIM_smoke.m`

- [ ] **Step 1: Write failing diagnostics tests**

In `test_reconstructLatticeSIM_smoke.m`, add after existing combine diagnostics checks:

```matlab
assert(isfield(result.diagnostics.combine, 'finalConfidenceMask'));
assert(isequal(size(result.diagnostics.combine.finalConfidenceMask), imageSize * 2));
assert(all(result.diagnostics.combine.finalConfidenceMask(:) >= 0));
assert(all(result.diagnostics.combine.finalConfidenceMask(:) <= 1));
assert(any(result.diagnostics.combine.finalConfidenceMask(:) > 0 & result.diagnostics.combine.finalConfidenceMask(:) < 1));
```

- [ ] **Step 2: Run and confirm failure**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_reconstructLatticeSIM_smoke"
```

Expected: failure because `finalConfidenceMask` does not exist.

- [ ] **Step 3: Implement confidence mask**

In `combineLatticeSpectrum.m`, replace the hard-only final mask application:

```matlab
combinedSpectrum(~finalSupportMask) = 0;
```

with a smooth confidence calculation:

```matlab
denMax = max(den(:));
supportMask = den > denMax * params.supportThreshold;
reliabilityMask = den > denMax * params.reliabilityThreshold;
finalSupportMask = supportMask & reliabilityMask;

low = denMax * params.supportThreshold;
high = denMax * params.reliabilityThreshold;
confidence = (den - low) ./ (high - low + eps);
confidence = min(max(confidence, 0), 1);
finalConfidenceMask = confidence .* confidence .* (3 - 2*confidence);
combinedSpectrum = combinedSpectrum .* finalConfidenceMask;
```

Store:

```matlab
diagnostics.finalConfidenceMask = finalConfidenceMask;
```

- [ ] **Step 4: Ensure apodization respects confidence**

If `applyLatticeApodization.m` currently only receives `finalSupportMask`, keep the signature stable and pass a weighted support:

```matlab
[combinedSpectrum, apodizationMask] = applyLatticeApodization(combinedSpectrum, finalConfidenceMask > 0, params);
combinedSpectrum = combinedSpectrum .* finalConfidenceMask;
```

If this double-applies too much attenuation, move the `finalConfidenceMask` multiplication after apodization and compare FFT output.

- [ ] **Step 5: Run targeted test**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_reconstructLatticeSIM_smoke"
```

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/applyLatticeApodization.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_reconstructLatticeSIM_smoke.m
git commit -m "fix: smooth final lattice spectrum support"
```

---

## Task 4: Normalize Band Contributions By Stable Weight Sum

**Files:**
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m`

- [ ] **Step 1: Add a stability test**

In `test_combineLatticeSpectrum_modulation.m`, add:

```matlab
valid = diagnostics.finalSupportMask & diagnostics.bandWeightSum > 0;
weightedSpectrum = abs(diagnostics.simSpectrum(valid));
assert(all(isfinite(weightedSpectrum)));
assert(prctile(weightedSpectrum, 99) / max(prctile(weightedSpectrum, 50), eps) < 1e6);
```

- [ ] **Step 2: Run and capture baseline**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_combineLatticeSpectrum_modulation"
```

Expected: may pass or fail; record the percentile ratio for comparison if needed.

- [ ] **Step 3: Replace denominator with explicit stable blend denominator**

In `combineLatticeSpectrum.m`, keep `den` as the physical Wiener denominator, but add a separate `blendDen`:

```matlab
blendDen = zeros(outputSize);
```

For each band:

```matlab
blendWeight = abs(shiftedOtf).^2 .* otfTaper;
blendDen = blendDen + blendWeight;
num = num + freq .* conj(shiftedOtf) .* blendWeight ./ max(sidebandAmplitude, eps);
```

For the central band:

```matlab
blendWeight0 = abs(otf0).^2 .* otf0Taper;
blendDen = blendDen + blendWeight0;
num = freq0 .* conj(otf0) .* blendWeight0;
```

Then reconstruct:

```matlab
combinedSpectrum = num ./ (blendDen + params.wiener);
```

Keep `diagnostics.wienerDenominator = blendDen` or add a new field:

```matlab
diagnostics.physicalOtfDenominator = den;
diagnostics.blendDenominator = blendDen;
```

- [ ] **Step 4: Run targeted tests**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('functions'); addpath('tests'); addpath('tests/helpers'); test_combineLatticeSpectrum_modulation; test_reconstructLatticeSIM_smoke;"
```

Expected: both pass.

- [ ] **Step 5: Commit**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m
git commit -m "fix: stabilize lattice band contribution blending"
```

---

## Task 5: Regenerate Test Pattern And Compare FFT Regions

**Files:**
- Modify only generated output files under `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/testpat_lattice_simulation_output`

- [ ] **Step 1: Run full tests**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('tests'); run_all_lattice_sim_tests"
```

Expected: `All Lattice-SIM tests passed.`

- [ ] **Step 2: Regenerate the visual output**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); simulate_testpat_lattice_sim"
```

Expected output:

```text
Lattice-SIM result: /Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/testpat_lattice_simulation_output/Lattice-SIM.tif
```

- [ ] **Step 3: Confirm diagnostics mode and output timestamp**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); load('testpat_lattice_simulation_output/diagnostics/result.mat'); fprintf('mode=%s min=%g max=%g mean=%g\\n', result.diagnostics.combine.outputImageMode, min(result.SIM(:)), max(result.SIM(:)), mean(result.SIM(:))); info=dir('testpat_lattice_simulation_output/Lattice-SIM.tif'); fprintf('saved=%s bytes=%d\\n', info.date, info.bytes);"
```

Expected:

```text
mode=real
```

- [ ] **Step 4: Compare FFT zones numerically**

Run a MATLAB script that measures the central square, outer axis strips, and corners on `abs(fftshift(fft2(result.SIM)))`. The central square should no longer have a sharp rectangular step relative to immediately adjacent regions. If the square remains visible, inspect `diagnostics.bandCoverageCount`, `diagnostics.bandWeightSum`, and `diagnostics.finalConfidenceMask`; the defect should align with one of these maps.

- [ ] **Step 5: Commit final reconstruction changes**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests
git commit -m "fix: reduce lattice spectrum stitching discontinuities"
```

Do not commit generated TIFF/MAT output unless explicitly requested.

---

## Task 6: If The Square Remains, Investigate Phase And Sideband Consistency

**Files:**
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/combineLatticeSpectrum.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions/estimateLatticeCarrier.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_estimateLatticeCarrier.m`
- Modify: `/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests/test_combineLatticeSpectrum_modulation.m`

- [ ] **Step 1: Add diagnostics for sideband phase references**

In `combineLatticeSpectrum.m`, store:

```matlab
diagnostics.sidebandPhaseReference(idx) = angle(phaseReference);
diagnostics.sidebandPhaseMagnitude(idx) = abs(phaseReference);
```

- [ ] **Step 2: Add tests for finite phase diagnostics**

In `test_combineLatticeSpectrum_modulation.m`, assert:

```matlab
assert(isfield(diagnostics, 'sidebandPhaseReference'));
assert(numel(diagnostics.sidebandPhaseReference) == 4);
assert(all(isfinite(diagnostics.sidebandPhaseReference)));
```

- [ ] **Step 3: Replace single-center phase reference if unstable**

If `sidebandPhaseMagnitude` is close to zero, estimate phase over a small reliable overlap window instead of one center pixel:

```matlab
overlap = finalSupportMask & abs(freq0) > max(abs(freq0(:))) * 1e-3 & abs(freq) > max(abs(freq(:))) * 1e-3;
phaseReference = sum(freq(overlap) .* conj(freq0(overlap)), 'all') ./ max(nnz(overlap), 1);
```

- [ ] **Step 4: Run tests and regenerate output**

Run:

```bash
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); addpath('tests'); run_all_lattice_sim_tests"
matlab -batch "cd('Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core'); simulate_testpat_lattice_sim"
```

Expected: tests pass; FFT no longer shows phase-driven cross or square discontinuity. If not, compare positive and negative sideband amplitudes separately.

- [ ] **Step 5: Commit**

```bash
git add Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/functions Lattice_SIM_Reconstruction_Program_matlab/lattice_sim_core/tests
git commit -m "fix: stabilize lattice sideband phase alignment"
```

---

## Execution Order

1. Task 1 first: make the artifact measurable.
2. Task 2 next: remove hard per-band OTF edges.
3. Task 3 next: remove hard final support boundary.
4. Task 4 next: stabilize overlap weighting.
5. Task 5 after each meaningful change: regenerate and inspect output.
6. Task 6 only if the square persists after support/overlap fixes.

## Success Criteria

- `run_all_lattice_sim_tests` passes.
- `simulate_testpat_lattice_sim` regenerates `Lattice-SIM.tif`.
- `result.diagnostics.combine.outputImageMode` remains `real`.
- The FFT of `Lattice-SIM.tif` no longer has a sharp central square boundary aligned to the original embedded spectrum.
- Diagnostic maps identify remaining defects clearly:
  - if the defect follows `bandCoverageCount`, it is a stitching coverage problem;
  - if it follows `bandWeightSum` or `blendDenominator`, it is an overlap weighting problem;
  - if it follows neither, inspect sideband phase consistency.

