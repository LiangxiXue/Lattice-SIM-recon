# Lattice-SIM Reconstruction Core Design

Date: 2026-05-28

## Purpose

This design turns the current Lattice-SIM MATLAB demo into a reusable, testable
2D five-frame Lattice-SIM reconstruction core. The first version focuses on a
function-based algorithm library. A GUI can be added later after the core data
model, parameter estimation, and reconstruction pipeline are stable.

The implementation should preserve the Lattice-SIM physical model from the
paper and use HiFi-SIM only as an engineering reference for robust
preprocessing, parameter handling, carrier search, OTF/Wiener filtering,
apodization, diagnostics, and result management.

## Scope

In scope for the first version:

- 2D Lattice-SIM reconstruction from five phase-shifted raw images.
- Inputs from a single five-page TIFF stack, five separate TIFF files, or an
  in-memory `height x width x 5` array.
- Whole-image reconstruction only.
- Functional MATLAB API with structured parameters and structured output.
- Wide-field and SIM outputs.
- Diagnostics that make reconstruction failures inspectable.
- Optional result-saving wrapper.

Out of scope for the first version:

- MATLAB GUI.
- ROI selection.
- Automatic large-image tiling and stitching.
- 3D-SIM.
- Multi-color workflows.
- Time-series batch reconstruction.
- Direct modification of the HiFi-SIM GUI code.

## Core Constraint: Lattice-SIM Is Not Traditional 2D-SIM

The parameter estimation must not copy the traditional 2D-SIM model directly.

Traditional 2D-SIM commonly uses:

```text
3 directions x 3 phases = 9 frames
```

Each direction can be separated and estimated from its own phase stack.

The target Lattice-SIM model uses:

```text
2 orthogonal directions jointly illuminated x 5 phase combinations = 5 frames
```

Each raw frame contains the combined contribution:

```text
dc + s+ + s- + t+ + t-
```

Therefore the algorithm must first perform five-frame Lattice-SIM demodulation
to recover the five frequency components:

```text
C0, Cs+, Cs-, Ct+, Ct-
```

Only after this joint demodulation should it estimate the two carrier vectors:

```text
ks = [ks_x, ks_y]
kt = [kt_x, kt_y]
```

HiFi-SIM can inform the carrier-search mechanics, such as low-frequency notch
suppression and local peak refinement, but not the acquisition model or band
separation layout.

## Public API

The main reconstruction entry point should be:

```matlab
result = reconstructLatticeSIM(input, params)
```

Supported `input` forms:

- A character/string path to a five-page TIFF stack.
- A 1x5 or 5x1 cell array/string array of TIFF file paths.
- A numeric `height x width x 5` image stack.

The output should be a structure:

```matlab
result.WF
result.SIM
result.params
result.diagnostics
```

Saving should be separate from reconstruction:

```matlab
saveLatticeSIMResult(result, outputDir)
```

This separation keeps the algorithm usable from scripts, tests, future batch
tools, and a future GUI.

## Module Breakdown

### `readLatticeSIMInput`

Responsibilities:

- Accept the public `input` forms.
- Return a numeric `height x width x 5` stack.
- Preserve original image precision information in metadata where useful.
- Validate that the input contains exactly five frames.
- Validate that all frames have identical dimensions.

### `defaultLatticeSIMParams`

Responsibilities:

- Provide conservative default parameters.
- Keep defaults close to the current demo and paper where appropriate.
- Allow user overrides through `params`.

Initial default candidates:

```matlab
params.pixelSizeNm = [];
params.emissionWavelengthNm = [];
params.NA = [];
params.modulationS = 0.7;
params.modulationT = 0.7;
params.wiener = 0.04;
params.apodizationEnabled = false;
params.apodizationStrength = 0.4;
params.notchScale = 0.5;
params.outputScaleMode = "none";
```

Required physical parameters for OTF-based reconstruction:

```matlab
pixelSizeNm
emissionWavelengthNm
NA
```

### `normalizeSIMFrames`

Responsibilities:

- Normalize each raw frame to a common mean intensity.
- Return both normalized stack and normalization factors.
- Avoid division by zero for empty or saturated frames.

### `separateLatticeBands`

Responsibilities:

- Implement the five-frame Lattice-SIM phase-demodulation matrix.
- Return the five complex components:

```matlab
bands.C0
bands.CsPlus
bands.CsMinus
bands.CtPlus
bands.CtMinus
```

The phase matrix should match the Lattice-SIM five phase combinations:

```text
[0, 0]
[0, 2*pi/3]
[0, 4*pi/3]
[2*pi/3, 0]
[4*pi/3, 2*pi/3]
```

The implementation should avoid using the traditional `nrDirs x nrPhases`
separation pattern from HiFi-SIM.

### `estimateLatticeCarrier`

Responsibilities:

- Estimate `ks` from `CsPlus/CsMinus`.
- Estimate `kt` from `CtPlus/CtMinus`.
- Use a notch or low-frequency suppression before peak selection.
- Optionally refine the peak locally using correlation-style refinement inspired
  by HiFi-SIM.
- Return carrier vectors in both pixel-frequency coordinates and radians per
  pixel where useful.

Diagnostics should include:

```matlab
diagnostics.carrierS
diagnostics.carrierT
diagnostics.carrierMagnitudeS
diagnostics.carrierMagnitudeT
diagnostics.carrierAngleDeg
diagnostics.orthogonalityErrorDeg
diagnostics.carrierSearchMaps
```

Warnings should be generated when:

- A carrier peak is too close to the DC region.
- The two carrier directions are not close to orthogonal.
- Carrier magnitudes differ unexpectedly.
- The selected peak is weak or ambiguous.

### `buildLatticeOTF`

Responsibilities:

- Build a simple theoretical OTF from `NA`, emission wavelength, and pixel size.
- Expose the frequency grid needed by spectrum combination.
- Keep the OTF model replaceable so that measured OTF support can be added
  later.

The first version can use a simplified OTF close to the current demo or a
HiFi-SIM-style theoretical OTF helper, but it should be isolated from the rest
of the pipeline.

### `combineLatticeSpectrum`

Responsibilities:

- Fourier-transform the separated components.
- Shift `Cs+`, `Cs-`, `Ct+`, and `Ct-` to their correct Fourier-space locations.
- Align phases relative to the DC term.
- Weight components with the OTF.
- Combine the spectrum using Wiener-style damping.
- Produce the wide-field spectrum and reconstructed SIM spectrum.

The Lattice-SIM component layout should remain:

```text
C0, Cs+, Cs-, Ct+, Ct-
```

No assumption should be introduced that there are three illumination directions.

### `applyLatticeApodization`

Responsibilities:

- Optionally apply apodization to the combined spectrum.
- Reduce ringing artifacts from abrupt spectral support edges.
- Be disabled by default for easier baseline comparison with the original demo.

### `saveLatticeSIMResult`

Responsibilities:

- Save `result.WF` and `result.SIM`.
- Optionally save selected diagnostic images.
- Avoid changing the algorithm result structure.

Suggested output names:

```text
Wide-field.tif
Lattice-SIM.tif
diagnostics/
```

## Data Flow

```text
input
  -> readLatticeSIMInput
  -> rawStack(height x width x 5)
  -> normalizeSIMFrames
  -> separateLatticeBands
  -> estimateLatticeCarrier
  -> buildLatticeOTF
  -> combineLatticeSpectrum
  -> applyLatticeApodization, optional
  -> inverse FFT with 2x zero filling
  -> result.WF, result.SIM, result.params, result.diagnostics
```

## Error Handling

The core should fail early with clear messages for:

- Input not containing exactly five frames.
- Inconsistent frame sizes.
- Non-numeric or empty image stack.
- Missing required physical parameters when OTF reconstruction is requested.
- Invalid modulation factors.
- Invalid Wiener parameter.
- Failure to estimate either carrier.
- `NaN` or `Inf` in intermediate spectra or final images.

Non-fatal quality issues should be placed in:

```matlab
result.diagnostics.warnings
```

Examples:

- Carrier peak is weak.
- Carrier directions are not close to orthogonal.
- Estimated modulation appears very low.
- Apodization altered a large fraction of the spectrum.

## Testing Strategy

The first test layer should verify code behavior rather than biological image
quality:

- Input loading test for a five-page TIFF stack.
- Input loading test for five separate TIFF files.
- Input loading test for an in-memory `height x width x 5` array.
- Error test for four frames and six frames.
- Error test for inconsistent frame sizes.
- Synthetic five-frame Lattice-SIM demodulation test using known phase
  combinations.
- Carrier-estimation test on synthetic sinusoidal components with known
  `ks/kt`.
- Baseline comparison against the current demo on the same data, when demo
  data are available.

The second test layer should verify reconstruction diagnostics:

- Carrier vectors are recorded.
- Orthogonality error is recorded.
- Warnings are emitted for invalid or ambiguous carrier maps.
- Output dimensions are doubled when 2x spectral zero filling is enabled.

## Migration Plan

1. Keep the original demo scripts unchanged as reference files.
2. Add a new function-based implementation under a clear folder such as
   `Lattice_SIM_Reconstruction_Program_matlab/functions`.
3. Build the reader and five-frame demodulation first.
4. Add Lattice-specific carrier estimation.
5. Add OTF/Wiener/apodization modules.
6. Add a small example script that calls `reconstructLatticeSIM`.
7. Add tests and compare against the original demo behavior.

## Open Decisions

These decisions should be made during implementation planning:

- Exact folder layout for new functions.
- Whether tests should use MATLAB's `matlab.unittest` or simpler assertion
  scripts.
- Whether to include a synthetic data generator in the public API or keep it
  test-only.
- Whether first-version OTF should match the current demo exactly or use the
  cleaner theoretical form from HiFi-SIM.

