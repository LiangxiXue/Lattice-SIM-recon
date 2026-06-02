# Spectrum Shift Before Wiener Test

This folder contains a standalone HiFi-SIM debug script that stops after
frequency shifting and direct spectrum combination, before Wiener filtering.

## Run

In MATLAB:

```matlab
cd('/Users/xueliangxi/Lattice-SIM-recon/Lattice_SIM_Reconstruction_Program_matlab/HiFi-SIM_v1.01/SpectrumShiftBeforeWienerTest')
debugResult = show_spectrum_shift_before_wiener;
```

By default, the script uses:

```text
../TestData/2D-SIM(3 angles_3 phases)_9 frames_Group1.tif
```

To use another 3-angle, 3-phase raw stack:

```matlab
debugResult = show_spectrum_shift_before_wiener('/absolute/path/to/raw_stack.tif');
```

## Output

The script displays two MATLAB figures:

- `After spectrum shift/combination, before Wiener`
- `log(1 + abs(fftDirectlyCombined))`

It also writes these files under `output/`:

- `before_wiener_spatial.png`
- `before_wiener_spectrum.png`
- `before_wiener_debug.mat`

The `.mat` file includes `fftDirectlyCombined`, `beforeWienerImage`,
`frequencyMagnitude`, `IIrawFFT`, `shiftedBands`, `param`, and `K0`.

## Debug Boundary

This script reproduces the same HiFi-SIM boundary as the commented block in
`Main_fun/HiFiSIM.m`:

```matlab
Temp1 = real(ifft2(fftshift((fftDirectlyCombined))));
```

It intentionally stops before the `Wk1`, `Wk2`, and final Wiener/apodization
steps so the direct spectrum placement can be inspected independently.
