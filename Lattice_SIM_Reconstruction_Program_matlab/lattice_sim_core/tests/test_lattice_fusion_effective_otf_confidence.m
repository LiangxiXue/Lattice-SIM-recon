function test_lattice_fusion_effective_otf_confidence()
%TEST_LATTICE_FUSION_EFFECTIVE_OTF_CONFIDENCE Guard final fusion damping.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;
params.hifiOtfA = 0.55;
params.otfAttenuationEnabled = true;
params.otfAttenuationStrength = 0.35;
params.otfAttenuationFwhm = 0.5;
params.supportThreshold = 1e-4;
params.reliabilityThreshold = 0.2;
params.apodizationEnabled = false;

bands = makeFrequencyBands([48, 48]);
carriers.ksPixel = [7, 0];
carriers.ktPixel = [0, 8];
otf = buildLatticeOTF(48, 48, params);

[~, diagnostics] = combineLatticeSpectrum(bands, carriers, otf, params);

assert(isfield(diagnostics, 'effectiveOtfValues'), ...
    'Fusion diagnostics should expose the damped effective OTF actually used for weighting.');
assert(max(abs(diagnostics.effectiveOtfValues(:) - diagnostics.physicalOtfValues(:))) > 1e-3, ...
    'Effective OTF should differ from the physical OTF when damping/attenuation is enabled.');
assert(all(abs(diagnostics.effectiveOtfValues(:)) <= ...
    abs(diagnostics.physicalOtfValues(:)) + 1e-12), ...
    'Effective OTF should not amplify the physical OTF support.');

transitionMask = diagnostics.supportMask & ~diagnostics.reliabilityMask;
assert(nnz(transitionMask) > 0, ...
    'Test setup should create a support-to-reliability transition region.');
transitionConfidence = diagnostics.finalConfidenceMask(transitionMask);
assert(any(transitionConfidence(:) > 0 & transitionConfidence(:) < 1), ...
    'Final confidence should taper smoothly between support and reliability thresholds.');
assert(all(diagnostics.finalConfidenceMask(~diagnostics.supportMask) == 0), ...
    'Unsupported frequencies should still be removed.');
assert(all(diagnostics.finalConfidenceMask(diagnostics.reliabilityMask) == 1), ...
    'Reliable frequencies should retain full confidence.');
end

function bands = makeFrequencyBands(imageSize)
[x, y] = meshgrid(linspace(-1, 1, imageSize(2)), linspace(-1, 1, imageSize(1)));
base = exp(-12 * (x .^ 2 + y .^ 2));
spectrum = FFT2D(base, false);

bands.C0 = spectrum;
bands.CsPlus = 0.35 * spectrum;
bands.CsMinus = 0.33 * spectrum;
bands.CtPlus = 0.31 * spectrum;
bands.CtMinus = 0.29 * spectrum;
bands.domain = 'frequency';
end
