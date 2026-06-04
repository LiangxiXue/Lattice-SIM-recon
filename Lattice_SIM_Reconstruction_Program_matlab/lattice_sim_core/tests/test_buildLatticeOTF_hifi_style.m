function test_buildLatticeOTF_hifi_style()
%TEST_BUILDLATTICEOTF_HIFI_STYLE Verify HiFi-style OTF damping is in the OTF itself.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;

otf = buildLatticeOTF(64, 64, params);
center = [33, 33];
assert(abs(otf.values(center(1), center(2)) - ...
    (1 - params.otfAttenuationStrength)) < 1e-12, ...
    'Default OTF should include HiFi-style attenuation in generated values.');
assert(isfield(otf, 'hifiOtfA'), 'OTF diagnostics should record the HiFi damping factor.');
assert(isfield(otf, 'attenuationStrength'), 'OTF diagnostics should record attenuation strength.');

idealParams = params;
idealParams.otfAttenuationStrength = 0;
idealParams.hifiOtfA = 1.0;
idealOtf = buildLatticeOTF(64, 64, idealParams);
expectedIdeal = idealIncoherentOtf(idealOtf.fxCyclesPerNm, idealOtf.fyCyclesPerNm, ...
    idealOtf.cutoffCyclesPerNm);
assert(max(abs(idealOtf.values(:) - expectedIdeal(:))) < 1e-12, ...
    'Zero attenuation with hifiOtfA=1 should recover the ideal incoherent OTF.');

dampedParams = idealParams;
dampedParams.hifiOtfA = 0.85;
dampedOtf = buildLatticeOTF(64, 64, dampedParams);
midMask = expectedIdeal > 0.1 & expectedIdeal < 0.8;
assert(any(dampedOtf.values(midMask) < idealOtf.values(midMask)), ...
    'hifiOtfA below one should damp nonzero high-frequency OTF values.');
end

function values = idealIncoherentOtf(fx, fy, cutoff)
rho = hypot(fx, fy) / cutoff;
values = zeros(size(rho));
inside = rho <= 1;
values(inside) = (2 / pi) * (acos(rho(inside)) - ...
    rho(inside) .* sqrt(1 - rho(inside).^2));
end
