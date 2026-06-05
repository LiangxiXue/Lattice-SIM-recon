function test_buildLatticeOTF_hifi_style()
%TEST_BUILDLATTICEOTF_HIFI_STYLE Verify OTF physics is separate from filters.

testDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(testDir);
addpath(fullfile(coreDir, 'functions'));

params = defaultLatticeSIMParams();
params.pixelSizeNm = 97.5;
params.emissionWavelengthNm = 561;
params.NA = 1.42;

otf = buildLatticeOTF(64, 64, params);
center = [33, 33];
expectedPhysical = idealIncoherentOtf(otf.fxCyclesPerNm, otf.fyCyclesPerNm, ...
    otf.cutoffCyclesPerNm);
assert(max(abs(otf.values(:) - expectedPhysical(:))) < 1e-12, ...
    'Physical OTF values should not include low-frequency attenuation or apodization.');
assert(abs(otf.values(center(1), center(2)) - 1) < 1e-12, ...
    'Physical OTF center should be normalized to one.');
assert(isfield(otf, 'hifiOtfA'), 'OTF diagnostics should record the HiFi damping factor.');
assert(isfield(otf, 'attenuationStrength'), 'OTF diagnostics should record attenuation strength.');
assert(isfield(otf, 'attenuationMask'), 'OTF diagnostics should expose attenuation separately.');
assert(otf.attenuationMask(center(1), center(2)) < otf.values(center(1), center(2)), ...
    'Separate attenuation mask should still record the optional low-frequency attenuation.');

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
assert(max(abs(dampedOtf.values(:) - idealOtf.values(:))) < 1e-12, ...
    'HiFi damping should not alter the physical OTF values.');
assert(any(dampedOtf.empiricalDampingMask(midMask) < 1), ...
    'HiFi damping should be available as a separate empirical filter.');
end

function values = idealIncoherentOtf(fx, fy, cutoff)
rho = hypot(fx, fy) / cutoff;
values = zeros(size(rho));
inside = rho <= 1;
values(inside) = (2 / pi) * (acos(rho(inside)) - ...
    rho(inside) .* sqrt(1 - rho(inside).^2));
end
