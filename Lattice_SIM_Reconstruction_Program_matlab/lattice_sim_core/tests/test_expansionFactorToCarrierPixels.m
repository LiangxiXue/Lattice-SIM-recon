function test_expansionFactorToCarrierPixels()
%TEST_EXPANSIONFACTORTOCARRIERPIXELS Verify SIM expansion-factor carrier setup.

imageSize = [256, 256];
pixelSizeNm = 65;
emissionWavelengthNm = 532;
NA = 1.2;

[ksPixel, ktPixel, carrierMagnitudePixels] = expansionFactorToCarrierPixels( ...
    1.5, imageSize, pixelSizeNm, emissionWavelengthNm, NA);

expectedMagnitude = (1.5 - 1) * (2 * NA / emissionWavelengthNm) ...
    * imageSize(2) * pixelSizeNm;

assert(abs(carrierMagnitudePixels - expectedMagnitude) < 1e-12);
assert(isequal(ksPixel, [38, 0]));
assert(isequal(ktPixel, [0, 38]));

try
    expansionFactorToCarrierPixels(2.1, imageSize, pixelSizeNm, emissionWavelengthNm, NA);
    error('Expected expansion-factor validation to fail.');
catch err
    assert(strcmp(err.identifier, 'LatticeSIM:InvalidExpansionFactor'));
end

fprintf('test_expansionFactorToCarrierPixels passed.\n');
end
