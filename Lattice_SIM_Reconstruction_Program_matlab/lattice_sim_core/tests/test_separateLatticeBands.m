function test_separateLatticeBands()
%TEST_SEPARATELATTICEBANDS Verify exact five-frame Lattice demodulation.

[rawStack, expectedBands] = makeSyntheticLatticeSIMStack([16, 18], ...
    'ksPixel', [3, 0], 'ktPixel', [0, 4], 'carrierAmplitude', 0);

bands = separateLatticeBands(rawStack);

assert(max(abs(bands.C0(:) - expectedBands.C0(:))) < 1e-10);
assert(max(abs(bands.CsPlus(:) - expectedBands.CsPlus(:))) < 1e-10);
assert(max(abs(bands.CsMinus(:) - expectedBands.CsMinus(:))) < 1e-10);
assert(max(abs(bands.CtPlus(:) - expectedBands.CtPlus(:))) < 1e-10);
assert(max(abs(bands.CtMinus(:) - expectedBands.CtMinus(:))) < 1e-10);
assert(isequal(size(bands.phaseMatrix), [5, 5]));

fprintf('test_separateLatticeBands passed.\n');
end
