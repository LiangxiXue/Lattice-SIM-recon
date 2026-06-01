function test_normalizeSIMFrames()
%TEST_NORMALIZESIMFRAMES Verify per-frame mean normalization.

stack = zeros(4, 5, 5);
for idx = 1:5
    stack(:, :, idx) = idx;
end

[normalized, info] = normalizeSIMFrames(stack);
means = squeeze(mean(mean(normalized, 1), 2));
assert(max(abs(means - mean(1:5))) < 1e-12);
assert(numel(info.frameMeans) == 5);
assert(numel(info.scaleFactors) == 5);

zeroStack = zeros(4, 5, 5);
expectError(@() normalizeSIMFrames(zeroStack), 'LatticeSIM:InvalidFrameMean');

fprintf('test_normalizeSIMFrames passed.\n');
end
