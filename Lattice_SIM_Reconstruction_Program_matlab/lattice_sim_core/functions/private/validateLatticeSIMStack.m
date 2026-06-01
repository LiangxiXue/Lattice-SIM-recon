function validateLatticeSIMStack(stack)
%VALIDATELATTICESIMSTACK Validate an H x W x 5 numeric stack.

if ~isnumeric(stack) || isempty(stack)
    error('LatticeSIM:InvalidInput', 'Image stack must be a non-empty numeric array.');
end
if ndims(stack) ~= 3
    error('LatticeSIM:InvalidInput', 'Image stack must have shape H x W x 5.');
end
if size(stack, 3) ~= 5
    error('LatticeSIM:InvalidFrameCount', 'Expected exactly five frames, got %d.', size(stack, 3));
end
if any(size(stack, 1:2) <= 1)
    error('LatticeSIM:InvalidInput', 'Image frames must be at least 2 x 2 pixels.');
end
if any(~isfinite(double(stack(:))))
    error('LatticeSIM:InvalidInput', 'Image stack contains NaN or Inf values.');
end
end
