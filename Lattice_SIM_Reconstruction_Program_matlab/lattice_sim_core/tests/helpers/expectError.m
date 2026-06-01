function expectError(fn, expectedIdentifier)
%EXPECTERROR Assert that a function handle throws a specific error.

didThrow = false;
try
    fn();
catch err
    didThrow = true;
    assert(strcmp(err.identifier, expectedIdentifier), ...
        'Expected error %s, got %s: %s', expectedIdentifier, err.identifier, err.message);
end

assert(didThrow, 'Expected error %s, but no error was thrown.', expectedIdentifier);
end
