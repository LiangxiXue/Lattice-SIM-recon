function result = reconstructLatticeSIM(input, params)
%RECONSTRUCTLATTICESIM Reconstruct 2D five-frame Lattice-SIM data.

if nargin < 2
    params = struct();
end
params = defaultLatticeSIMParams(params);
validateLatticeSIMParams(params, true);

[rawStack, inputMetadata] = readLatticeSIMInput(input);
[stack, normalizationInfo] = normalizeSIMFrames(rawStack, params);
params.normalizationInfo = normalizationInfo;

bands = separateLatticeBands(stack, params);
[carriers, carrierDiagnostics] = estimateLatticeCarrier(bands, params);
otf = buildLatticeOTF(size(stack, 1), size(stack, 2), params);
[SIM, combineDiagnostics] = combineLatticeSpectrum(bands, carriers, otf, params);

result.WF = upsampleWidefield(bands.C0, 2);
result.SIM = SIM;
result.params = params;
result.diagnostics = carrierDiagnostics;
result.diagnostics.input = inputMetadata;
result.diagnostics.normalization = normalizationInfo;
result.diagnostics.otf = rmfield(otf, {'values'});
result.diagnostics.combine = combineDiagnostics;
end

function image = upsampleWidefield(image, scale)
image = imresize(abs(image), scale);
end
