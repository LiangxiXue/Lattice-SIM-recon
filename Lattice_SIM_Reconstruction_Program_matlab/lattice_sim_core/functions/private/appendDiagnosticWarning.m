function diagnostics = appendDiagnosticWarning(diagnostics, message)
%APPENDDIAGNOSTICWARNING Append a warning string to diagnostics.

if ~isfield(diagnostics, 'warnings') || isempty(diagnostics.warnings)
    diagnostics.warnings = {};
end
diagnostics.warnings{end + 1} = message;
end
