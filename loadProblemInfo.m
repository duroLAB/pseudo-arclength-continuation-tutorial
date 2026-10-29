function info = loadProblemInfo(jsonPath)
%LOADPROBLEMINFO Load a problem descriptor from a JSON file.
%
%   info = LOADPROBLEMINFO(jsonPath)
%
%   This is the one generic loader used for ANY problem - a new problem
%   needs only a new JSON file (see e.g.
%   Example3_CSTR_3Reactions.json), not a new .m descriptor file.
%   The JSON is expected to have the fields:
%     rhs        - name (string) of the right-hand-side function,
%                  F = rhs(X), with X = [states, parameters]
%     stateNames - array of state variable names
%     paramNames - array of continuation parameter names
%     X0         - the starting point, one flat object with a field per
%                  state/parameter name
%
%   One JSON file = one continuation task (one parameter, one starting
%   point, one set of continuationSettings). To continue in a different
%   physical parameter, edit the rhs file and use a separate JSON file
%   with its own X0 and settings.
%
%   Older JSON files with X0presets/activePreset are still accepted: the
%   active preset (or the first one) is used as X0 and a warning is shown.
%     continuationSettings - default run settings (termination criteria,
%                  step size control, stability analysis on/off)
%
%   info is returned exactly as runContinuation() expects it, with two
%   additions:
%     info.rhs      - resolved from the JSON's function-name string to an
%                      actual function handle via str2func
%     info.jsonPath - the file this was loaded from, so the "Continuation
%                      Run Setup" dialog can save edits back to it

raw = jsondecode(fileread(jsonPath));

info = raw;
info.rhs = str2func(raw.rhs);
info.jsonPath = jsonPath;

% jsondecode turns a JSON array into a column cell array (Nx1), not the
% row (1xN) shape that {'a','b',...} gives in hand-written MATLAB code -
% and for a single-element array the "column" is just 1x1, which still
% doesn't horzcat against a 7x1. Force both to row vectors here, once,
% so [info.stateNames, info.paramNames] (used by runContinuation,
% plotContinuationResults) always works regardless of how many of each
% there are.
info.stateNames = reshape(info.stateNames, 1, []);
info.paramNames = reshape(info.paramNames, 1, []);

% Backward compatibility: older files had several named starting points
% (X0presets + activePreset) instead of a single X0.
if ~isfield(info,'X0') && isfield(info,'X0presets') && ~isempty(fieldnames(info.X0presets))
    presetNames = fieldnames(info.X0presets);
    useName = presetNames{1};
    if isfield(info,'activePreset') && any(strcmp(presetNames, info.activePreset))
        useName = info.activePreset;
    end
    info.X0 = info.X0presets.(useName);
    warning('loadProblemInfo:LegacyPresets', ...
        ['%s uses the old X0presets/activePreset format - preset "%s" is used as X0. ', ...
         'Press "Save to JSON" in the setup dialog to convert the file to the new format.'], ...
        jsonPath, useName);
end
if isfield(info,'X0presets')
    info = rmfield(info,'X0presets');
end
if isfield(info,'activePreset')
    info = rmfield(info,'activePreset');
end
