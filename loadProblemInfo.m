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
%     X0presets  - object of named starting points, one flat object per
%                  preset with a field per state/parameter name
%     activePreset - name of the X0presets field to start from
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
