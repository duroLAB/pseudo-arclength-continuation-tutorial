function runContinuation(info, settingsOverride)
%RUNCONTINUATION Generic natural-parameter continuation engine - works for any
%problem described by an "info" descriptor. A problem is normally
%defined entirely by a JSON file (see e.g.
%Example3_CSTR_3Reactions.json) loaded with loadProblemInfo - the
%engine itself has no knowledge of any specific problem, not even CSTR.
%
%   RUNCONTINUATION(info) runs continuation for the problem described by info.
%   Before the run starts, an editable, tabbed dialog shows:
%     - "Initial Guess" - state variables + parameters, pre-filled from
%       info.X0presets.(info.activePreset) when info provides one, left
%       blank for manual entry when it doesn't
%     - "Continuation Settings" - which direction(s) to continue in,
%       termination criteria, step size control, Newton step adaptation
%       factors, which two variables to show in the live plot, and
%       whether to run the stability analysis, pre-filled from
%       info.continuationSettings when info provides one, else from the
%       built-in defaults below
%   so the numbers actually used are always visible and can be adjusted
%   without editing any file. A "Save to JSON" button in that dialog
%   writes the edited values back to info.jsonPath (when info was loaded
%   from one via loadProblemInfo), so edits can be kept for next time.
%   RUNCONTINUATION(info, settingsOverride) also overrides any subset of the
%   continuation settings (on top of info.continuationSettings, before
%   the dialog is shown), e.g.:
%       runContinuation(loadProblemInfo('Example3_CSTR_3Reactions.json'), struct('MaxContIter', 2000))
%   RUNCONTINUATION() with no info lets you pick a problem definition JSON file.
%
%   info must have the fields:
%     rhs        - handle to the right-hand-side function, F = rhs(X)
%     stateNames - cell array of state variable names
%     paramNames - cell array of continuation parameter names
%     X0presets  - struct of named starting points, one flat struct per
%                  preset with a field per state/parameter name
%     activePreset - name of the X0presets field to start from
%   info may also have the fields:
%     continuationSettings - struct overriding any subset of the
%                  built-in default run settings below
%     jsonPath   - path of the JSON file info was loaded from (set
%                  automatically by loadProblemInfo); enables the dialog's
%                  "Save to JSON" button
%
%   Natural-parameter / pseudo-arclength continuation:
%   1) in each step, the Jacobian JAC (NoEQ x (NoEQ+1)) of the system
%      with respect to all variables (including alpha) is computed
%      numerically by central differences,
%   2) SVD is used to find the tangent of the solution curve (the last
%      right singular vector), and the "active" variable k = the one
%      with the largest tangent component is chosen - it is temporarily
%      removed (ejected) from the system in this step and fixed to a
%      value (predictor),
%   3) the remaining variables are solved by fsolve (corrector),
%   4) the step size dz is adapted based on the number of Newton
%      iterations used by fsolve.
%   settings.continuationDirection picks which way this is run from the
%   converged starting point: 'forward' (+1), 'backward' (-1, matches the
%   original single hard-coded direction), or 'both' - which runs both
%   directions from the same starting point and stitches them into one
%   continuous, correctly ordered branch (see runContinuationBranch).
%
%   After continuation finishes, stability is evaluated at every point
%   (eigenvalues of the Jacobian df/dx at fixed alpha) and fold/Hopf
%   bifurcation candidates are estimated (analyzeBifurcations) - this can
%   be switched off in the settings below.
%
%   Results are saved to Res_<problem>_<timestamp>.mat (not always the
%   same Res.mat, so older runs are kept instead of overwritten) - Res,
%   Bif, info and settings are all stored, so the run can be reopened
%   later, in a fresh MATLAB session, with:
%       plotContinuationResults('Res_myproblem_20260927_143015.mat')
%   or plotContinuationResults() with no arguments to pick the file
%   interactively.

close all
clc

if nargin < 1 || isempty(info)
    info = getProblemInfo();
end
if nargin < 2
    settingsOverride = struct();
end

%% ---- run settings ----
verbose  = false;  % true = print progress (Newton iterations, dz) to the console
livePlot = true;   % true = plot the branch live while continuation is running

% Continuation settings - direction, termination criteria, step size
% control, and whether to run the stability analysis. Built-in fallback
% defaults below; a problem's info.m can override any subset via
% info.continuationSettings, and settingsOverride overrides on top of
% that. Whatever comes out of this is then shown - and can still be
% edited - in the "Continuation Settings" tab of the dialog below.
settings.MaxContIter                = 1200;
settings.MinAlpha                   = -50;
settings.MaxAlpha                   = 80;
settings.InitialIntegrationStepSize = 0.1;
settings.MinIntegrationStepSize     = 0.02;
settings.MaxIntegrationStepSize     = 150;
settings.runStabilityAnalysis       = true;
settings.stepIncreaseFactor         = 1.3;  % dz multiplier after an easy Newton solve
settings.stepDecreaseFactor         = 0.5;  % dz multiplier after a hard Newton solve
settings.maxNewtonIters              = 3;   % iteration count above which a step counts as "hard"
settings.livePlotXVar                = '';  % variable name for the live-plot X axis, '' = default (alpha)
settings.livePlotYVar                = '';  % variable name for the live-plot Y axis, '' = default (4th state, or the last one)
settings.continuationDirection       = 'backward'; % 'forward' | 'backward' | 'both' - 'backward' matches the original, single hard-coded direction (InitialDirection = -1)

if isfield(info,'continuationSettings')
    infoFields = fieldnames(info.continuationSettings);
    for ii = 1:numel(infoFields)
        settings.(infoFields{ii}) = info.continuationSettings.(infoFields{ii});
    end
end

overrideFields = fieldnames(settingsOverride);
for ii = 1:numel(overrideFields)
    settings.(overrideFields{ii}) = settingsOverride.(overrideFields{ii});
end

ProblemFCN = info.rhs;

%% ---- initial guess X0 + continuation settings dialog ----
% One tabbed dialog, shown before every run:
%   "Initial Guess"         - pre-filled from info.X0presets.(info.activePreset)
%                              when info provides one, blank otherwise
%   "Continuation Settings" - pre-filled from the settings struct built above
% Both tabs are editable. runContinuation() itself still has no
% problem-specific data of its own - it only displays and lets you
% adjust whatever info (and settingsOverride) supplied.
[X, settings] = getRunConfig(info, settings);
if isempty(X)
    disp('Continuation cancelled by the user.');
    return
end

InitialIntegrationStepSize = settings.InitialIntegrationStepSize;
MinIntegrationStepSize     = settings.MinIntegrationStepSize;
MaxIntegrationStepSize     = settings.MaxIntegrationStepSize;
MaxContIter                = settings.MaxContIter;

MinAlpha = settings.MinAlpha;
MaxAlpha = settings.MaxAlpha;

runStabilityAnalysis = settings.runStabilityAnalysis;

% Adaptive continuation step size control based on fsolve convergence -
% now part of settings (info.continuationSettings / the dialog) instead
% of being hard-coded here.
stepIncreaseFactor = settings.stepIncreaseFactor;
stepDecreaseFactor = settings.stepDecreaseFactor;
maxNewtonIters     = settings.maxNewtonIters;

continuationDirection = settings.continuationDirection;

NoEQ = numel(info.stateNames);

% Which columns to show in the live plot during continuation, chosen by
% variable name (settings.livePlotXVar / livePlotYVar) and resolved here
% against [stateNames, paramNames] - the same column order as Res, so
% this works for any problem. Falls back to alpha (X) / the 4th state,
% or the last state if there are fewer than 4 (Y) when the name is blank
% or doesn't match anything in info.
varNames = [info.stateNames, info.paramNames];
livePlotXCol = find(strcmp(varNames, settings.livePlotXVar), 1);
if isempty(livePlotXCol)
    livePlotXCol = NoEQ+1;
end
livePlotYCol = find(strcmp(varNames, settings.livePlotYVar), 1);
if isempty(livePlotYCol)
    livePlotYCol = min(4, NoEQ);
end

ISsquare = 1;
eject = NoEQ+1;
constVal = X(NoEQ+1);

FsolveSettingsStruct=optimset('Display','none','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',5000);

X = X';

% Initial solution of the nonlinear equations (the single starting point
% both directions - when settings.continuationDirection is 'both' - walk
% away from)
[X, ~, exitflag] = fsolve(@continuationResidual, X(1:NoEQ), FsolveSettingsStruct, ProblemFCN, ISsquare, eject, constVal);
if exitflag <= 0
    warning('runContinuation:InitialSolveFailed', 'Initial fsolve solution did not converge (exitflag=%d).', exitflag);
end
% Build X0 explicitly to a fixed (NoEQ+1)x1 column, rather than growing X
% by index (the old "X(NoEQ+1) = constVal;"). For NoEQ==1 (a single-state
% problem, e.g. CubicFold), fsolve's output X is a bare scalar (1x1), and
% MATLAB's default array-growth-by-linear-index on a scalar produces a
% ROW vector, not a column - X(NoEQ+1)=constVal would silently turn X0
% into a 1x2 row instead of a 2x1 column. That corrupted every
% single-state-variable problem downstream: runContinuationBranch's very
% first step (X = X + DX.*dz, with X0 a row and DX a column) then
% broadcasts into a (NoEQ+1)x(NoEQ+1) matrix instead of a (NoEQ+1)x1
% vector, which only surfaces much later as a shape-mismatch error at
% Res(Nop,:) = X. Building X0 explicitly like this is exact for every
% NoEQ, including NoEQ==1, since it never relies on growing an
% ambiguously-shaped array.
X0 = zeros(NoEQ+1,1);
X0(1:NoEQ) = X;
X0(NoEQ+1) = constVal;
X = X0; % converged starting point, column vector, length NoEQ+1

branchOpts = struct( ...
    'MaxContIter',MaxContIter, 'MinAlpha',MinAlpha, 'MaxAlpha',MaxAlpha, ...
    'InitialIntegrationStepSize',InitialIntegrationStepSize, ...
    'MinIntegrationStepSize',MinIntegrationStepSize, ...
    'MaxIntegrationStepSize',MaxIntegrationStepSize, ...
    'stepIncreaseFactor',stepIncreaseFactor, 'stepDecreaseFactor',stepDecreaseFactor, ...
    'maxNewtonIters',maxNewtonIters, 'verbose',verbose, 'livePlot',livePlot);

% Run one or both directions from X0. Each direction is generated by the
% exact same predictor-corrector loop as before (runContinuationBranch) -
% running 'both' does not change a single direction's numbers, it just
% runs the loop twice and stitches the two branches together afterwards.
switch continuationDirection
    case 'forward'
        [Res, ResJac] = runContinuationBranch(X0, 1, ProblemFCN, NoEQ, branchOpts, livePlotXCol, livePlotYCol, varNames);
    case 'backward'
        [Res, ResJac] = runContinuationBranch(X0, -1, ProblemFCN, NoEQ, branchOpts, livePlotXCol, livePlotYCol, varNames);
    case 'both'
        [ResBackward, ResJacBackward] = runContinuationBranch(X0, -1, ProblemFCN, NoEQ, branchOpts, livePlotXCol, livePlotYCol, varNames);
        [ResForward,  ResJacForward ] = runContinuationBranch(X0,  1, ProblemFCN, NoEQ, branchOpts, livePlotXCol, livePlotYCol, varNames);
        % Stitch the two branches back into one continuous, correctly
        % ordered curve: the backward branch was generated walking away
        % from X0, so reverse it to walk back up TO X0, then continue
        % with the forward branch walking away from X0 the other way.
        % (X0 itself is not duplicated - it isn't stored in either
        % branch's Res, matching the original single-direction behavior
        % where the starting point was likewise never stored, only the
        % points after each accepted step.)
        Res    = [flip(ResBackward,1);    ResForward];
        ResJac = cat(1, flip(ResJacBackward,1), ResJacForward);
    otherwise
        error('runContinuation:InvalidDirection', ...
            ['settings.continuationDirection must be ''forward'', ''backward'' ', ...
             'or ''both'' (got ''%s'').'], continuationDirection);
end

if runStabilityAnalysis
    Bif = analyzeBifurcations(ResJac);
else
    Bif = [];
end

% Save results under a name that won't collide with the previous run -
% Res_<problem>_<timestamp>.mat instead of always Res.mat - so older runs
% are kept side by side instead of silently overwritten. info and
% settings are stored alongside Res/Bif, so the run can be reopened later
% - possibly in a different MATLAB session, without still having info in
% the workspace - with plotContinuationResults (see its help).
if isfield(info,'jsonPath') && ~isempty(info.jsonPath)
    [~, problemName] = fileparts(info.jsonPath);
else
    problemName = 'problem';
end
resultsFileName = sprintf('Res_%s_%s.mat', problemName, datestr(now,'yyyymmdd_HHMMSS'));
save(resultsFileName, 'Res', 'Bif', 'info', 'settings');
fprintf('Continuation results saved to %s\n', resultsFileName);

plotContinuationResults(Res, Bif, ProblemFCN, info)


function info = getProblemInfo()
%GETPROBLEMINFO Fallback used when runContinuation() is called with no
%problem descriptor: lets you pick a problem definition JSON file and
%loads it with loadProblemInfo.
%
%   Prefer calling runContinuation(loadProblemInfo('myProblem.json'))
%   directly - this is just a convenience so runContinuation() doesn't
%   have to be given an argument at all.

[fileName, pathName] = uigetfile('*.json', 'Select a problem definition (JSON)');

if isequal(fileName,0)
    error('runContinuation:NoProblemInfo', ...
        ['runContinuation() needs a problem descriptor. Call it as runContinuation(info), ', ...
         'e.g. runContinuation(loadProblemInfo(''Example3_CSTR_3Reactions.json'')).']);
end

info = loadProblemInfo(fullfile(pathName, fileName));


function [Res, ResJac] = runContinuationBranch(X0, direction, ProblemFCN, NoEQ, opts, livePlotXCol, livePlotYCol, varNames)
%RUNCONTINUATIONBRANCH Runs the predictor-corrector continuation loop in
%one direction (DIRECTION = +1 or -1) starting from the already-converged
%point X0 (column vector, length NoEQ+1). This is the exact same
%predictor-corrector loop runContinuation always ran - factored out here
%so it can be run twice (once per direction, from the same X0) when
%settings.continuationDirection is 'both'.
%
%   Returns Res (one row per accepted continuation step) and ResJac (the
%   Jacobian used to compute each step's predictor, i.e. the Jacobian at
%   the point the step was taken FROM - one index "behind" Res, same
%   convention the algorithm always used). X0 itself is not included in
%   Res, matching the original behavior, where the starting point was
%   likewise never stored, only the points after each accepted step.

MaxContIter                = opts.MaxContIter;
MinAlpha                   = opts.MinAlpha;
MaxAlpha                   = opts.MaxAlpha;
InitialIntegrationStepSize = opts.InitialIntegrationStepSize;
MinIntegrationStepSize     = opts.MinIntegrationStepSize;
MaxIntegrationStepSize     = opts.MaxIntegrationStepSize;
stepIncreaseFactor         = opts.stepIncreaseFactor;
stepDecreaseFactor         = opts.stepDecreaseFactor;
maxNewtonIters              = opts.maxNewtonIters;
verbose                     = opts.verbose;
livePlot                    = opts.livePlot;

ISsquare = 1;
FsolveSettingsStruct = optimset('Display','none','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',5000);

dz = InitialIntegrationStepSize;
X = X0;

% Preallocate arrays for computational efficiency
DX   = zeros(NoEQ+1,1);
beta = zeros(NoEQ+1,1);
N    = ones(NoEQ+1,1)*direction;
JAC  = zeros(NoEQ, NoEQ+1);
LFH  = zeros(NoEQ, NoEQ);

Res    = zeros(MaxContIter, NoEQ+1);
ResJac = zeros(MaxContIter, NoEQ, NoEQ);

Nop = 0;

for cyk=1:MaxContIter

    delta = sqrt(eps);

    FF = continuationResidual(X,ProblemFCN,0,0,0);
    for g = 1:NoEQ+1
        nF_plus  = X;
        nF_minus = X;

        h = delta * max(1, abs(X(g)));

        nF_plus(g)  = nF_plus(g)  + h;
        nF_minus(g) = nF_minus(g) - h;

        F_plus  = continuationResidual(nF_plus,ProblemFCN,0,0,0);
        F_minus = continuationResidual(nF_minus,ProblemFCN,0,0,0);

        JAC(:,g) = (F_plus - F_minus) / (2*h);
    end

    % Choose k as the variable with the largest component in the tangent
    % vector (null space of JAC) for stable continuation
    [~,~,V] = svd(JAC);
    t = V(:,end);
    t = t / norm(t);
    [~,k] = max(abs(t));

    LFH(:,1:k-1) = JAC(:,1:k-1);
    LFH(:,k:NoEQ) = JAC(:,k+1:NoEQ+1);

    beta = zeros(NoEQ+1,1);
    freeIdx = [1:k-1, k+1:NoEQ+1];

    beta(freeIdx) = -(LFH \ JAC(:,k));
    beta(k) = 1;

    DXk = N(k) / norm(beta);
    DX(k) = DXk;
    DX([1:k-1, k+1:end]) = beta([1:k-1, k+1:end]) * DXk;

    N(DX>0)  = 1;
    N(DX<0)  = -1;
    N(DX==0) = 0;

    X = X+DX.*dz;

    eject = k;

    constVal = X(eject);
    keepIdx = [1:eject-1, eject+1:NoEQ+1];
    Xstv = X(keepIdx)';

    [Xstv, ~, exitflag, output] = fsolve(@continuationResidual, Xstv, FsolveSettingsStruct, ProblemFCN, ISsquare, eject, constVal);
    if exitflag <= 0
        warning('runContinuation:StepSolveFailed', 'Step %d: fsolve did not converge (exitflag=%d).', cyk, exitflag);
    end

    % Adaptive continuation step size control based on fsolve convergence
    if output.iterations >= maxNewtonIters
        dz = dz * stepDecreaseFactor;
    else
        dz = dz * stepIncreaseFactor;
    end

    if verbose
        fprintf('step %d: Newton iter = %d, dz = %.4g\n', cyk, output.iterations, dz);
    end

    X(1:eject-1)      = Xstv(1:eject-1);
    X(eject+1:NoEQ+1) = Xstv(eject:NoEQ);
    X(eject)          = constVal;

    if dz < MinIntegrationStepSize
        dz = MinIntegrationStepSize;
    end
    if dz > MaxIntegrationStepSize
        dz = MaxIntegrationStepSize;
    end

    Nop = Nop+1;

    Res(Nop,:) = X;
    ResJac(Nop,:,:) = JAC(1:NoEQ,1:NoEQ);

    if livePlot
        plot(Res(1:Nop,livePlotXCol),Res(1:Nop,livePlotYCol),'-+');
        xlabel(varNames{livePlotXCol}); ylabel(varNames{livePlotYCol});
        drawnow;
    end

    if X(NoEQ+1) > MaxAlpha
        break
    end
    if X(NoEQ+1) < MinAlpha
        break
    end
end

Res    = Res(1:Nop,:);
ResJac = ResJac(1:Nop,:,:);


function [X, settings] = getRunConfig(info, settings)
%GETRUNCONFIG Tabbed, editable dialog shown before every run.
%   Tab "Initial Guess"         - one numeric field per entry of
%                                  [info.stateNames, info.paramNames].
%                                  If info.X0presets is available, a
%                                  dropdown lets you pick which preset
%                                  pre-fills the fields (default
%                                  info.activePreset); the fields start
%                                  blank if info has no X0presets.
%   Tab "Continuation Settings"  - which direction(s) to continue in
%                                  (dropdown: forward / backward / both),
%                                  termination criteria, step size
%                                  control, Newton step adaptation
%                                  factors, which two variables the live
%                                  plot shows (dropdowns over
%                                  [info.stateNames, info.paramNames]),
%                                  and whether to run the stability
%                                  analysis, pre-filled from the settings
%                                  struct passed in.
%   Both tabs are editable before the run starts. "Run" uses the current
%   field values for this run only. "Save to JSON" additionally writes
%   them back to info.jsonPath (the selected preset and
%   continuationSettings), so they are still there next time - it does
%   nothing to the run itself and the dialog stays open.
%
%   Returns the flat vector X (in the order given by
%   [info.stateNames, info.paramNames]) and the (possibly edited)
%   settings struct, or X = [] if the dialog is cancelled.

varNames = [info.stateNames, info.paramNames];
nVar = numel(varNames);

hasPresets = isfield(info,'X0presets') && ~isempty(fieldnames(info.X0presets));
if hasPresets
    presetNames = fieldnames(info.X0presets);
    activeIdx = find(strcmp(presetNames, info.activePreset), 1);
    if isempty(activeIdx)
        activeIdx = 1;
    end
else
    presetNames = {};
    activeIdx = 1;
end

rowH = 26;
guessRows    = nVar + hasPresets;
settingsRows = 22;  % 5 section headers + 9 edit fields + 1 (tall) checkbox + 3 dropdowns, with margin
tabH = max(guessRows, settingsRows)*rowH + 20;
figW = 380;
figH = tabH + 100;

fig = figure('Name','Continuation Run Setup','NumberTitle','off', ...
    'MenuBar','none','ToolBar','none','WindowStyle','modal', ...
    'Position',[450 150 figW figH]);
set(fig,'CloseRequestFcn',@(src,~) cancelRun(src));

tg = uitabgroup('Parent',fig,'Units','pixels','Position',[10 60 figW-20 tabH]);
tabGuess    = uitab('Parent',tg,'Title','Initial Guess');
tabSettings = uitab('Parent',tg,'Title','Continuation Settings');

% ============= Tab: Initial Guess =============
y = tabH - 30;

if hasPresets
    uicontrol(tabGuess,'Style','text','String','Preset:', ...
        'Units','pixels','Position',[10 y 80 20],'HorizontalAlignment','left');
    popupPreset = uicontrol(tabGuess,'Style','popupmenu','String',presetNames, ...
        'Units','pixels','Position',[95 y+2 200 24],'Value',activeIdx, ...
        'Callback',@fillFromPreset);
    y = y - rowH;
else
    popupPreset = [];
end

editHandles = gobjects(nVar,1);
for i = 1:nVar
    uicontrol(tabGuess,'Style','text','String',varNames{i}, ...
        'Units','pixels','Position',[10 y 80 20],'HorizontalAlignment','left');
    if hasPresets
        valStr = num2str(info.X0presets.(presetNames{activeIdx}).(varNames{i}), '%.10g');
    else
        valStr = '';
    end
    editHandles(i) = uicontrol(tabGuess,'Style','edit','String',valStr, ...
        'Units','pixels','Position',[95 y+2 200 24]);
    y = y - rowH;
end

% ============= Tab: Continuation Settings =============
y = tabH - 30;

% ---- Direction ----
uicontrol(tabSettings,'Style','text','String','Direction','FontWeight','bold', ...
    'Units','pixels','Position',[10 y 300 20],'HorizontalAlignment','left');
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Continue:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
directionOptions = {'forward','backward','both'};
directionIdx = find(strcmp(directionOptions, settings.continuationDirection), 1);
if isempty(directionIdx)
    directionIdx = 2; % default: backward, matches the original single hard-coded direction
end
popupDirection = uicontrol(tabSettings,'Style','popupmenu','String',directionOptions, ...
    'Units','pixels','Position',[190 y+2 110 24],'Value',directionIdx);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Termination criteria','FontWeight','bold', ...
    'Units','pixels','Position',[10 y 300 20],'HorizontalAlignment','left');
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Max. continuation steps:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMaxContIter = uicontrol(tabSettings,'Style','edit','String',num2str(settings.MaxContIter), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Min. alpha:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMinAlpha = uicontrol(tabSettings,'Style','edit','String',num2str(settings.MinAlpha), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Max. alpha:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMaxAlpha = uicontrol(tabSettings,'Style','edit','String',num2str(settings.MaxAlpha), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Continuation step size','FontWeight','bold', ...
    'Units','pixels','Position',[10 y 300 20],'HorizontalAlignment','left');
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Initial step size:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editInitStep = uicontrol(tabSettings,'Style','edit','String',num2str(settings.InitialIntegrationStepSize), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Min. step size:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMinStep = uicontrol(tabSettings,'Style','edit','String',num2str(settings.MinIntegrationStepSize), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Max. step size:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMaxStep = uicontrol(tabSettings,'Style','edit','String',num2str(settings.MaxIntegrationStepSize), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Newton step adaptation','FontWeight','bold', ...
    'Units','pixels','Position',[10 y 300 20],'HorizontalAlignment','left');
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Step increase factor:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editStepIncrease = uicontrol(tabSettings,'Style','edit','String',num2str(settings.stepIncreaseFactor), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Step decrease factor:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editStepDecrease = uicontrol(tabSettings,'Style','edit','String',num2str(settings.stepDecreaseFactor), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Max. Newton iterations:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
editMaxNewtonIters = uicontrol(tabSettings,'Style','edit','String',num2str(settings.maxNewtonIters), ...
    'Units','pixels','Position',[190 y+2 110 24]);
y = y - rowH;

chkStability = uicontrol(tabSettings,'Style','checkbox', ...
    'String','Perform stability analysis (eigenvalues, fold/Hopf detection)', ...
    'Units','pixels','Position',[10 y-10 340 40],'Value',settings.runStabilityAnalysis);
y = y - 58; % checkbox is taller than a normal row (spans y-10 down to y-50)

% ============= Live plot variable choice =============
% Which two variables (by name, from [stateNames, paramNames]) the plot
% shown while continuation is running uses for X/Y. Falls back to alpha
% (X) / the 4th state, or the last one if there are fewer than 4 (Y) when
% settings.livePlotXVar/YVar is blank or doesn't match any variable name.
uicontrol(tabSettings,'Style','text','String','Live plot','FontWeight','bold', ...
    'Units','pixels','Position',[10 y 300 20],'HorizontalAlignment','left');
y = y - rowH;

nStatesForDefault = numel(info.stateNames);

defXIdx = find(strcmp(varNames, settings.livePlotXVar), 1);
if isempty(defXIdx)
    defXIdx = nVar; % default: alpha (last entry, single-parameter case)
end
defYIdx = find(strcmp(varNames, settings.livePlotYVar), 1);
if isempty(defYIdx)
    defYIdx = min(4, nStatesForDefault);
end

uicontrol(tabSettings,'Style','text','String','X variable:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
popupLiveX = uicontrol(tabSettings,'Style','popupmenu','String',varNames, ...
    'Units','pixels','Position',[190 y+2 110 24],'Value',defXIdx);
y = y - rowH;

uicontrol(tabSettings,'Style','text','String','Y variable:', ...
    'Units','pixels','Position',[10 y 170 20],'HorizontalAlignment','left');
popupLiveY = uicontrol(tabSettings,'Style','popupmenu','String',varNames, ...
    'Units','pixels','Position',[190 y+2 110 24],'Value',defYIdx);

% ============= Save to JSON / Run / Cancel (shared by both tabs) =============
uicontrol(fig,'Style','pushbutton','String','Save to JSON','Units','pixels', ...
    'Position',[10 15 110 30],'Callback',@saveToJson);
uicontrol(fig,'Style','pushbutton','String','Run','Units','pixels', ...
    'Position',[figW-180 15 80 30],'Callback',@acceptRun);
uicontrol(fig,'Style','pushbutton','String','Cancel','Units','pixels', ...
    'Position',[figW-90 15 80 30],'Callback',@(src,~) cancelRun(ancestor(src,'figure')));

data.varNames        = varNames;
data.editHandles     = editHandles;
data.info            = info;
data.hasPresets      = hasPresets;
data.presetNames     = presetNames;
data.popupPreset     = popupPreset;
data.directionOptions   = directionOptions;
data.popupDirection     = popupDirection;
data.editMaxContIter = editMaxContIter;
data.editMinAlpha    = editMinAlpha;
data.editMaxAlpha    = editMaxAlpha;
data.editInitStep    = editInitStep;
data.editMinStep     = editMinStep;
data.editMaxStep     = editMaxStep;
data.editStepIncrease    = editStepIncrease;
data.editStepDecrease    = editStepDecrease;
data.editMaxNewtonIters  = editMaxNewtonIters;
data.chkStability    = chkStability;
data.popupLiveX      = popupLiveX;
data.popupLiveY      = popupLiveY;
data.resultX         = [];
data.resultSettings  = [];
guidata(fig,data);

uiwait(fig);

data = guidata(fig);
X        = data.resultX;
settings = data.resultSettings;
delete(fig);


function fillFromPreset(src,~)
figHandle = ancestor(src,'figure');
data = guidata(figHandle);
selName = data.presetNames{get(data.popupPreset,'Value')};
preset = data.info.X0presets.(selName);
for i = 1:numel(data.varNames)
    set(data.editHandles(i),'String', num2str(preset.(data.varNames{i}), '%.10g'));
end


function [vals, result, ok] = readFormValues(data)
%READFORMVALUES Read and validate both tabs' current field values.
%   Shared by the "Run" and "Save to JSON" buttons so the two stay in
%   sync. ok is false (and an errordlg already shown) if a field isn't a
%   valid number - vals/result are then unusable.

ok = true;
vals = [];
result = [];

vals = zeros(1,numel(data.varNames));
for i = 1:numel(data.varNames)
    v = str2double(get(data.editHandles(i),'String'));
    if isnan(v)
        errordlg(sprintf('Please enter a valid numeric value for %s.', data.varNames{i}));
        ok = false;
        return
    end
    vals(i) = v;
end

result.continuationDirection      = data.directionOptions{get(data.popupDirection,'Value')};
result.MaxContIter                = round(str2double(get(data.editMaxContIter,'String')));
result.MinAlpha                   = str2double(get(data.editMinAlpha,'String'));
result.MaxAlpha                   = str2double(get(data.editMaxAlpha,'String'));
result.InitialIntegrationStepSize = str2double(get(data.editInitStep,'String'));
result.MinIntegrationStepSize     = str2double(get(data.editMinStep,'String'));
result.MaxIntegrationStepSize     = str2double(get(data.editMaxStep,'String'));
result.stepIncreaseFactor         = str2double(get(data.editStepIncrease,'String'));
result.stepDecreaseFactor         = str2double(get(data.editStepDecrease,'String'));
result.maxNewtonIters              = round(str2double(get(data.editMaxNewtonIters,'String')));
result.runStabilityAnalysis       = logical(get(data.chkStability,'Value'));
result.livePlotXVar                = data.varNames{get(data.popupLiveX,'Value')};
result.livePlotYVar                = data.varNames{get(data.popupLiveY,'Value')};

fn = fieldnames(result);
for ii = 1:numel(fn)
    val = result.(fn{ii});
    if ischar(val)
        continue % continuationDirection/livePlotXVar/YVar are names, not numbers
    end
    if ~islogical(val) && isnan(val)
        errordlg('Please enter valid numeric values for all continuation settings.');
        ok = false;
        return
    end
end


function acceptRun(src,~)
figHandle = ancestor(src,'figure');
data = guidata(figHandle);

[vals, result, ok] = readFormValues(data);
if ~ok
    return
end

data.resultX = vals;
data.resultSettings = result;
guidata(figHandle,data);
uiresume(figHandle);


function saveToJson(src,~)
%SAVETOJSON Write the current field values back to info.jsonPath -
%the selected (or newly named) preset, plus continuationSettings.
%Leaves the dialog open; does not affect Run/Cancel.

figHandle = ancestor(src,'figure');
data = guidata(figHandle);

if ~isfield(data.info,'jsonPath')
    errordlg('This problem was not loaded from a JSON file (no jsonPath) - nothing to save to.');
    return
end

[vals, result, ok] = readFormValues(data);
if ~ok
    return
end

if data.hasPresets
    presetName = data.presetNames{get(data.popupPreset,'Value')};
else
    answer = inputdlg('New preset name:', 'Save Initial Guess', [1 50]);
    if isempty(answer) || isempty(strtrim(answer{1}))
        return
    end
    presetName = matlab.lang.makeValidName(strtrim(answer{1}));
end

presetStruct = struct();
for i = 1:numel(data.varNames)
    presetStruct.(data.varNames{i}) = vals(i);
end

raw = jsondecode(fileread(data.info.jsonPath));
raw.X0presets.(presetName) = presetStruct;
raw.activePreset = presetName;
raw.continuationSettings = result;

fid = fopen(data.info.jsonPath,'w');
fwrite(fid, jsonencode(raw,'PrettyPrint',true));
fclose(fid);

msgbox(sprintf('Saved to %s\n(preset "%s").', data.info.jsonPath, presetName), 'Saved');


function cancelRun(figHandle)
data = guidata(figHandle);
data.resultX = [];
data.resultSettings = [];
guidata(figHandle,data);
uiresume(figHandle);


function Bif = analyzeBifurcations(ResJac, tol)
%ANALYZEBIFURCATIONS Stability evaluation and fold/Hopf candidate
%detection along the continuation branch.
%   ResJac - Nop x NoEQ x NoEQ, Jacobian (df/dx at fixed alpha) at every point
%   tol    - tolerance for a "zero" real part of an eigenvalue (default 1e-6)
%
%   Bif.stabStatus - true/false for every point (true = unstable)
%   Bif.stabChange - indices of points where the number of unstable
%                    eigenvalues changed
%   Bif.fold       - subset of stabChange classified as fold (no complex pair)
%   Bif.hopf       - subset of stabChange classified as Hopf (with a complex pair)
%   Bif.eigOut     - Jacobian eigenvalues at every point (nX x Nop)
%
%   NOTE: the fold/Hopf classification is heuristic - it is based on
%   whether ANY complex pair exists among the eigenvalues at that point,
%   not necessarily the one that actually crossed zero. A more accurate
%   classification would need to track the specific eigenvalue that
%   changes stability.

if nargin < 2
    tol = 1e-6;
end

Nop = size(ResJac,1);
nX  = size(ResJac,2);

nUnstable  = zeros(1,Nop);
maxRe      = zeros(1,Nop);
hasComplex = false(1,Nop);
eigOut     = zeros(nX,Nop);

% --- Eigenvalues at every point ---
for k = 1:Nop
    J  = squeeze(ResJac(k,:,:));
    e  = eig(J);
    re = real(e);
    im = imag(e);

    eigOut(:,k)   = e;
    maxRe(k)      = max(re);
    nUnstable(k)  = sum(re > tol);
    hasComplex(k) = any(abs(im) > tol);
end

stabStatus = nUnstable > 0;

% --- Detection of stability changes / fold / Hopf ---
stabChange = [];
fold = [];
hopf = [];

for k = 2:Nop
    if nUnstable(k) ~= nUnstable(k-1)
        stabChange(end+1) = k; %#ok<AGROW>
        if hasComplex(k)
            hopf(end+1) = k; %#ok<AGROW>
        else
            fold(end+1) = k; %#ok<AGROW>
        end
    end
end

% --- output ---
Bif.stabStatus = stabStatus;
Bif.stabChange = stabChange;
Bif.fold       = fold;
Bif.hopf       = hopf;
Bif.maxRe      = maxRe;
Bif.nUnstable  = nUnstable;
Bif.eigOut     = eigOut;
