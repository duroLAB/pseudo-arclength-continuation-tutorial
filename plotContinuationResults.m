function plotContinuationResults(Res, Bif, ProblemFCN, info)
% PLOTCONTINUATIONRESULTS - visualization of continuation results + bifurcation table
%
%   PLOTCONTINUATIONRESULTS(Res, Bif, ProblemFCN, info) - plots a
%   continuation run's results already in the workspace (this is how
%   runContinuation calls it right after a run finishes).
%
%   PLOTCONTINUATIONRESULTS(matFilePath) - opens a previously saved
%   continuation results .mat file (as runContinuation saves after every
%   run, one file per run - see its help) and plots that instead, so a
%   run can be revisited without still having its Res/Bif/info in the
%   workspace, in a fresh MATLAB session, or a different run entirely.
%
%   PLOTCONTINUATIONRESULTS() with no arguments, or PLOTCONTINUATIONRESULTS('')
%   with an empty path, opens a file picker to choose the .mat file.
%
% inputs:
%   Res       - Nop x (X + 1) matrix (last column = parameter alpha), OR
%               the path to a saved results .mat file (char/string), OR
%               omitted/empty to pick that file interactively
%   Bif       - struct from analyzeBifurcations (optional, ignored when
%               Res is a file path - Bif is then loaded from the file)
%   ProblemFCN- handle to the model function F = ProblemFCN([state; alpha]),
%               used by the "Dynamic Test" button
%               (optional, default = @Example3_CSTR_3Reactions;
%               ignored when Res is a file path - resolved from info.rhs
%               in the file instead)
%   info      - problem descriptor (see loadProblemInfo and e.g.
%               Example3_CSTR_3Reactions.json), used to label
%               columns by their real name (cA, cB, ..., alpha) instead
%               of "Col 1", "Col 2", ... (optional, ignored when Res is a
%               file path - loaded from the file instead)

if nargin == 0 || ischar(Res) || isstring(Res)
    if nargin == 0
        matFilePath = '';
    else
        matFilePath = char(Res);
    end
    [Res, Bif, ProblemFCN, info] = loadContinuationResults(matFilePath);
else
    if nargin < 2
        Bif = [];
    end
    if nargin < 3
        ProblemFCN = @Example3_CSTR_3Reactions;
    end
    if nargin < 4
        info = [];
    end
end

[nRows,nCols] = size(Res);

if ~isempty(info)
    colNames = [info.stateNames, info.paramNames];
    alphaColName = info.paramNames{1};
else
    colNames = arrayfun(@(i) sprintf('Col %d',i), 1:nCols, 'UniformOutput', false);
    alphaColName = 'Alpha';
end

colX = 1;
colY = 2;

fig = figure('Position',[200 200 1300 700],...
             'Name','Results Viewer','NumberTitle','off');
ax = axes('Parent',fig,'Position',[0.32 0.36 0.65 0.59]);

% ============= Controls =============
uicontrol(fig,'Style','text','String','X Column','Units','normalized',...
    'Position',[0.05 0.88 0.22 0.04]);
popupX = uicontrol(fig,'Style','popupmenu','String',colNames,...
    'Units','normalized','Position',[0.05 0.84 0.22 0.05],'Value',colX,'Callback',@changeAxes);

uicontrol(fig,'Style','text','String','Y Column','Units','normalized',...
    'Position',[0.05 0.76 0.22 0.04]);
popupY = uicontrol(fig,'Style','popupmenu','String',colNames,...
    'Units','normalized','Position',[0.05 0.72 0.22 0.05],'Value',colY,'Callback',@changeAxes);

btnDyn = uicontrol(fig,'Style','pushbutton','String','Dynamic Test','Units','normalized',...
    'Position',[0.05 0.62 0.22 0.06],'Callback',@runDynamicsTest);

% ============= Bifurcation table =============
tblBif = uitable(fig,'Data',{},'ColumnName',{'Row',alphaColName,'Type'},...
    'Units','normalized','Position',[0.05 0.36 0.22 0.22],'RowName',[]);

% ============= Data grid (all continuation points) =============
% Shows every point, one row per continuation step. Clicking a row
% selects that point - it is highlighted on the plot and becomes the
% starting point for the "Dynamic Test" button. This replaces the old
% "Row" dropdown.
tableData = arrayfun(@(v) sprintf('%.8f', v), Res, 'UniformOutput', false);
rowNames  = cellstr(num2str((1:nRows)'));
tbl = uitable(fig,'Data',tableData,'ColumnName',colNames,'RowName',rowNames,...
    'Units','normalized','Position',[0.05 0.05 0.92 0.28],...
    'CellSelectionCallback',@selectTableRow);

% ============= Store data =============
data.Res = Res;
data.Bif = Bif;
data.ProblemFCN = ProblemFCN;
data.colNames = colNames;
data.colX = colX;
data.colY = colY;
data.currentIndex = 1;
data.ax = ax;
data.popupX = popupX;
data.popupY = popupY;
data.btnDyn = btnDyn;
data.tbl = tbl;
data.tblBif = tblBif;

guidata(fig,data);
updatePlot(fig);

% =====================================================
% ===================== Functions =====================
% =====================================================

function changeAxes(src,~)
    data = guidata(src);
    data.colX = get(data.popupX,'Value');
    data.colY = get(data.popupY,'Value');
    if data.colX == data.colY
        errordlg('X and Y columns must be different.');
        return
    end
    data.currentIndex = 1;
    guidata(src,data);
    updatePlot(fig);
end

function selectTableRow(src,eventdata)
    % Fires when the user clicks a cell in the data grid - selects that
    % row as the current point (highlighted on the plot, and the
    % starting point for "Dynamic Test").
    if isempty(eventdata.Indices)
        return
    end
    data = guidata(src);
    data.currentIndex = eventdata.Indices(1,1);
    updateHighlight(data);
    guidata(src,data);
end

function updatePlot(src)
    data = guidata(src);
    x = data.Res(:,data.colX);
    y = data.Res(:,data.colY);

    cla(data.ax); hold(data.ax,'on');

    if ~isempty(data.Bif)
        stab = data.Bif.stabStatus;
        % stable/unstable
        scatter(data.ax, x(~stab), y(~stab), 40, 'b', 'filled'); % stable
        scatter(data.ax, x(stab),  y(stab),  40, 'r', 'filled'); % unstable

        % stability changes
        if isfield(data.Bif,'stabChange') && ~isempty(data.Bif.stabChange)
            scatter(data.ax, x(data.Bif.stabChange), y(data.Bif.stabChange), 100, 'k', 'o', 'LineWidth',2);
        end

        % fold
        if isfield(data.Bif,'fold') && ~isempty(data.Bif.fold)
            scatter(data.ax, x(data.Bif.fold), y(data.Bif.fold), 80, 'g', 's', 'filled');
        end

        % Hopf
        if isfield(data.Bif,'hopf') && ~isempty(data.Bif.hopf)
            scatter(data.ax, x(data.Bif.hopf), y(data.Bif.hopf), 80, 'y', 'd', 'filled');
        end

        % ============= Update bifurcation table ===========
        BifRows = {};  % <- cell array

        if isfield(data.Bif,'stabChange')
            for i = 1:length(data.Bif.stabChange)
                idx = data.Bif.stabChange(i);
                BifRows(end+1,:) = {idx, data.Res(idx,end), 'StabChange'}; %#ok<AGROW>
            end
        end
        if isfield(data.Bif,'fold')
            for i = 1:length(data.Bif.fold)
                idx = data.Bif.fold(i);
                BifRows(end+1,:) = {idx, data.Res(idx,end), 'Fold'}; %#ok<AGROW>
            end
        end
        if isfield(data.Bif,'hopf')
            for i = 1:length(data.Bif.hopf)
                idx = data.Bif.hopf(i);
                BifRows(end+1,:) = {idx, data.Res(idx,end), 'Hopf'}; %#ok<AGROW>
            end
        end
        set(data.tblBif,'Data',BifRows);

    else
        % fallback - points only
        scatter(data.ax, x, y, 40, 'b', 'filled');
        set(data.tblBif,'Data',{}); % empty table
    end

    % highlight the selected point
    data.hHighlight = plot(data.ax, x(data.currentIndex), y(data.currentIndex),...
        'ko','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','c');

    xlabel(data.ax, data.colNames{data.colX});
    ylabel(data.ax, data.colNames{data.colY});
    title(data.ax, ['Row: ', num2str(data.currentIndex)]);

    data.x = x; data.y = y;
    guidata(src,data);
end

function updateHighlight(data)
    set(data.hHighlight,'XData',data.x(data.currentIndex),...
        'YData',data.y(data.currentIndex));
    title(data.ax,['Row: ', num2str(data.currentIndex)]);
end

function runDynamicsTest(src,~)
    % Run a dynamic test from the currently selected point:
    %   - initial condition = state variables at data.currentIndex
    %   - a dialog asks for the integration time and a step change of
    %     parameter alpha, given as a percentage
    %   - profile: alpha0 (first 1/3 of the time) -> step to
    %     alpha0*(1+p/100) (middle 1/3) -> back to alpha0 (last 1/3),
    %     the same pattern as in runDynamicTest.m
    %   - the result (selected state variable and the alpha profile that
    %     caused it) is plotted in a new figure

    data = guidata(src);
    idx = data.currentIndex;

    NoEQ = size(data.Res,2) - 1;
    IC = data.Res(idx,1:NoEQ)';
    alpha0 = data.Res(idx,end);

    prompt = {sprintf('Current point: row %d, alpha0 = %.6g\nIntegration time [s]:', idx, alpha0), ...
              'Step change of parameter alpha [%]:'};
    answer = inputdlg(prompt, 'Dynamic Test Settings', [1 60], {'30000','10'});
    if isempty(answer)
        return
    end

    tEnd = str2double(answer{1});
    pct  = str2double(answer{2});
    if isnan(tEnd) || tEnd <= 0 || isnan(pct)
        errordlg('Invalid input - enter a positive number for the integration time and a number for the percentage.');
        return
    end

    alphaNew = alpha0 * (1 + pct/100);
    t1 = tEnd/3;
    t2 = 2*tEnd/3;

    [t, X] = ode15s(@(tt,XX) dynStepWrapper(tt,XX,data.ProblemFCN,alpha0,alphaNew,t1,t2), ...
                     linspace(0,tEnd,1000), IC);

    % Which state variable to plot - use the currently selected Y column
    % if it is a state variable (not alpha itself)
    if data.colY <= NoEQ
        plotVarIdx = data.colY;
    else
        plotVarIdx = min(4,NoEQ);
    end

    alphaProfile = alpha0*ones(size(t));
    alphaProfile(t>t1 & t<=t2) = alphaNew;

    varLabel = data.colNames{plotVarIdx};

    figDyn = figure('Name', sprintf('Dynamic test - row %d (alpha %.4g -> %.4g)', idx, alpha0, alphaNew), ...
                     'NumberTitle','off');
    axDyn = axes('Parent',figDyn);
    yyaxis(axDyn,'left');
    plot(axDyn, t, X(:,plotVarIdx), 'b-');
    ylabel(axDyn, varLabel);
    xlabel(axDyn, 't [s]');
    yyaxis(axDyn,'right');
    plot(axDyn, t, alphaProfile, 'r-');
    ylabel(axDyn, '\alpha');
    title(axDyn, sprintf('Initial point = row %d, \\alpha_0 = %.4g, step = %+.3g %%', idx, alpha0, pct));
    legend(axDyn, {varLabel, '\alpha'}, 'Location','best');
end

function dxdt = dynStepWrapper(t, X, ProblemFCN, alpha0, alphaNew, t1, t2)
    if t > t1 && t <= t2
        alpha = alphaNew;
    else
        alpha = alpha0;
    end
    dxdt = ProblemFCN([X; alpha]);
    dxdt = dxdt(:);
end

function [Res, Bif, ProblemFCN, info] = loadContinuationResults(matFilePath)
    % Loads a previously saved continuation results .mat file (as
    % runContinuation saves after every run) and unpacks it into
    % Res/Bif/ProblemFCN/info. Opens a file picker if matFilePath is
    % empty. ProblemFCN is resolved from info.rhs (already a function
    % handle, restored by loadProblemInfo/runContinuation before saving) -
    % it isn't stored separately in the .mat file.
    if isempty(matFilePath)
        [fileName, pathName] = uigetfile('*.mat', 'Select a continuation results file');
        if isequal(fileName,0)
            error('plotContinuationResults:NoFile', 'No results file selected.');
        end
        matFilePath = fullfile(pathName, fileName);
    end

    S = load(matFilePath, 'Res', 'Bif', 'info');
    if ~isfield(S,'Res')
        error('plotContinuationResults:InvalidFile', ...
            '%s does not contain a "Res" variable - is this a continuation results file saved by runContinuation?', matFilePath);
    end
    Res = S.Res;

    if isfield(S,'Bif')
        Bif = S.Bif;
    else
        Bif = [];
    end

    if isfield(S,'info') && ~isempty(S.info)
        info = S.info;
        ProblemFCN = info.rhs;
    else
        info = [];
        ProblemFCN = @Example3_CSTR_3Reactions;
    end
end

end
