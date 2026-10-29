function runDynamicTest
%RUNDYNAMICTEST Dynamic response of Example 3 (CSTR with 3 reactions) to
%a step change in the continuation parameter alpha = FAf, the propylene
%oxide feed flow rate [mol/s].
%
%   Starts from the steady state at alpha = alphaO, steps alpha to
%   alpha(i) at t = tStep(1) and back to alphaO at t = tStep(2), for each
%   value in the list, and plots the reactor temperature T.
%
%   The starting steady state is computed here with fsolve from the
%   initial guess IC, so the test always starts from a true steady state
%   of the current model file (Example3_CSTR_3Reactions.m), even after
%   its fixed parameters (feed temperature, U*A, ...) have been changed.
%   The model function is taken from the JSON file (info.rhs).

IC     = [1618.74834739 14.94513495 904.24894566 458.33792607 383.15112349 4391.85398925 256.35306397];
alphaO = 13.01407695;      % FAf [mol/s] of the starting steady state
alpha  = [14,17.5,18,22];  % FAf [mol/s] during the step
tStep  = [4000 18000];     % step up at tStep(1), back at tStep(2) [s]
tEnd   = 30000;            % [s]

% Look up the reactor temperature column by name instead of a hard-coded
% index, so a mix-up between T (reactor) and Tc (coolant) - as previously
% happened here - cannot silently happen again.
info = loadProblemInfo('Example3_CSTR_3Reactions.json');
rhs  = info.rhs;
tIdx = find(strcmp(info.stateNames,'T'));

% Starting steady state at alphaO
opts = optimset('Display','none','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',5000);
[ICss,~,exitflag] = fsolve(@(x) reshape(rhs([x(:); alphaO]),[],1), IC(:), opts);
if exitflag <= 0
    warning('runDynamicTest:SteadyState', ...
        'Steady state at alpha = %g not found (exitflag %d) - starting from IC as given.', ...
        alphaO, exitflag);
    ICss = IC(:);
else
    fprintf('Starting steady state at alpha = %g: T = %.4f K (IC guess T = %.4f K)\n', ...
        alphaO, ICss(tIdx), IC(tIdx));
end

aLim = [min([alphaO alpha]) max([alphaO alpha])];
aLim = aLim + [-0.1 0.1]*max(diff(aLim), 1);

figure('Name','Dynamic test - Example 3','NumberTitle','off');
for i = 1:length(alpha)
    [t,X] = ode15s(@(tt,XX) DynWrapper(tt,XX,rhs,alphaO,alpha(i),tStep), ...
                   linspace(0,tEnd,1000), ICss);
    T = X(:,tIdx); % reactor temperature, looked up by name (see info above)

    subplot(2,2,i);
    yyaxis left
    plot(t,T);
    ylabel('T [K]');
    xlabel('time [s]');
    yyaxis right
    plot([0 tStep(1) tStep(1) tStep(2) tStep(2) tEnd], ...
         [alphaO alphaO alpha(i) alpha(i) alphaO alphaO],'r');
    ylim(aLim);
    ylabel('\alpha = F_A [mol/s]');
    title(sprintf('F_A: %.4g \\rightarrow %.4g mol/s', alphaO, alpha(i)));
end

function dxdt = DynWrapper(t,X,rhs,alphaO,alphaNew,tStep)
alpha = alphaO;
if t > tStep(1) && t <= tStep(2)
    alpha = alphaNew;
end
dxdt = reshape(rhs([X; alpha]),[],1);
