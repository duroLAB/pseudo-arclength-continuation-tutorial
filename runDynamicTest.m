function runDynamicTest
%RUNDYNAMICTEST Dynamic response of the CSTR model to a step change in
%parameter alpha.
%
%   Simulates a transition from the steady-state value alphaO to
%   alpha(i) at t=4000 s and back to alphaO at t=18000 s, for each value
%   of alpha in the list.

IC = [1618.74834739	14.94513495	904.24894566	458.33792607	383.15112349	4391.85398925	256.35306397];
alphaO = 13.01407695;
alpha = [14,17.5,18,22];

% Look up the reactor temperature column by name instead of a hard-coded
% index, so a mix-up between T (reactor) and Tc (coolant) - as previously
% happened here - cannot silently happen again.
info = loadProblemInfo('Example3_CSTR_3Reactions.json');
tIdx = find(strcmp(info.stateNames,'T'));

for i = 1:length(alpha)
    [t,X] = ode15s(@DynWrapper,linspace(0,30000,1000),IC,[],alpha(i));
    T = X(:,tIdx); % reactor temperature, looked up by name (see info above)
    subplot(2,2,i);
    yyaxis left
    ylabel('T [K]');
    xlabel('time [s]');
    plot(t,T);
    yyaxis right
    plot([0 4e3 4e3 18e3 18e3 30e3],[alphaO,alphaO,alpha(i),alpha(i),alphaO,alphaO],'r');
    ylim([12 25]);  % limit right axis
    ylabel('\alpha - Fa');
end

function dxdt = DynWrapper(t,X,alphaNew)
alphaOP = 13.01407695;
alpha = alphaOP;
if(t>4000)
    alpha = alphaNew;
end
if(t>18000)
    alpha = alphaOP;
end

dxdt = Example3_CSTR_3Reactions([X; alpha]);
dxdt = dxdt';
