function F = Example3_CSTR_3Reactions(X)
%EXAMPLE3_CSTR_3REACTIONS CSTR model with 3 reactions and a cooling
%jacket, for use in continuation with respect to parameter alpha.
%
%   F = EXAMPLE3_CSTR_3REACTIONS(X), where X = [par(1:7), alpha]
%
%   State variable order (par):
%     1: cA   - concentration of A (propylene oxide)   [mol/m3]
%     2: cB   - concentration of B (water)              [mol/m3]
%     3: cC   - concentration of C (propylene glycol)   [mol/m3]
%     4: T    - reactor temperature                     [K]
%     5: Tc   - coolant temperature                      [K]
%     6: cD   - concentration of D (dipropylene glycol) [mol/m3]
%     7: cE   - concentration of E (tripropylene glycol)[mol/m3]
%
%   alpha - continuation parameter, shifts the feed temperature
%           Tf = 273.15 + alpha  [K]
%
%   Output F = [dcA/dt dcB/dt dcC/dt dT/dt dTc/dt dcD/dt dcE/dt]  (in this order)
%
%   Source / inspiration:
%     This example is fully inspired by the model system of
%
%     Svandova, Z., Jelemensky, L., Markos, J., Molnar, A.: Steady States
%     Analysis and Dynamic Simulation as a Complement in the Hazop Study
%     of Chemical Reactors. Process Safety and Environmental Protection
%     83(5) (2005), pp. 463-471. ISSN 0957-5820.
%     https://doi.org/10.1205/psep.04262
%     https://www.sciencedirect.com/science/article/pii/S0957582005712829
%
%   Summary of the paper:
%     The paper couples mathematical modelling and simulation of chemical
%     reactors with the HAZOP (Hazard and Operability) study and proposes
%     a methodology for hazard investigation. The model system is the
%     exothermic hydrolysis of propylene oxide to mono-propylene glycol,
%     with consecutive reactions to higher glycols (di- and tripropylene
%     glycol), in jacket-cooled CSTRs (two in series in the paper). The
%     safety analysis covers multiplicity of steady states and their
%     stability, safe operating conditions, trajectories that move the
%     reactor from one steady state to another, and parametric studies
%     of failures of the reactant feed flow rates (propylene oxide,
%     water) and of the cooling medium.
%     Keywords: CSTR; HAZOP; nonlinear behaviour; multiple steady states;
%     safety analysis; stability analysis.

F = reactorEquations(X(1:7), X(8));


%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PS = reactorEquations(par, alpha)
cA=par(1);  cB=par(2); cC=par(3); T=par(4); Tc=par(5); cD=par(6); cE=par(7);

%% ---- constants ----
%  Propylene oxide   Water           Propylene glycol   Dipropylene glycol  Tripropylene glycol  Methanol
Mh=[58.08e-3	   18.02e-3        76.11e-3          134.17e-3          192.26e-3      32.04e-3]; % kg/mol

roa=[1.4855,           0,          1.015             0.65838            0.44619        2.3267 ];   % kmol/m3
rob=[0.2763            0           0.25071           0.265              0.26591        0.27073];
roc=[482.25            0           676.4             654                677              512.5];
rod=[0.29365           0           0.23083           0.2857             0.24367        0.24713];

ro_water=[17.863 58.606 -95.396 213.89 -141.26]*1000; % mol/m3
Tc_water=647.13; % kelvin

cpa=[167910         276370            0             141370              129010           256040]/1000;  %J/mol/K
cpb=[-696.6	        -2090.1           0              619.5              996.27          -2741.4]/1000;
cpc=[2.45           8.125             0                0                0.1742           14.777]/1000;
cpd=[-0.0021734     -0.01412          0                0                  0            -0.03508]/1000;
cpe=[0              9.37e-6           0                0                  0             3.27e-5]/1000;

cp_propylene_glycol=[265720 1754.6 20442 -281090 0]/1000; %J/mol/K
Tc_propylene_glycol=626; % kelvin

%% ---- input data ----
VR=1; %m3
N_coolant=391; %coolant medium holdup, mol
n_coolant=105; %mol/s

Tc_feed=288.15; % K 288.15 K

% Heat transfer U*A. In the original file the physically meaningful
% values were immediately overwritten by test values (practically
% infinite heat transfer, T approx. equal to Tc). Both variants are kept
% here and switched via useTestHeatTransfer - confirm which regime you
% actually want to use.
useTestHeatTransfer = false;
if useTestHeatTransfer
    U = 1;       % (J s^-1 m^-2 K^-1) - test value
    A = 10000;   % m2                 - test value
else
    U = 1.65e3;  % (J s^-1 m^-2 K^-1) - physical value
    A = 6.7;     % m2                 - physical value
end

%feed
n_methanol=0;             %mol/s
cF=0;

FAf = 11.94444;
FBf = 6.0;

Ff=[FAf, FBf, 0, 0, 0, n_methanol]; %mol/s
Tf = 273.15 + alpha; % Tf depends on the continuation parameter alpha

%A: n_PO (propylene oxide) 43->11.94444
%B: n_H2O (water) 21.6->6 kmol/h

%% ---- kinetic parameters ----
kvoo1=96000;      %m3/mol/s
E1=75362.0;       %J/mol
drHref1=-91360.0; %J/mol; ref. at 25 C

kvoo2=9600;      %m3/mol/s
E2=82899.0;       %J/mol
drHref2=-111000.0; %J/mol; ref. at 25 C

kvoo3=960;      %m3/mol/s
E3=91189;       %J/mol
drHref3= -1.1e5; % -244600.0;  J/mol; ref. at 25 C

%% ---- feed calculations ----
m_mix=sum(Ff.*Mh);                                 % kg/s

ro_feed = getPureComponentMolarDensity(Tf);         % mol/m3 (correlation, temperature clamped to valid range)

Vfi=Ff./(ro_feed); % volumetric flow of each component  % m3/s
Vf=sum(Vfi); % total feed volumetric flow
cf=Ff./Vf;   % feed concentrations                       mol/m3
cAf=cf(1);  cBf=cf(2); cCf=cf(3); cDf=cf(4);  cEf=cf(5); cFf=cf(6);

%% ---- heat exchanged ----
% for T
Q=U*A*(Tc-T);

%% ---- enthalpies and heat capacities ----
% T and Tf at Tref = 298
Tref=298;
cp_T=cpa+cpb*(T)+cpc*(T)^2+cpd*(T)^3+cpe*(T)^4;
tau_propylene_glycol=1-T/Tc_propylene_glycol;
cp_propylene_glycol_T=cp_propylene_glycol(1)+cp_propylene_glycol(2)/tau_propylene_glycol+cp_propylene_glycol(3)*tau_propylene_glycol+cp_propylene_glycol(4)*tau_propylene_glycol^2+cp_propylene_glycol(5)*tau_propylene_glycol^3; %J/kmol/K
cp_T(3)=cp_propylene_glycol_T;
CiCpi=cp_T(1)*cA+cp_T(2)*cB+cp_T(3)*cC+cp_T(4)*cD+cp_T(5)*cE+cp_T(6)*cF;


tau_propylene_glycol_ref=1-Tref/Tc_propylene_glycol;
tau_propylene_glycol=1-T/Tc_propylene_glycol;
h_propylene_glycol=-Tc_propylene_glycol*[cp_propylene_glycol(1)*(tau_propylene_glycol-tau_propylene_glycol_ref)+cp_propylene_glycol(2)*(log(tau_propylene_glycol)-log(tau_propylene_glycol_ref))+cp_propylene_glycol(3)/2*(tau_propylene_glycol^2-tau_propylene_glycol_ref^2)+cp_propylene_glycol(4)/3*(tau_propylene_glycol^3-tau_propylene_glycol_ref^3)+cp_propylene_glycol(5)/4*(tau_propylene_glycol^4-tau_propylene_glycol_ref^4)];

tau_propylene_glycol_feed=1-Tf/Tc_propylene_glycol;
h_propylene_glycol_feed=-Tc_propylene_glycol*[cp_propylene_glycol(1)*(tau_propylene_glycol_feed-tau_propylene_glycol_ref)+cp_propylene_glycol(2)*(log(tau_propylene_glycol_feed)-log(tau_propylene_glycol_ref))+cp_propylene_glycol(3)/2*(tau_propylene_glycol_feed^2-tau_propylene_glycol_ref^2)+cp_propylene_glycol(4)/3*(tau_propylene_glycol_feed^3-tau_propylene_glycol_ref^3)+cp_propylene_glycol(5)/4*(tau_propylene_glycol_feed^4-tau_propylene_glycol_ref^4)];


h_T= cpa*(T-Tref) +cpb/2*(T^2-Tref^2) +cpc/3*(T^3-Tref^3) +cpd/4*(T^4-Tref^4) +cpe/5*(T^5-Tref^5);
h_T(3)=h_propylene_glycol;
h_Tf=cpa*(Tf-Tref)+cpb/2*(Tf^2-Tref^2)+cpc/3*(Tf^3-Tref^3)+cpd/4*(Tf^4-Tref^4)+cpe/5*(Tf^5-Tref^5);
h_Tf(3)=h_propylene_glycol_feed;

% Tc and Tc_feed for the cooling medium (water)
Cp_coolant=cpa(2)+cpb(2)*Tc+cpc(2)*Tc^2+cpd(2)*Tc^3+cpe(2)*Tc^4;
h_coolant= cpa(2)*Tc +cpb(2)/2*Tc^2 +cpc(2)/3*Tc^3 +cpd(2)/4*Tc^4 +cpe(2)/5*Tc^5;
h_coolant_feed=cpa(2)*Tc_feed+cpb(2)/2*Tc_feed^2+cpc(2)/3*Tc_feed^3+cpd(2)/4*Tc_feed^4+cpe(2)/5*Tc_feed^5;

% reaction enthalpy change at Tref = 298
Tref=298;
h_propylene_glycol_Hr=-Tc_propylene_glycol*[cp_propylene_glycol(1)*(tau_propylene_glycol-tau_propylene_glycol_ref)+cp_propylene_glycol(2)*(log(tau_propylene_glycol)-log(tau_propylene_glycol_ref))+cp_propylene_glycol(3)/2*(tau_propylene_glycol^2-tau_propylene_glycol_ref^2)+cp_propylene_glycol(4)/3*(tau_propylene_glycol^3-tau_propylene_glycol_ref^3)+cp_propylene_glycol(5)/4*(tau_propylene_glycol^4-tau_propylene_glycol_ref^4)];
h_Hr= cpa*(T-Tref) +cpb/2*(T^2-Tref^2) +cpc/3*(T^3-Tref^3) +cpd/4*(T^4-Tref^4) +cpe/5*(T^5-Tref^5);
h_Hr(3)=h_propylene_glycol_Hr;
drH1=drHref1+(-h_Hr(1)-h_Hr(2)+h_Hr(3));
drH2=drHref2+(-h_Hr(1)-h_Hr(3)+h_Hr(4));
drH3=drHref3+(-h_Hr(1)-h_Hr(4)+h_Hr(5));
drH=[drH1 drH2 drH3];

%% ---- volumetric flow ----
% for T, cA, cB, cC, cD, cE
ro = getPureComponentMolarDensity(T);               % mol/m3 (correlation, temperature clamped to valid range)

c_total=cA+cB+cC+cD+cE+cF;
x(1)=cA/c_total; x(2)=cB/c_total; x(3)=cC/c_total; x(4)=cD/c_total; x(5)=cE/c_total; x(6)=cF/c_total;

ro_mix_molar=(sum(x./ro))^-1;      %molar mixture density      mol/m3
Mh_mix=sum(x.*Mh);                 %molar mass of the mixture  kg/mol
ro_mix_mass=ro_mix_molar*Mh_mix;   %mass density of the mixture kg/m3
V=m_mix/ro_mix_mass;                %volumetric flow, from mass balance

%% ---- kinetics ----
% for T, cA, cB, cC, cD, cE
kv1=kvoo1*exp(-E1/8.314/T);
ksi1=kv1*cA*cB;

kv2=kvoo2*exp(-E2/8.314/T);
ksi2=kv2*cA*cC;

kv3=kvoo3*exp(-E3/8.314/T);
ksi3=kv3*cA*cD;

%% ---- differential equations ----
dcAdt=(Vf*cAf-V*cA-VR*ksi1-VR*ksi2-VR*ksi3)/VR;
dcBdt=(Vf*cBf-V*cB-VR*ksi1)/VR;
dcCdt=(Vf*cCf-V*cC+VR*ksi1-VR*ksi2)/VR;
dcDdt=(Vf*cDf-V*cD+VR*ksi2-VR*ksi3)/VR;
dcEdt=(Vf*cEf-V*cE+VR*ksi3)/VR;
dTdt =(sum((Vf.*cf).*(h_Tf-h_T))+VR*(-drH(1))*ksi1+VR*(-drH(2))*ksi2+VR*(-drH(3))*ksi3+Q)/VR/CiCpi;
dTcdt=(n_coolant*(h_coolant_feed-h_coolant)-Q)/N_coolant/Cp_coolant;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PS=[dcAdt dcBdt dcCdt dTdt dTcdt dcDdt dcEdt];

%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%%C%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ro = getPureComponentMolarDensity(T)
% A - propylene oxide, valid range 161-482 K
% B - water, valid range 273 - 647 K
% C - propylene glycol, valid range 213 - 676 K
% D - dipropylene glycol, valid range 233 - 654 K
% E - tripropylene glycol, valid range 228 - 677 K
% F - methanol, valid range 175 - 512.50 K
% test
%     300         450
% A 1.41341E+04 9614.541
% B 5.51762E+04 5.38866E+04
% C 1.35548E+04 1.18573E+04
% D 7572.06 6437.22
% E 5291.428 4629.935
% F 2.45876E+04 1.86885E+04
Tthreshold = [161 482;
              273 647;
              213 676;
              233 654;
              228 677;
              175 512.5];

T = ones(6,1) * T';
T = min(max(T, Tthreshold(:,1)), Tthreshold(:,2));

roa=[1.4855            0          1.015             0.65838            0.44619        2.3267]';   % kmol/m3
rob=[0.2763            0           0.25071           0.265              0.26591        0.27073]';
roc=[482.25            0           676.4             654                677              512.5]';
rod=[0.29365           0           0.23083           0.2857             0.24367        0.24713]';

ro_water=[17.863 58.606 -95.396 213.89 -141.26]*1000; % mol/m3
Tc_water=647.13; % kelvin
tau_water=1-(T(2))/(Tc_water);
ro_water=ro_water(1)+ro_water(2)*tau_water^0.35+ro_water(3)*tau_water^(2/3)+ro_water(4)*tau_water+ro_water(5)*tau_water^(4/3);   % mol/m3
ro=(roa./(rob.^(1+(1-T./roc).^rod)))*1000;
ro(2)=ro_water;   % mol/m3
ro = ro';
