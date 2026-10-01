function F = Example2_CSTR_cascade(X)
%EXAMPLE2_CSTR_CASCADE Cascade of two continuous stirred tank reactors
%(CSTRs) with recycle, first-order exothermic reaction.
%
%   F = EXAMPLE2_CSTR_CASCADE(X), X = [x1, x2, theta1, theta2, Da1]
%
%   Source: parametric study of a cascade of two CSTRs with recycle,
%   section 13.4, eqs. (13.24a-d), p. 313 of
%
%     Kubicek, M., Jacak, V., Marek, M.: Numericke algoritmy reseni
%     chemicko-inzenyrskych uloh. Praha: SNTL / Alfa, 1983.
%     (SNTL publication no. 04-614-83)
%
%   Citation:
%     Kubicek M., Jacak V., Marek M. (1983). Numericke algoritmy reseni
%     chemicko-inzenyrskych uloh, Sect. 13.4, Eq. (13.24), p. 313.
%     SNTL / Alfa, Praha. SNTL 04-614-83.
%
%   Dynamic model (4 ODEs):
%
%   dx1/dt     = (1-Lambda)*x2 - x1 + Da1*(1-x1)*exp(theta1/(1+theta1/gamma))
%
%   dtheta1/dt = (1-Lambda)*theta2 - theta1
%                + B*Da1*(1-x1)*exp(theta1/(1+theta1/gamma))
%                - beta1*(theta1 - theta_c1)
%
%   dx2/dt     = (Da1/Da2) * [ x1 - x2 + Da2*(1-x2)*exp(theta2/(1+theta2/gamma)) ]
%
%   dtheta2/dt = (Da1/Da2) * [ theta1 - theta2
%                + B*Da2*(1-x2)*exp(theta2/(1+theta2/gamma))
%                - beta2*(theta2 - theta_c2) ]
%
%   Variables and parameters:
%     x1, x2           - conversions in the 1st and 2nd reactor
%     theta1, theta2   - dimensionless temperatures in the 1st and 2nd reactor
%     Da1, Da2         - Damkohler numbers in the 1st and 2nd reactor
%     beta1, beta2     - heat transfer coefficients
%     theta_c1, theta_c2 - coolant temperatures
%     B                - heat of reaction (exothermicity) parameter
%     Lambda           - recycle parameter
%     gamma            - activation energy parameter
%
%   Notes on this implementation:
%     - State ordering here is X = [x1, x2, theta1, theta2], i.e. it
%       differs from the textbook vector y = (x1, theta1, x2, theta2).
%       F(1) = dx1/dt, F(2) = dx2/dt, F(3) = dtheta1/dt, F(4) = dtheta2/dt.
%     - The continuation parameter is X(5) = Da1, with Da2 = Da1.
%     - kapa = Da2/Da1, so the factor Da1/Da2 in (13.24c,d) is written
%       as division by kapa (kapa = 1 here).
%     - Fixed values: gamma = 1000, B = 22, beta1 = beta2 = 2,
%       theta_c1 = theta_c2 = 0, Lambda = 1.
%     - For steady states, x1 and x2 could be eliminated from (13.24a,c)
%       to get 2 equations, but the full 4-equation system is kept, which
%       is simpler to differentiate and better suited to continuation.


gama=1000;
B=22;
beta1=2;
beta2=2;
fic1=0;
fic2=0;
lam=1;


Da1=X(5);
Da2=Da1;    

kapa=1;

F(1)=(1-lam)*X(2)-X(1)+Da1*(1-X(1))*exp(X(3)/(1+X(3)/gama));
F(3)=(1-lam)*X(4)-X(3)+Da1*B*(1-X(1))*exp(X(3)/(1+X(3)/gama))-beta1*(X(3)-fic1);
F(2)=(X(1)-X(2)+Da2*(1-X(2))*exp(X(4)/(1+X(4)/gama)))/kapa;
F(4)=(X(3)-X(4)+Da2*B*(1-X(2))*exp(X(4)/(1+X(4)/gama))-beta2*(X(4)-fic2))/kapa;