function F = continuationResidual(Xval, problem, ISsquare, eject, constVal)
%CONTINUATIONRESIDUAL Adapter between fsolve / finite differences and the model function.
%
%   F = CONTINUATIONRESIDUAL(Xval, problem, ISsquare, eject, constVal)
%
%   During continuation, one variable (index EJECT) is temporarily
%   removed from the system in each step and fixed to the value
%   CONSTVAL (e.g. the continuation parameter alpha itself, or another
%   variable chosen based on the tangent of the solution curve). The
%   remaining variables are solved by fsolve.
%
%   Inputs:
%     Xval     - row vector of unknowns
%                (if ISsquare==1: without the ejected variable;
%                 if ISsquare==0: full vector including the parameter)
%     problem  - handle to the model function, F = problem(X)
%     ISsquare - 0: Xval is passed to problem() unmodified (used e.g.
%                   when computing the Jacobian by finite differences)
%                1: Xval does not contain the ejected variable; it is
%                   inserted at position EJECT with value CONSTVAL
%                   before calling problem()
%     eject    - index of the variable that is fixed (only if ISsquare==1)
%     constVal - value to which variable EJECT is fixed
%
%   Output:
%     F - residuals of the model function problem()

if ISsquare == 0
    F = problem(Xval);
else
    NoEQ = size(Xval, 2);

    X = zeros(1, NoEQ + 1);
    X(1:eject-1)      = Xval(1:eject-1);
    X(eject+1:NoEQ+1) = Xval(eject:NoEQ);
    X(eject)          = constVal;

    F = problem(X);
end
