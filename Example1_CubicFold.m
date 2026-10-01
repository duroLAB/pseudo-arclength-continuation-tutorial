function F = CubicFold(X)
%CUBICFOLD Minimal 1-state, 1-parameter fold/hysteresis benchmark, useful
%as a fast, EXACT sanity check for the continuation engine itself
%(fold locations, 'both'-direction stitching) without any physical model
%complexity getting in the way.
%
%   F = CUBICFOLD(X), X = [x, alpha]
%
%   dx/dtau = alpha + x - x^3
%
%   Steady states: alpha = x^3 - x, an S-shaped (backward-folding) curve
%   with exact, analytically known fold points at
%     x = +1/sqrt(3) = 0.5774, alpha = -2/(3*sqrt(3)) = -0.3849
%     x = -1/sqrt(3) = -0.5774, alpha = +2/(3*sqrt(3)) = +0.3849
%   The outer branches (|x| > 1/sqrt(3)) are stable (df/dx < 0), the
%   middle branch is unstable (df/dx > 0) - the classic cusp/hysteresis
%   picture, just with a single state variable and a closed-form answer.
%
%   Good uses:
%     - run with continuationDirection = 'both' from the interior
%       preset x=0, alpha=0 (the default) and check that the traced curve
%       turns at exactly the two fold x/alpha values above,
%     - run with continuationDirection = 'forward'/'backward' from one of
%       the outer-branch presets and check MinAlpha/MaxAlpha termination,
%     - check that Bif.fold (not Bif.hopf) is where the stability change
%       is flagged, since all eigenvalues here are real.

x = X(1); alpha = X(2);

F = alpha + x - x^3;
