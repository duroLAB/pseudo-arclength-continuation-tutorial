# Continuation Methods Playground

> **This is a learning/testing project, not a production library.** It
> started as a natural-parameter continuation script for one specific
> chemical reactor model and grew into a small, problem-agnostic
> pseudo-arclength continuation engine, used here mainly to understand and
> exercise the method itself - including a handful of deliberately minimal
> benchmark problems (see below) added specifically to stress-test the
> engine's edge cases. Treat the code as a teaching example / sandbox to
> read, run and adapt, not as a maintained, general-purpose toolbox.

A small MATLAB engine for natural-parameter / pseudo-arclength continuation
of steady-state branches of nonlinear systems `F(x, alpha) = 0`, with
stability evaluation and fold/Hopf bifurcation detection along the branch.

It was built around a 3-reaction CSTR (continuous stirred-tank reactor)
model, but the engine itself is problem-agnostic: any system is described
by a small JSON file plus a MATLAB function `F = problem(X)`, `X = [states, parameters]`.

## Requirements

- MATLAB (developed against a recent release; no version-specific syntax
  used, but not tested on old releases).
- **Optimization Toolbox** (`fsolve`, `optimset` are used by the corrector
  step and are required).

## Repository layout

```
engine/     Generic continuation engine - works for ANY problem
              runContinuation.m         - main entry point, run this
              plotContinuationResults.m - results viewer GUI
              continuationResidual.m    - fsolve/Jacobian adapter (internal)
              loadProblemInfo.m         - loads a <problem>.json descriptor
examples/   Problem definitions - one <problem>.m (the equations) + one
            <problem>.json (state/parameter names, initial guesses,
            continuation settings) per problem
              CSTR_3Reactions_Coolant_Test_SS.*  - the original CSTR model
              runDynamicTest.m                   - standalone dynamic step-response
                                                    test for the CSTR model
              UppalRayPooreCSTR.*                - classic fold/hysteresis benchmark
              CubicFold.*                         - minimal exact fold sanity-check
              Brusselator.*                        - exact Hopf-bifurcation benchmark
```

## Quickstart

```matlab
addpath('engine');
addpath('examples');

% Simplest possible example - single state, single parameter, exact
% analytically-known fold points (see below):
runContinuation(loadProblemInfo('CubicFold.json'))
```

This opens a "Continuation Run Setup" dialog (initial guess + continuation
settings, both editable and pre-filled from the JSON), then runs the
continuation, saves the results to a timestamped `.mat` file in the
current directory, and opens the results viewer.

Or try the Brusselator, an exact Hopf-bifurcation benchmark (see below):

```matlab
runContinuation(loadProblemInfo('Brusselator.json'))
```

A saved run can be reopened later, independent of any workspace state:

```matlab
plotContinuationResults('Res_CubicFold_20260927_154212.mat')   % specific file
plotContinuationResults()                                       % file picker
```

`setupPath.m` at the repo root is an optional convenience that does the
two `addpath` calls above for you:

```matlab
setupPath
```

## How it works

Given `F(x, alpha) = 0` (states `x`, a scalar continuation parameter
`alpha`, packed together as one vector), each continuation step:

1. computes the Jacobian of `F` with respect to `[x, alpha]` by central
   finite differences,
2. finds the tangent to the solution curve via SVD (the null-space
   direction),
3. temporarily "ejects" the variable with the largest tangent component,
   fixes it, and solves the remaining (square) system with `fsolve`
   (the corrector step),
4. adapts the step size based on how many Newton iterations `fsolve`
   needed.

This is what lets the engine trace curves around fold points, where a
plain "increase alpha and resolve" approach gets stuck.

After the run, stability is evaluated at every point (eigenvalues of the
Jacobian `dF/dx` at fixed `alpha`), and points where the number of
unstable eigenvalues changes are flagged as fold (real crossing) or Hopf
(complex pair crossing) candidates. This classification is heuristic - it
flags *that* stability changed and whether a complex pair was present, not
necessarily the exact eigenvalue that crossed zero.

`settings.continuationDirection` (`'forward'`, `'backward'`, or `'both'`)
controls which way (or both ways) the branch is traced from the converged
starting point; `'both'` stitches the two half-branches into one ordered
curve.

## What's worth reading this for

A few things in here are specifically useful if you're learning
continuation methods (or testing an implementation of your own), not just
running the CSTR model:

- `CubicFold` is a 1-state, 1-parameter problem with **exact, closed-form**
  fold points - a good reference case to check any continuation
  implementation against, your own or this one.
- `Brusselator` is used here for an **exact, closed-form Hopf point**
  rather than its usual limit-cycle demo role - a reference case for
  testing fold-vs-Hopf classification specifically.
- The single-state `CubicFold` problem is also what caught a real bug in
  this engine: building the initial point as a growing array rather than a
  fixed-size one silently produced a row vector instead of a column vector
  whenever there was exactly one state variable, because that's how MATLAB
  resolves the orientation of a 1x1 array. The 7-state CSTR model could
  never have triggered it. See the "Known limitations" / git history for
  the fix - it's a reasonable illustration of why minimal benchmark
  problems are worth having even when a "real" model already runs fine.

## Example problems

- **`CSTR_3Reactions_Coolant_Test_SS`** - the original 7-state, 3-reaction
  cooled CSTR model, continued with respect to feed temperature (`alpha`).
  Several named initial-guess presets (`X0presets`) are included.
  `runDynamicTest.m` runs a standalone step-response simulation of the same
  model outside the continuation framework.
- **`UppalRayPooreCSTR`** - the classic 2-state, 1-reaction exothermic CSTR
  model (Uppal, Ray & Poore, *Chem. Eng. Sci.* **29** (1974) 967-985 and
  **31** (1976) 205-214), continued with respect to the Damköhler number.
  With the parameters used here (`B=8`, `beta=0.3`, `gamma=20`) it produces
  a classic S-shaped hysteresis curve with two fold points near
  `Da ≈ 0.0662` and `Da ≈ 0.0799`.
  **Note:** those three parameter values were verified numerically (they do
  reliably produce the hysteresis window) rather than copied from a
  specific published numeric table in the original papers - treat them as
  "a parameter set known to work", not as a literal reproduction of a
  published figure.
- **`CubicFold`** - `alpha + x - x^3 = 0`, a single state and single
  parameter with **exact, closed-form** fold points at
  `x = ±1/√3, alpha = ∓2/(3√3)`. No physical meaning - it exists purely as
  a fast, exact correctness check for the engine (fold detection, and
  `continuationDirection = 'both'` in particular).
- **`Brusselator`** - the classic 2-variable autocatalytic oscillator
  (Prigogine & Lefever, 1968), used here as an **exact Hopf-bifurcation**
  benchmark rather than for its usual limit-cycle behavior. The steady
  state `x = a` is constant in the continuation parameter `b`, so the
  traced curve is a straight line by design - the point of this example is
  that `analyzeBifurcations` should classify the stability change at
  `b = 1 + a^2` as `Bif.hopf`, not `Bif.fold`. To see the actual
  oscillation, pick a point in the results viewer's data grid and use the
  "Dynamic Test" button with a parameter step that pushes `b` above the
  Hopf value.

## Known limitations

- Fold/Hopf classification in `analyzeBifurcations` is heuristic (based on
  whether a complex eigenvalue pair is present at a stability change, not
  on tracking the specific eigenvalue that crosses zero).
- `UppalRayPooreCSTR`'s parameters are self-derived (see above), not a
  literal reproduction of a specific published table.

## License

MIT - see `LICENSE`. (Replace the placeholder name in `LICENSE` with your own before publishing.)
