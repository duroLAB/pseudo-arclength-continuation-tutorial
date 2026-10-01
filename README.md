# Continuation Algorithm — Documentation

MATLAB toolbox for one-parameter continuation of steady states.

## Contents

- [Overview](#overview)
- [Origins: DERPAR by Milan Kubíček](#origins-derpar-by-milan-kubíček)
- [Numerics in brief](#numerics-in-brief)
- [Files](#files)
- [Running it](#running-it)
- [Input parameters](#input-parameters)
- [Outputs](#outputs)
- [Adding a new problem](#adding-a-new-problem)
- [Examples](#examples)
- [References](#references)

## Overview

This MATLAB package traces a branch of steady states of a nonlinear system f(x, α) = 0 as one parameter α changes, evaluates their stability and marks turning points (folds) and Hopf points. Three examples are bundled, from a one-equation fold benchmark to a seven-state CSTR model (see [Examples](#examples)).

The core is problem-independent: a new problem = one .m function with the right-hand sides + one JSON file describing it. The core code does not change. This document is written for a user who wants to run the algorithm on their own problem.

## Origins: DERPAR by Milan Kubíček

This toolbox is explicitly inspired by **DERPAR**, the continuation code by Professor Milan Kubíček (University of Chemistry and Technology, Prague), published as [Algorithm 502: Dependence of Solution of Nonlinear Systems on a Parameter](https://dl.acm.org/doi/abs/10.1145/355666.355675), ACM Transactions on Mathematical Software 2(1), 1976.

Professor Kubíček's idea is elegant and remarkably durable. Differentiate f(x, α) = 0 with respect to arc length, and let the variable that changes fastest along the curve serve as the local parameter. With that one step, the method follows a solution branch through turning points, where plain parameter stepping in α breaks down. Half a century later the scheme is still the backbone of practical continuation in chemical engineering, and it is the backbone of this toolbox too. Together with M. Marek he also wrote *Computational Methods in Bifurcation Theory and Dissipative Structures* (Springer, 1983), a standard reference in the field.

The toolbox keeps DERPAR's core idea and changes mainly the implementation. The local parameter k is taken from the tangent obtained by SVD of the full Jacobian, the predictor is a single Euler step along that tangent, and the corrector is MATLAB's `fsolve` instead of a Newton iteration coded in the routine. The Jacobian is always computed by central differences, the step length follows the number of `fsolve` iterations, and the branch can be traced in both directions from the starting point. Stability analysis with fold/Hopf detection, a JSON problem file and a graphical interface are additions on top of the original algorithm.

## Numerics in brief

The engine follows the solution curve of f(z) = 0, z = [x, α], with NoEQ equations and NoEQ + 1 unknowns. It uses a predictor–corrector scheme with local parametrization: in every step, the variable that changes fastest along the curve is held fixed. That variable can be α or any state, so the branch passes through turning points.

One continuation step:

1. **Jacobian.** J = ∂f/∂z (NoEQ × (NoEQ+1)) by central differences, h = √eps · max(1, |zᵢ|).
2. **Tangent.** t = last right singular vector of J (its null space). The index k = argmax |tᵢ| becomes the local parameter.
3. **Predictor.** Solve for the direction β below, then step z ← z + dz · Δz. The sign of each component is carried over from the previous step, so the branch does not turn back.
4. **Corrector.** zₖ stays fixed; `fsolve` solves the remaining NoEQ unknowns (TolFun = TolX = 1e-12, MaxFunEvals = 5000) through `continuationResidual`.
5. **Step control.** If `fsolve` needed ≥ maxNewtonIters iterations, dz is multiplied by stepDecreaseFactor, otherwise by stepIncreaseFactor. dz is then clamped to [MinIntegrationStepSize, MaxIntegrationStepSize].
6. **Stop** when α leaves [MinAlpha, MaxAlpha] or after MaxContIter steps.

$$
J_{:,-k}\,\beta_{-k} = -J_{:,k}, \qquad \beta_k = 1, \qquad \Delta z = \pm\frac{\beta}{\lVert \beta \rVert}
$$

**Stability.** After the run, the eigenvalues of ∂f/∂x (α fixed) are computed at every point. A point is unstable if any eigenvalue has Re λ > 1e-6. Where the number of unstable eigenvalues changes between two neighbouring points, the point is flagged as a fold if no complex pair is present, otherwise as Hopf. This classification is heuristic: it does not track which eigenvalue actually crossed the imaginary axis.

## Files

For a new problem you write only a model function and its JSON file; the core files are not edited.

| File | Kind | What it does |
| --- | --- | --- |
| `runContinuation.m` | core, **main** | Continuation engine. Shows the setup dialog, traces the branch, runs the stability analysis (`analyzeBifurcations`), saves the results and opens the Results Viewer. |
| `loadProblemInfo.m` | core | Reads a problem JSON into the `info` struct and turns the `rhs` name into a function handle. |
| `continuationResidual.m` | core | Adapter between `fsolve` and the model: puts the fixed (ejected) variable back into the vector before calling `rhs`. |
| `plotContinuationResults.m` | core | Results Viewer: plot with selectable X/Y columns, table of all points, table of fold/Hopf points, "Dynamic Test" button. Also reopens a saved results file. |
| `Example1_CubicFold.m` / `.json` | example | One-equation fold benchmark with known exact answer. |
| `Example2_CSTR_cascade.m` / `.json` | example | Two CSTRs in series with recycle, from Kubíček, Jacák and Marek (1983). |
| `Example3_CSTR_3Reactions.m` / `.json` | example | CSTR with three reactions and a cooling jacket, after Švandová et al. (2005). |
| `runDynamicTest.m` | optional script | Stand-alone dynamic simulation (`ode15s`) of Example 3 for step changes of α; values are hard-coded. |

## Running it

The entry point is `runContinuation.m`. Put all files in one folder (or on the MATLAB path) and call:

```matlab
runContinuation(loadProblemInfo('Example3_CSTR_3Reactions.json'))

% override some settings for this run only
runContinuation(loadProblemInfo('Example3_CSTR_3Reactions.json'), struct('MaxContIter', 2000))

% no argument: a file picker asks for the problem JSON
runContinuation()
```

Requires MATLAB with the Optimization Toolbox (`fsolve`). The Dynamic Test uses `ode15s` from core MATLAB.

What happens after the call:

1. **Setup dialog** with two tabs. *Initial Guess*: a preset dropdown and one field per variable. *Continuation Settings*: all settings from the table below. Buttons: **Run** uses the values for this run only; **Save to JSON** writes them back to the JSON file and keeps the dialog open; **Cancel** stops.
2. The first point is corrected by `fsolve` with α fixed, so the initial guess only has to be close to a steady state.
3. The branch is traced forward, backward or both ways from that point (setting `continuationDirection`). With `both`, the two halves are joined into one ordered branch.
4. **Live plot** of the branch during the run (axes set by livePlotXVar / livePlotYVar).
5. Stability analysis, the results file is saved, and the **Results Viewer** opens.

## Input parameters

All inputs live in the problem JSON; every value can still be changed in the setup dialog before a run. Settings are applied in three layers: built-in defaults in `runContinuation.m` → `continuationSettings` from the JSON → the optional `settingsOverride` struct.

### The model function (`rhs`)

- Signature `F = myModel(X)`, where X is a row vector `[states, alpha]` in the order of `stateNames` then `paramNames`.
- Returns NoEQ values, one per state: the time derivatives dx/dt. Steady states are F = 0; the Dynamic Test integrates the same function with `ode15s`.
- Exactly one continuation parameter is supported, and it is the last element of X.

### JSON fields

| Field | Type | Meaning |
| --- | --- | --- |
| `rhs` | string | Name of the model function (its file name, without `.m`). |
| `stateNames` | array of strings | State names, in the order used in X. Their count sets NoEQ. |
| `paramNames` | array of 1 string | Name of the continuation parameter, e.g. `"alpha"`. |
| `X0presets` | object | Named starting points; each is an object with one value per state and parameter name. |
| `activePreset` | string | Preset pre-selected in the dialog. |
| `continuationSettings` | object | Run settings, table below. Any field left out takes the built-in default. |

### continuationSettings

| Setting | Default | Meaning |
| --- | --- | --- |
| `MaxContIter` | 1200 | Maximum number of continuation steps (per direction). |
| `MinAlpha` | -50 | Stop when α falls below this value. |
| `MaxAlpha` | 80 | Stop when α rises above this value. |
| `continuationDirection` | `'backward'` | `'forward'`, `'backward'` or `'both'`: which way to trace the branch from the starting point. |
| `InitialIntegrationStepSize` | 0.1 | Initial continuation step dz (arc length along the curve, not a time step). |
| `MinIntegrationStepSize` | 0.02 | Lower bound for dz. |
| `MaxIntegrationStepSize` | 150 | Upper bound for dz. |
| `stepIncreaseFactor` | 1.3 | dz multiplier after an easy corrector solve. |
| `stepDecreaseFactor` | 0.5 | dz multiplier after a hard corrector solve. |
| `maxNewtonIters` | 3 | `fsolve` iteration count from which a solve counts as hard. |
| `runStabilityAnalysis` | true | Compute eigenvalues and fold/Hopf candidates after the run. |
| `livePlotXVar` | `''` (→ α) | Variable name on the live-plot X axis. |
| `livePlotYVar` | `''` (→ 4th state) | Variable name on the live-plot Y axis. |

The example JSON files set their own values; see [Examples](#examples) below.

### Hard-coded in runContinuation.m

- `verbose = false`: set to true to print Newton iterations and dz for each step.
- `livePlot = true`: set to false to switch off the live plot (faster).
- `fsolve` tolerances: TolFun = TolX = 1e-12, MaxFunEvals = 5000.

## Outputs

Each run writes a new file `Res_<problem>_<timestamp>.mat` to the current folder, so earlier runs are kept, and opens the Results Viewer.

| Variable | Size / fields | Content |
| --- | --- | --- |
| `Res` | Nop × (NoEQ+1) | One row per continuation point: states in `stateNames` order, α in the last column. |
| `Bif.stabStatus` | 1 × Nop, logical | true = the point is unstable. |
| `Bif.nUnstable`, `Bif.maxRe` | 1 × Nop | Number of eigenvalues with Re > 0 and the largest real part. |
| `Bif.eigOut` | NoEQ × Nop | All eigenvalues of ∂f/∂x at every point. |
| `Bif.stabChange` | row indices | Points where stability changes. |
| `Bif.fold`, `Bif.hopf` | row indices | The same points, split into fold and Hopf candidates. |
| `info`, `settings` | structs | The problem descriptor and the settings actually used. |

`Bif` is empty when runStabilityAnalysis is false. To reopen a run later, also in a new MATLAB session: `plotContinuationResults('Res_Example3_CSTR_3Reactions_20260927_143015.mat')`, or `plotContinuationResults()` to pick the file.

**Results Viewer:** choose the X and Y columns from dropdowns; click a row in the data grid to highlight that point; the bifurcation table lists fold/Hopf rows with their α. **Dynamic Test** simulates from the selected point: α is kept for the first third of the time, stepped by the given percentage for the middle third, then returned; the response opens in a new window.

## Adding a new problem

1. Write the model function, e.g. `myModel.m` with `F = myModel(X)`, X = [states, parameter] (see [The model function](#the-model-function-rhs)).
2. Copy one of the example JSON files to `myModel.json` and edit `rhs`, `stateNames`, `paramNames`, at least one preset and the α range.
3. Run `runContinuation(loadProblemInfo('myModel.json'))`, check the values in the dialog, press Run.

Tips:

- If the run stops after a few steps, the branch probably left the α range: switch `continuationDirection` (or use `'both'`) or widen MinAlpha / MaxAlpha.
- Repeated "fsolve did not converge" warnings: lower MaxIntegrationStepSize or InitialIntegrationStepSize.
- Close to a fold or Hopf point, compare against a finer run; the classification is a candidate, not a proof.

## Examples

The three examples go from a check with a known answer to a full engineering model. Each is run as `runContinuation(loadProblemInfo('<name>.json'))`.

### Example 1 – Cubic fold (`Example1_CubicFold`)

One equation, one parameter: dx/dτ = α + x − x³. The steady states lie on the S-shaped curve α = x³ − x, with exact turning points at x = ±1/√3 ≈ ±0.5774, α = −2/(3√3) ≈ −0.3849 (for x > 0) and +0.3849 (for x < 0). The outer branches are stable and the middle branch is unstable.

- Start: preset `middle_unstable_start` (x = 0, α = 0), direction `both`, α in [−1, 1], dz 0.01 (0.0005–0.05), at most 300 steps per direction.
- Expected result: the branch turns exactly at the two points above, and both are reported in `Bif.fold`, none in `Bif.hopf` (all eigenvalues are real).
- Use it to test the engine itself after any change to the code. The presets `lower_stable_branch` and `upper_stable_branch` (x = −2, α = −6 and x = 2, α = 6) lie outside the stored α range; widen MinAlpha / MaxAlpha before using them.

### Example 2 – Cascade of two CSTRs with recycle (`Example2_CSTR_cascade`)

Two stirred reactors in series with a first-order exothermic reaction, taken from Kubíček, Jacák and Marek (1983), Section 13.4, Eq. (13.24), p. 313. Four states: conversions x1, x2 and dimensionless temperatures θ1, θ2. The continuation parameter is the Damköhler number Da1 = Da2.

- In the JSON the states are named x1–x4 (x3 = θ1, x4 = θ2) and the parameter is named `alpha` (= Da). The order differs from the book's vector (x1, θ1, x2, θ2).
- Fixed values: γ = 1000, B = 22, β1 = β2 = 2, θc1 = θc2 = 0, Λ = 1.
- Start: all states 0 at Da = 0, which is an exact steady state; direction `forward`, Da in [−1e−5, 0.0595], dz 0.2 (0.1–0.5), factors 1.2 / 0.5, at most 600 steps. The live plot shows Da against θ2.
- Expected result: the temperature–Da curve folds back on itself, so several steady states exist for the same Da. Check the turning points and stability changes in `Bif`.

### Example 3 – CSTR with three reactions and cooling (`Example3_CSTR_3Reactions`)

The model system of Švandová, Jelemenský, Markoš and Molnár (2005): exothermic hydrolysis of propylene oxide to propylene glycol, with consecutive reactions to di- and tripropylene glycol, in a jacket-cooled CSTR. The paper uses this system for steady-state multiplicity and stability analysis within a HAZOP study.

- Seven states: cA (propylene oxide), cB (water), cC (propylene glycol), cD (dipropylene glycol), cE (tripropylene glycol) in mol/m³, reactor temperature T and coolant temperature Tc in K.
- Parameter α shifts the feed temperature: Tf = 273.15 + α K.
- Five presets; the default `feed_temperature` starts at α ≈ −28.6. Direction `backward`, α in [−50, 80], dz 0.1 (0.02–150), at most 1200 steps. The live plot shows α against T.
- Heat transfer: the switch `useTestHeatTransfer` in the model file is false, which uses the physical values U = 1.65e3 W m⁻² K⁻¹ and A = 6.7 m².
- `runDynamicTest.m` simulates this model for step changes of α starting from a fixed steady state.

## References

- M. Kubíček: [Algorithm 502: Dependence of Solution of Nonlinear Systems on a Parameter](https://dl.acm.org/doi/abs/10.1145/355666.355675). ACM Transactions on Mathematical Software 2(1), 1976.
- M. Kubíček, V. Jacák, M. Marek: *Numerické algoritmy řešení chemicko-inženýrských úloh*. SNTL / Alfa, Praha, 1983 (SNTL 04-614-83).
- M. Kubíček, M. Marek: *Computational Methods in Bifurcation Theory and Dissipative Structures*. Springer, 1983.
- Z. Švandová, Ľ. Jelemenský, J. Markoš, A. Molnár: [Steady States Analysis and Dynamic Simulation as a Complement in the HAZOP Study of Chemical Reactors](https://doi.org/10.1205/psep.04262). Process Safety and Environmental Protection 83(5), 2005, pp. 463–471.
