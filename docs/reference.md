# Pharmacokinetics reference

Library version **0.1** · [Package overview](../README.md)

Use one consistent system of amount, volume and time units. States in ODE
functions are **amounts**. Closed-form functions and `nonlinear_concentration`
return **central concentrations**, amount divided by central volume. Clearance
and intercompartmental clearance have units volume/time. Absorption has units
1/time. Bioavailability is dimensionless.

Closed-form functions take elapsed time since administration and return zero
before administration. Boluses include the dose at time zero. Pass a vector
of times directly: `pk.one_compartment_oral(times, dose, p)`. Each function
broadcasts its numeric inputs and keeps the parameter record fixed.

| Function | Purpose |
|---|---|
| [`one_compartment_bolus`](compartments.md#one-compartment-bolus) | One-compartment bolus concentration |
| [`one_compartment_oral`](compartments.md#one-compartment-oral) | First-order oral absorption |
| [`one_compartment_infusion`](compartments.md#one-compartment-infusion) | Finite constant-rate infusion |
| [`two_compartment_bolus`](compartments.md#two-compartment-bolus) | Central concentration after a bolus |
| [`two_compartment_oral`](compartments.md#two-compartment-oral) | Two-compartment oral absorption |
| [`two_compartment_infusion`](compartments.md#two-compartment-infusion) | Two-compartment finite infusion |
| [`dose_superposition`](#dose_superposition) | Sum repeated linear bolus or oral responses |
| [`infusion_superposition`](#infusion_superposition) | Sum repeated linear infusions |
| [`one_compartment_rhs`](#one_compartment_rhs) | Linear central-amount dynamics |
| [`two_compartment_rhs`](#two_compartment_rhs) | Linear central/peripheral-amount dynamics |
| [`michaelis_menten_rhs`](#michaelis_menten_rhs) | Saturable central elimination |
| [`nonlinear_concentration`](#nonlinear_concentration) | Solve saturable elimination on a chosen grid |
| [`additive_error`](#additive_error) | Gaussian concentration errors |
| [`combined_error`](#combined_error) | Additive and proportional Gaussian errors |
| [`lognormal_error`](#lognormal_error) | Positive multiplicative errors |

## `dose_superposition`

`dose_superposition(response, time, dose_times, doses, p)` returns

$$C(t)=\sum_j R(t-t_j,d_j,p).$$

Pass a bolus or oral response function from the table. `dose_times` and `doses`
are aligned nonempty vectors. Their order does not matter. The responses return
zero before each dose. The model must be linear: do not use superposition for
Michaelis–Menten elimination.

## `infusion_superposition`

`infusion_superposition(response, time, starts, rates, durations, p)` sums
`response(time - starts[j], rates[j], durations[j], p)`. Pass either finite
infusion response. The three vectors are aligned and nonempty. Infusions may
overlap; rates are amounts/time and durations are positive times.

## `one_compartment_rhs`

`one_compartment_rhs(t, y, p)` returns the derivative of `[central_amount]`:

$$\dot A=R-\frac{CL}{V}A.$$

The fields of `p` are `input_rate`, `clearance`, `volume`. The input rate is
constant and nonnegative. Clearance and volume are positive. Time `t` is accepted
for the ODE interface. With a changing infusion rate, solve separate intervals
with the appropriate rate and carry the final state forward.

## `two_compartment_rhs`

`two_compartment_rhs(t, y, p)` takes `[central_amount, peripheral_amount]`:

$$\dot A_c=R-CL\frac{A_c}{V_c}-Q\left(\frac{A_c}{V_c}-\frac{A_p}{V_p}\right),$$
$$\dot A_p=Q\left(\frac{A_c}{V_c}-\frac{A_p}{V_p}\right).$$

Fields are `input_rate`, `clearance`, `central_volume`, `peripheral_volume`,
`intercompartmental_clearance`. The latter is $Q$. All except the nonnegative
input rate are positive. Exchange conserves total amount; elimination and input
change it. Neither compartment has an absorption depot in this RHS.

## `michaelis_menten_rhs`

`michaelis_menten_rhs(t, y, p)` takes `[central_amount]` and returns

$$\dot A=R-\frac{V_{\max} A}{K_m V+A}.$$

Fields are `input_rate`, `vmax`, `km`, `volume`. Input rate is nonnegative;
the other fields are positive. `vmax` has units amount/time, `km` has units
concentration. The elimination rate reaches half of `vmax` at concentration
`km`. This model has no additional linear clearance term.

## `nonlinear_concentration`

`nonlinear_concentration(solver, initial_amount, times, p)` calls the supplied
ODE solver with `michaelis_menten_rhs`. It returns a concentration vector at
`times`, including `initial_amount / p.volume` at `times[1]`. Use nonnegative
initial amount and a strictly increasing grid with at least two times.

For example, supply `ode.vern7` or `ode.rosenbrock23` from
[`flatppl-numerics-ode`](https://github.com/flatppl-community/flatppl-numerics-ode).
The caller chooses the solver and grid. There are no automatic dose events,
positivity correction, adaptive steps or inferred priors. Refine the grid and
check posterior stability. Intermediate ODE stages must stay in the RHS domain;
a coarse grid may produce unphysical amounts. See the ODE reference for the
stiff methods' finite-difference and Newton limits.

## `additive_error`

`additive_error(concentration, sigma)` returns independent
$Y_i\sim\mathrm{Normal}(C_i,\sigma)$ observation laws. `concentration` is a
vector and `sigma` is positive. Gaussian errors permit negative measurements.

## `combined_error`

`combined_error(concentration, sigma_add, sigma_prop)` returns independent
Gaussian laws with means $C_i$ and standard deviations

$$s_i=\sqrt{\sigma_{\rm add}^2+(\sigma_{\rm prop} C_i)^2}.$$

Both noise parameters are nonnegative, and every $s_i$ must be positive.
`sigma_add` has concentration units; `sigma_prop` is dimensionless. This adds
variances, not standard deviations. Both noise parameters are scalars.

## `lognormal_error`

`lognormal_error(concentration, sigma_log)` returns independent
$\log Y_i\sim\mathrm{Normal}(\log C_i,\sigma_{\log})$ laws. Every $C_i$ and
`sigma_log` must be positive. Here $C_i$ is the **median**, not the mean:
$E[Y_i]=C_i\exp(\sigma_{\log}^2/2)$. Zero observations require a different
observation model or an explicit censoring model.

## Examples and validation

The [oral-dose example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-pharma-pk%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)
infers clearance and volume from six synthetic measurements. It imports the
library from a pinned GitHub commit. The [nonlinear example](../examples/nonlinear.flatppl)
infers a saturable elimination rate with
[flatppl-numerics-ode](https://github.com/flatppl-community/flatppl-numerics-ode).

The library provides one- and two-compartment bolus, oral and finite-infusion
solutions, repeated linear doses, amount-based ODE functions and concentration
observation laws. Parameters and priors are explicit. Use consistent units:
ODE states are amounts; prediction functions return central concentrations.
See the [reference](reference.md) for domains, formulas and solver limits.

Run the mathematical checks with Node 24 or newer and a JS engine checkout:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js node tests/check.cjs
```
