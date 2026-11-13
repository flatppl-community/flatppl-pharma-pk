# flatppl-pharma-pk

Compartment pharmacokinetics in FlatPPL. Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
pk = load_module("main.flatppl")
p = record(clearance = 5.0, volume = 20.0, absorption = 1.2, bioavailability = 0.8)
concentrations = pk.one_compartment_oral([0.5, 1.0, 2.0, 4.0], 100.0, p)
```

The [oral-dose example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-pharma-pk%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)
infers clearance and volume from six synthetic measurements. It imports the
library from a pinned GitHub commit. The [nonlinear example](examples/nonlinear.flatppl)
infers a saturable elimination rate with
[flatppl-numerics-ode](https://github.com/flatppl-community/flatppl-numerics-ode).

The library provides one- and two-compartment bolus, oral and finite-infusion
solutions, repeated linear doses, amount-based ODE functions and concentration
observation laws. Parameters and priors are explicit. Use consistent units:
ODE states are amounts; prediction functions return central concentrations.
See the [reference](docs/reference.md) for domains, formulas and solver limits.

Run the mathematical checks with Node 24 or newer and a JS engine checkout:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js node tests/check.cjs
```
