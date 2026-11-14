# flatppl-pharma-pk

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-pharma-pk%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

One- and two-compartment pharmacokinetics for FlatPPL, including bolus, oral, finite-infusion, and repeated-dose models. Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
pk = load_module("main.flatppl")
p = record(clearance = 5.0, volume = 20.0, absorption = 1.2, bioavailability = 0.8)
concentrations = pk.one_compartment_oral([0.5, 1.0, 2.0, 4.0], 100.0, p)
```

ODE states represent amounts. Prediction functions return central concentrations. Parameters, observation laws, and priors remain explicit. Units must agree.

The [oral-dose example](examples/playground.flatppl) infers clearance and volume from six synthetic measurements. The [nonlinear example](examples/nonlinear.flatppl) infers a saturable elimination rate.

See the reference for definitions, domains, numerical limits, and validation commands.
