# Closed-form compartment models

Every function returns **central concentration**, in amount/volume, after one
dose into an initially empty system. `time` is elapsed time from the dose or
infusion start. All functions return zero for negative time. At zero, a bolus
returns `dose / volume` (or `dose / central_volume`); oral dosing and infusion
return zero. Infusion concentration is continuous when the infusion ends.

Use consistent units: `time` and `duration` have time units; `dose` has amount
units; `rate` has amount/time units; clearance and intercompartmental clearance
have volume/time units; volumes have volume units; absorption has 1/time units.
Bioavailability is dimensionless. Inputs must be finite, volumes, clearances,
absorption, and duration strictly positive, dose and rate nonnegative, and
`0 < bioavailability <= 1`. These are caller preconditions; the functions do
not validate them. No endogenous baseline or previous dose is implicit.

Call these functions normally with scalar or vector numeric inputs; their
numeric arguments broadcast together. For example,
`pk.one_compartment_oral(times, 100, parameters)` returns one concentration
per time. `p` is one record of scalar parameters, held across those evaluations.
Do not dot-call the public function with `p` as a broadcast operand or wrap `p`
in an array: records are not broadcast operands or array elements in FlatPPL.

| Function | Parameter record `p` | Definition |
| --- | --- | --- |
| `one_compartment_bolus(time, dose, p)` | `clearance`, `volume` | [One-compartment bolus](#one-compartment-bolus) |
| `one_compartment_oral(time, dose, p)` | Above, `absorption`, `bioavailability` | [One-compartment oral](#one-compartment-oral) |
| `one_compartment_infusion(time, rate, duration, p)` | `clearance`, `volume` | [One-compartment infusion](#one-compartment-infusion) |
| `two_compartment_bolus(time, dose, p)` | `clearance`, `central_volume`, `peripheral_volume`, `intercompartmental_clearance` | [Two-compartment bolus](#two-compartment-bolus) |
| `two_compartment_oral(time, dose, p)` | Above, `absorption`, `bioavailability` | [Two-compartment oral](#two-compartment-oral) |
| `two_compartment_infusion(time, rate, duration, p)` | `clearance`, `central_volume`, `peripheral_volume`, `intercompartmental_clearance` | [Two-compartment infusion](#two-compartment-infusion) |

The definitions below use $t\geq0$, dose $D$, infusion rate $R$, duration $T$,
absorption rate $k_a$, and bioavailability $F$. Write

$$
E(a,b,t)=\int_0^t e^{-a(t-s)}e^{-bs}\,ds
=\begin{cases}
\dfrac{e^{-at}-e^{-bt}}{b-a},&a\ne b,\\
t e^{-at},&a=b,
\end{cases}
$$

and, with $u=\min(t,T)$ and $v=\max(t-T,0)$,

$$
I(k,t,T)=e^{-kv}\frac{1-e^{-ku}}{k}.
$$

## One-compartment bolus

For $k=CL/V$, central amount obeys $A'=-kA$, $A(0)=D$. Thus
$C(t)=D e^{-kt}/V$.

## One-compartment oral

Depot amount obeys $A_d'=-k_a A_d$, $A_d(0)=FD$, and central amount obeys
$A'=k_a A_d-kA$, $A(0)=0$. Thus $C(t)=FDk_a E(k,k_a,t)/V$.
At $k_a=k$ this is $FDkt e^{-kt}/V$; differentiation treats $k_a$ and $k$
as independent parameters before evaluating that limit.

## One-compartment infusion

Central input is $R$ from zero to $T$, then zero. The solution is
$C(t)=R I(k,t,T)/V$. In particular, $C(0)=0$ and
$C(T)=R(1-e^{-kT})/(Vk)$.

## Two-compartment bolus

Let $V_c$ and $V_p$ be the central and peripheral volumes and $Q$ the
intercompartmental clearance. Define

$$
k_{10}=CL/V_c,\quad k_{12}=Q/V_c,\quad k_{21}=Q/V_p,
$$

$$
\Delta=\sqrt{(k_{10}+k_{12}-k_{21})^2+4k_{12}k_{21}},\quad
\alpha=\frac{k_{10}+k_{12}+k_{21}+\Delta}{2},\quad
\beta=\frac{k_{10}k_{21}}{\alpha},
$$

$$
w_\alpha=\frac{k_{10}+k_{12}-\beta}{\alpha-\beta},\qquad
w_\beta=\frac{k_{21}-\beta}{\alpha-\beta}.
$$

With central and peripheral amounts $A_c,A_p$, the equations are

$$
A_c'=-(k_{10}+k_{12})A_c+k_{21}A_p,\qquad
A_p'=k_{12}A_c-k_{21}A_p.
$$

Starting at $(A_c,A_p)=(D,0)$ gives
$C(t)=D(w_\alpha e^{-\alpha t}+w_\beta e^{-\beta t})/V_c$.
The weights sum to one and are evaluated without subtracting nearly equal
positive rates. Strictly positive $Q$ makes $\Delta>0$.

## Two-compartment oral

Add the depot $A_d'=-k_aA_d$, $A_d(0)=FD$, with input $k_aA_d$ to the
central equation and initially empty central and peripheral compartments:

$$
C(t)=\frac{FDk_a}{V_c}\left[
w_\alpha E(\alpha,k_a,t)+w_\beta E(\beta,k_a,t)\right].
$$

The same $E$ definition covers $k_a=\alpha$ and $k_a=\beta$.

## Two-compartment infusion

Add input $R$ to the central equation from zero to $T$, with initially empty
compartments:

$$
C(t)=\frac{R}{V_c}\left[
w_\alpha I(\alpha,t,T)+w_\beta I(\beta,t,T)\right].
$$

## Numerical evaluation and references

The implementation uses `expm1` and a sixth-degree local series for
$(1-e^{-x})/x$ near zero. Near coincident rates, $E$ uses the symmetric form
$t e^{-(a+b)t/2}\sinh((a-b)t/2)/((a-b)t/2)$ with an even sixth-degree series.
This preserves the distinct parameter derivatives at equality. Polynomial
arguments and inactive denominators are bounded because FlatPPL `ifelse`
does not require short-circuit evaluation. Infusion uses elapsed input time
and decay since stopping, avoiding subtraction of two step responses.

These are finite-precision formulas, not arbitrary-precision solvers. Large
rate/time products may underflow to zero; extreme scales can overflow
intermediate products, and nearly coincident two-compartment eigenvalues can
lose derivative precision. The $Q=0$ boundary and exact repeated eigenvalues
are outside the domain. Use the one-compartment functions when there is no
peripheral exchange.

Torsten documents the corresponding compartment conventions and central
amount-to-concentration conversion in its
[one-compartment](https://metrumresearchgroup.github.io/Torsten/function/one-cpt/)
and [two-compartment](https://metrumresearchgroup.github.io/Torsten/function/two-cpt/)
model references. The functions here solve one dose or finite infusion;
they do not implement Torsten's event scheduler or steady-state protocol.
