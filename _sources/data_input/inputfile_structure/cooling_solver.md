# CoolingSolver

## Overview

This entry contains all relevant information regarding star's cooling PDE {cite}`Gudmundsson_crust`. In particular, one specifies here radius discretization, Newton step tolerances, simulation time limits, initial temperature profile etc. In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"NewtonTolerance"` (double) **:** Defines Newton step tolerance (desired relative difference between subsequent iterations) for cooling PDE. Defaults to $10^{-5}$
- `"NewtonMaxIter"` (uint) **:** Defines maximum number of Newton iterations for cooling PDE. Defaults to 50 
- `"TimeUnits"` (string/double, required*) **:** Conversion factor from time to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "S", "Yr", "Myr"], or as an actual multiplier. Used for "TimeBaseStep", "TimeInit", "TimeEnd" and "UponReachingTime".
- `"TimeInit"` (double) **:** Defines initial time for cooling PDE in "TimeUnits". Defaults to 0 "TimeUnits".
- `"TimeEnd"` (double, required*) **:** Defines final time for cooling PDE in "TimeUnits".
- `"TimeBaseStep"` (double, required*) **:** Defines base time step for cooling PDE in "TimeUnits".
- `"NumberPointsEstimate"` (uint, required*) **:** Defines an estimation (disregarding algorithms adaptivity and imprecision) for number of points for cooling PDE. 
```{note}
This setting is used for reserving memory for the solution, as well as it could be used for estimating the time grid.
```
- `"ExpansionRate"` (double/string) : Defines an inflation parameter for the time step, i.e. 

$$
\Delta t_{n} = \Delta t_{\text{base}} \cdot a^{n},
$$ 

where a is expansion rate. This is a very efficient solution to cover multiple timescales in rather finite number of time steps. One either may supply their own expansion rate or provide it as "Deduce" string, which is also the default behaviour.

```{note}
Deducing expansion rate is done via approximately limiting the summarized time grid:

$$
\sum_{n=0}^{N} \Delta t_{n} = \sum_{n=0}^{N} \Delta t_{\text{base}} \cdot a^{n} = \Delta t_{\text{base}} \cdot \frac{a^{N+1}-1}{a-1} \approx T_{\text{end}} - T_{\text{init}},
$$

where $T_{\text{end}}$ and $T_{\text{init}}$ are final and initial times respectively, and $N$ is number of points estimate. This equation is solved under assumptions $T_{\text{end}}$ - $T_{\text{init}} \gg \Delta t_{\text{base}}$ and $N \gg 1$, which in second order yields:

$$
a_0 = \left(\frac{T_{\text{end}} - T_{\text{init}}}{\Delta t_{\text{base}}} \right)^{1/N} \Rightarrow \\
a = a_0 \cdot (a_0-1)^{1/N}.
$$

Actual precision of this estimate is unpredictable (given that our PDE solvers are time-adaptive, which is not taken into account), but it is expected to suffice in most applications.
```
- `"TemperatureUnits"` (string/double, required) **:** Conversion factor from temperature to natural units (GeV powers). It must either be supplied as a choice from ["Gev", "MeV", "K"], or as an actual multiplier. Used for "TemperatureProfile" arguments. 
- `"TemperatureProfile"` (array[string, args..], required*) **:** Initial temperature profile (radial dependence) settings.

If array[0] is supplied as "InfiniteFlat", then the initial temperature profile $T^{\infty}$ is set to be constant and equal to the value supplied in array[1], i.e.

$$
T^{\infty}(r, t=0) = \text{array[1]} \cdot \text{conversion}.
$$

If array[0] is supplied as "InfiniteFlatSurfaceRedshifted", then the profile is set via $\Phi$ metric function on the surface and constant value supplied in array[1], i.e.:

$$
T^{\infty}(r, t=0) = \text{array[1]} \cdot \exp\left(\Phi(R)\right) \cdot \text{conversion}.
$$

```{admonition} devnote
It should be made possible to define "LocalFlat" profile.
```

- `"LengthUnits"` (string/double, required) **:** Conversion factor from length to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "Km", "M", "Cm"], or as an actual multiplier. Used for "RadiusStep"
- `"RadiusStep"` (double, required*) **:** Defines radius step for cooling PDE. Units are defined by "LengthUnits" entry.
- `"EnableEquilibrium"` **:** Settings reflecting code's ability to switch to equilibrium solver. 
    ```{note}
    One may want to switch for few reasons:
    - To substantially speed up the simulation, as equilibrium solver is much faster than generic nonequilibrium one.
    - At late times (usually $t \ge 10^{3-4}$ yrs) the star is expected to be in equilibrium, so it is reasonable.
    - To avoid numerical instabilities in nonequilibrium solver, which may occur (technically, at any time).
    - To enable rotochemical heating simulation, which is only (as of now!) formulated in terms of equilibrium solver quantities.
    ```
    - `"Mode"` (string) **:** Choose the way you'd like equilibrium solver to be enabled. Choose from ["Immediately", "Never", "Conditional"], with "Never" being a default. As the names suggest, "Never" mode never enables equilibrium solver, "Immediately" mode enables it immediately after the start of the simulation, and "Conditional" mode relies on external conditions (see below).
    - `"Conditions"` **:** Choice of conditions upon which cooling mode will get switched.
        - `"UponReachingTime"` (double) **:** Provide time in "TimeUnits", upon reaching which the equilibrium solver will be enabled. Disabled by default.
        - `"UponProfileFlattening"` (double) **:** Provide desired flattening ratio $\left|\frac{T^{\infty}(R) - T^{\infty}(0)}{T^{\infty}(R)}\right|$ (dimensionless), upon reaching which the equilibrium solver will be enabled. Disabled by default.
        ```{note}
        If none conditions are supplemented, then the equilibrium solver will get invoked immediately.
        ```

