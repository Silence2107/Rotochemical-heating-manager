# TOVSolver

## Overview

This entry contains all relevant information regarding Tolman-Oppenheimer-Volkov (TOV, {cite}`PhysRev.55.364`) system of PDEs. In particular, one specifies here radius discretization, central energy density for IVP, $P(\rho)$ EoS function discretization etc.
In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"EoSDiscretization"` (uint) **:** Defines discretization step for $P(\rho)$ dependence. Defaults to 1000.
- `"AdaptionLimit"` (uint) **:** Defines how many times TOV solver will try to adapt the radial step before concluding early termination. Defaults to 20.
```{warning}
Early termination happens silently and is only well noticeable if the final star's radius is very low.
```
- `"EoSInterpolation` (string) **:** Interpolation kind to be used for discretized $P(\rho)$ dependence. Choose from ["Linear", "Cubic"], with "Linear" being default. 
- `"BarionicDensityInterpolation"` (string) **:** Interpolation kind to be used for discretized $n_b(r)$ dependence. Choose from ["Linear", "Cubic"], with "Linear" being default. 
```{note}
This setting is not used by TOV solver itself (rather for cooling functionality), but since it originates from TOV solver, it is placed here. Something to consider moving.
```
- `"LengthUnits"` (string/double, required) **:** Conversion factor from length to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "Km", "M", "Cm"], or as an actual multiplier. Used for "RadiusStep"

- `"RadiusStep"` (double, required) **:** Defines radius discretization step. Units are defined by "LengthUnits" entry.

- `"CenterDensity"` **:** Defines central energy density for initial value problem.
    - `"ProvidedAs"` (string, required) **:** The way the "CenterDensity" is provided. Choose from ["Same", "LinspacedMinToMax"]. "LinspacedMinToMax" linearly maps $[0,1] \rightarrow [\rho_{min}, \rho_{max}]$, while "Same" assumes same units as for "EnergyDensity".
    - `"Value"` (double, required) **:** Value of central energy density. Interpretation depends on "ProvidedAs" entry.
```{note}
Though this setting is a key value for TOV solver, it is not used for producing M-R curves, since M-R curves are _parametrized_ by central density. It, of course, does not undermine its importance during any simulation on a given central density.
```
- `"SurfacePressure"` **:** Defines surface pressure as an ultimate exit condition for TOV.
    - `"ProvidedAs"` (string, required) **:** The way the "SurfacePressure" is provided. Choose from ["Same", "LinspacedMinToMax"]. "LinspacedMinToMax" linearly maps $[0,1] \rightarrow [P_{min}, P_{max}]$, while "Same" assumes same units as for "Pressure".
    - `"Value"` (double, required) **:** Value of surface pressure. Interpretation depends on "ProvidedAs" entry.