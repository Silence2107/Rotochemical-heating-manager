# TOVSolver

## Overview

This entry contains all relevant information regarding Tolman-Oppenheimer-Volkov (TOV, {cite}`PhysRev.55.364`) system of PDEs. In particular, one specifies here radius discretization, central energy density for IVP, $P(\rho)$ EoS function discretization etc.
In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"EoSDiscretization"` (uint) **:** Defines discretization step for $P(\rho)$ dependence. Defaults to 1000.
- `"EoSInterpolation` (string) **:** Interpolation kind to be used for discretized $P(\rho)$ dependence. Choose from ["Linear", "Cubic"], with "Linear" being default. 
- `"BarionicDensityInterpolation"` (string) **:** Interpolation kind to be used for discretized $n_b(r)$ dependence. Choose from ["Linear", "Cubic"], with "Linear" being default. 
```{note}
This setting is not used by TOV solver itself (rather for cooling functionality), but since it originates from TOV solver, it is placed here. Something to consider moving.
```
- `"LengthUnits"` (string/double, required) **:** Conversion factor from length to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "Km", "M", "Cm"], or as an actual multiplier. Used for "RadiusStep"

- `"RadiusStep"` (double, required) **:** Defines radius discretization step. Units are defined by "LengthUnits" entry.

- `"DensityUnits"` (double, required) **:** Conversion factor from energy density to natural units (GeV powers). It must either be supplied as a choice from ["Gev4", "Same", "RelativeToMax"], or as an actual multiplier. Used for "CenterDensity" and "DensityStep". "RelativeToMax" scales against `["EnergyDensity"]["Upp"]` entry, while "Same" assumes same units as for "EnergyDensity".

- `"CenterDensity"` (double, required) **:** Defines central energy density for initial value problem. Units are defined by "DensityUnits" entry.
```{note}
Though this setting is a key value for TOV solver, it is not used for producing M-R curves, since M-R curves are parametrized by central density. It, of course, does not undermine its importance during any simulation on a given central density.
```
- `"DensityStep"` (double, required) **:** Defines energy density at which star's radius is calculated. This same quantity is used for $P(\rho)$ differentiation, as well as it enforces TOV recaching if center density is changed by that much during one program's execution. Units are defined by "DensityUnits" entry.
```{note}
Appears to be very cumbersome variable. I should consider relaxing its responsibility.
```