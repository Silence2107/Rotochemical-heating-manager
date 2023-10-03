# CoolingSolver

## Overview

This entry contains all relevant information regarding star's cooling PDE {cite}`Gudmundsson_crust`. In particular, one specifies here radius discretization, Newton step tolerances, simulation time limits, initial temperature profile etc. In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"NewtonTolerance"` (double) **:** Defines Newton step tolerance (desired relative difference between subsequent iterations) for cooling PDE. Defaults to $10^{-5}$
- `"NewtonMaxIter"` (uint) **:** Defines maximum number of Newton iterations for cooling PDE. Defaults to 50 
- `"TimeInit"` (double) **:** Defines initial time for cooling PDE in years. Defaults to 0 years.
- `"TimeEnd"` (double, required*) **:** Defines final time for cooling PDE in years.
- `"TimeBaseStep"` (double, required*) **:** Defines base time step for cooling PDE in years.
- `NumberPointsEstimate` (uint, required*) **:** Defines an estimation (disregarding algorithms adaptivity and imprecision) for number of points for cooling PDE. 
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