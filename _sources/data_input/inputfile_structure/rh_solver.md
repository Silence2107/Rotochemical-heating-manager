# RHSolver

## Overview

This entry contains all relevant settings for running cooling simulation with rotochemical heating {cite}`Fernandez:2005cg`.

## Description

- `"TimeUnits"` (double, [<span style="color:red">RH</span>]) **:** Conversion factor from time to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "S", "Yr", "Ms"], or as an actual multiplier. Used for "RotationalOmegaSquareDot" parameters.
- `"RotationalOmegaSquareDot"` (array, [<span style="color:red">RH</span>]) **:** Rotational frequency settings. The first element is a string with mode, while the others denote parameters for this mode. 
```{note}
The only supported mode is "BeyondMagneticDipole", which is described in {cite}`yanagi2020thermal` (Eq. 2.102) and is given by

$$
\Omega(t) = \dfrac{2\pi}{P_0}\left(1 + \dfrac{(n-1)\dot P_0}{P_0}t\right)^{-\frac{1}{n-1}}.
$$

We denote $[n, P_0, \dot P_0]$ via array[1-3]. $P_0$ is measured in "TimeUnits".
```