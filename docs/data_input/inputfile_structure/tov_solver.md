# TOVSolver

## Overview

This entry contains all relevant information regarding Tolman-Oppenheimer-Volkov (TOV, {cite}`PhysRev.55.364`) system of PDEs. In particular, one specifies here radius discretization, central energy density for IVP, $P(\rho)$ EoS function discretization etc.
In what follows we discuss all the possible settings under this entry.

This section is essential for any simulation and is therefore mandatory to include.

In order to refer to an example, see <span style="color:blue">_presupplied/APR4/RHMconfig.json_</span> and <span style="color:blue">_presupplied/APR4/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"AdaptionLimit"` (uint, [<span style="color:red">TOV, COOL, RH</span>]) **:** Defines how many times TOV solver will try to adapt the radial step before concluding early termination. Defaults to 20.
```{note}
Early termination can be overwatched in "Debug" log.
```
- `"EoSInterpolation"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** Interpolation kind to be used for $\rho(P)$ dependence. Choose from ["Linear", "Cubic"], with "Cubic" being default. 
- `"DensityInterpolation"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** Interpolation kind to be used for $n_{bar}(P)$ dependence. Choose from ["Linear", "Cubic"], with "Linear" being default.
```{warning}
"Cubic" mode strongly distorts the $n_{bar}-M$ relation at phase transitions, even weak ones. Select on your own risk.
```
- `"RadialInterpolation"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** Interpolation kind to be used for discretized $m(r), P(r), \rho(r), ...$ dependences. Choose from ["Linear", "Cubic"], with "Cubic" being default. 
- `"LengthUnits"` (string/double, [<span style="color:red">TOV, COOL, RH</span>]) **:** Conversion factor from length to natural units (GeV powers). It must either be supplied as a choice from ["Gev-1", "Km", "M", "Cm"], or as an actual multiplier. Used for "RadiusStep"

- `"RadiusStep"` (double, [<span style="color:red">TOV, COOL, RH</span>]) **:** Defines radius discretization step. Units are defined by "LengthUnits" entry.

- `"CenterPressure"` **:** Defines central pressure for initial value problem.
    - `"ProvidedAs"` (string, [<span style="color:red">COOL, RH</span>]) **:** The way the "CenterPressure" is provided. Choose from ["Same", "LinspacedMinToMax", "MassCached"]. "LinspacedMinToMax" linearly maps $[0,1] \rightarrow [P_{min}, P_{max}]$, while "Same" assumes same units as for "Pressure". "MassCached" instead is fed with a star's desired mass (solar mass units) and the corresponding pressure is calculated from cached TOV output.
    - `"Value"` (double, [<span style="color:red">COOL, RH</span>]) **:** Value of central pressure. Interpretation depends on "ProvidedAs" entry.
    - `"CachePath"` (string, [<span style="color:red">COOL, RH</span>]) **:** Path to the cached TOV output. Only used if "ProvidedAs" is set to "MassCached". In order to produce cache, run
    ```bash
    ./bin/project/M-R_diagram/m_r_diagram.out --inputfile <path_to_inputfile> [OPTIONS] > <path_to_cache>
    ```
    , where options are inspectable via `--help` flag.
    ```{warning}
    1) One is not enforced to have TOV output ordered by mass, since this may be complicated to achieve. Because of this, the code looks for the first interval containing the desired mass. To limit the possible range, use `--left_fraction` and `--right_fraction` flags during cache creation.
    2) The resulting mass is not guaranteed to coincide with the desired mass, since the value is linearly interpolated. Pass a bigger `--selection_size` to increase precision.
    3) Avoid mixing different EoSs with the same cache. This can either lead to out-of-range errors or to wrong mass values.
    ```

```{note}
Though this setting is a key value for TOV solver, it is not used for producing M-R curves, since M-R curves are _parametrized_ by pressure. It, of course, does not undermine its importance during any simulation on a given pressure.
```
- `"SurfacePressure"` **:** Defines surface pressure as an ultimate exit condition for TOV.
    - `"ProvidedAs"` (string, [<span style="color:red">TOV, COOL, RH</span>]) **:** The way the "SurfacePressure" is provided. Choose from ["Same", "LinspacedMinToMax"]. "LinspacedMinToMax" linearly maps $[0,1] \rightarrow [P_{min}, P_{max}]$, while "Same" assumes same units as for "Pressure".
    - `"Value"` (double, [<span style="color:red">TOV, COOL, RH</span>]) **:** Value of surface pressure. Interpretation depends on "ProvidedAs" entry.