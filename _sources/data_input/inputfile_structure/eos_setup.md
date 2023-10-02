# EoSSetup

## Overview

This entry contains all relevant information regarding how program is supposed to interpret the supplied EoS. In particular, one specifies here EoS datafile format, interpolation kind, column to quantity mapping, units, particle species etc. In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

<!-- \begin{itemize}
        \item "Datafile" : Contains all relevant knowledge regarding EoS datafile structure.
        \begin{itemize}
            \item "Path" (string, required) : valid path to the EoS datafile.
            \item "Rows" (uint pair) : [first, last) row in the inputfile to consider (counting from 0th). Assumes whole file, if not specified.
            \item "Columns" (uint pair) : [first, last) row in the inputfile to consider (counting from 0th). Deduces from the first row, if not specified.
            \item "Interpolation" (string) : Interpolation kind to be used across all columns. Choose from ["Linear", "Cubic"], with "Linear" being default. Using "Cubic" mode may lead to issues with positivity of some functions, so it is at the moment advised to instead populate the datafile with more points, yet to use "Linear" mode.
        \end{itemize}
        \item "Particles" (string array, required*) : Names for the particles, available in the EoS. Must strictly match with predefined list ["Neutron", "Proton", "Electron", "Muon", "Tau", "Uquark", "Dquark", "Squark"]. If provided, the instantiator assumes user would like to invoke cooling functionality, hence it requires more provided variables to run.
        \item "Quantities" : A linker between EoS columns and the way the supply data to RHM.
        \begin{itemize}
            \item "BaryonicDensity" : total baryonic density at the point.
            \begin{itemize}
                \item "Units" (string/double, required) : Conversion factor to natural units (GeV powers). It must either be supplied as a choice from ["Gev3", "Fm-3"], or as an actual multiplier.
                Example : if your units are $m^{-3}$, the factor must be set to $\frac{1}{\text{m}^3\cdot \text{GeV}^3} \approx 7.72\text{E-}48$, since $n_b[\text{GeV}^3] = \dfrac{n_b[\text{m}^{-3}] \cdot \text{m}^{-3}}{\text{GeV}^3} = n_b[\text{m}^{-3}] \cdot \frac{1}{\text{m}^3\cdot \text{GeV}^3}$
                \item "Column" (uint, required) : Column number with baryonic density, counting from 0th.
                \item "Low" (double, required) : Smallest accessible $n_b$ within the datafile.
                \item "Upp" (double, required) : Biggest accessible $n_b$ within the datafile.
                \item "CoreLimit" (double, required) : First $n_b$ upon entering core.
                \item "CrustLimit" (double) : Last $n_b$ upon leaving crust. If not supplied, "CoreLimit" value is substituted.
            \end{itemize}
            \item "EnergyDensity" : Total energy density at the point.
            \begin{itemize}
                \item "Units" (string/double, required) : Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevOverFm3", "GOverCm3"], or specify an actual multiplier.
                \item "Column" (uint, required) : Column number with energy density, counting from 0th.
                \item "Low" (double/string) : Smallest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
                \item "Upp" (double/string) : Biggest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
                \item "CoreLimit" (double/string) : First $\rho_E$ upon entering core. Deduced automatically if specified as "Deduce" or if left blank.
            \end{itemize}
            \item "Pressure" : Total pressure at the point.
            \begin{itemize}
                \item "Units" (string/double, required) : Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevOverFm3", "DyneOverCm2"], or specify an actual multiplier.
                \item "Column" (uint, required) : Column number with pressure, counting from 0th.
                \item "Low" (double/string) : Smallest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
                \item "Upp" (double/string) : Biggest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
            \end{itemize}
            \item "BarionicDensities" : Number densities per each particle at the point.
            \begin{itemize}
                \item "PARTICLE\_NAME" : Properties, related to a specific particle. Substitute the entry name with actual (fermionic only) species among provided in "Particles" array one by one.
                \item "Units" (string/double) : Conversion factor to natural units (GeV powers). Choose from ["Gev3", "Fm-3", "DimLess", "Gev", "Fm-1"], or specify an actual multiplier.
            \end{itemize}
        \end{itemize}
    \end{itemize} -->

- `"Datafile"` **:** Contains all relevant knowledge regarding EoS datafile structure.
    - `"Path"` (string, required) **:** valid path to the EoS datafile.
    - `"Rows"` (uint pair) **:** [first, last) row in the inputfile to consider (counting from 0th). Assumes whole file, if not specified.
    - `"Columns"` (uint pair) **:** [first, last) row in the inputfile to consider (counting from 0th). Deduces from the first row, if not specified.
    - `"Interpolation"` (string) **:** Interpolation kind to be used across all columns. Choose from ["Linear", "Cubic"], with "Linear" being default. Using "Cubic" mode may lead to issues with positivity of some functions, so it is at the moment advised to instead populate the datafile with more points, yet to use "Linear" mode.
- `"Particles"` (string array, required*) **:** Names for the particles, available in the EoS. Must strictly match with predefined list ["Neutron", "Proton", "Electron", "Muon", "Tau", "Uquark", "Dquark", "Squark"]. If provided, the instantiator assumes user would like to invoke cooling functionality, hence it requires more provided variables to run.
- `"Quantities"` **:** A linker between EoS columns and the way the supply data to RHM.
    - `"BaryonicDensity"` **:** total baryonic density at the point.
        - `"Units"` (string/double, required) **:** Conversion factor to natural units (GeV powers). It must either be supplied as a choice from ["Gev3", "Fm-3"], or as an actual multiplier.
        ```{note}
        Example : if your units are $m^{-3}$, the factor must be set to $\frac{1}{\text{m}^3\cdot \text{GeV}^3} \approx 7.72\text{E-}48$, since $n_b[\text{GeV}^3] = \dfrac{n_b[\text{m}^{-3}] \cdot \text{m}^{-3}}{\text{GeV}^3} = n_b[\text{m}^{-3}] \cdot \frac{1}{\text{m}^3\cdot \text{GeV}^3}$
        ```
        - `"Column"` (uint, required) **:** Column number with baryonic density, counting from 0th.
        - `"Low"` (double, required) **:** Smallest accessible $n_b$ within the datafile.
        - `"Upp"` (double, required) **:** Biggest accessible $n_b$ within the datafile.
        - `"CoreLimit"` (double, required) **:** First $n_b$ upon entering core.
        - `"CrustLimit"` (double) **:** Last $n_b$ upon leaving crust. If not supplied, "CoreLimit" value is substituted.
    - `"EnergyDensity"` **:** Total energy density at the point.
        - `"Units"` (string/double, required) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevOverFm3", "GOverCm3"], or specify an actual multiplier.
        - `"Column"` (uint, required) **:** Column number with energy density, counting from 0th.
        - `"Low"` (double/string) **:** Smallest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"Upp"` (double/string) **:** Biggest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"CoreLimit"` (double/string) **:** First $\rho_E$ upon entering core. Deduced automatically if specified as "Deduce" or if left blank.
    - `"Pressure"` **:** Total pressure at the point.
        - `"Units"` (string/double, required) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevOverFm3", "DyneOverCm2"], or specify an actual multiplier.
        - `"Column"` (uint, required) **:** Column number with pressure, counting from 0th.
        - `"Low"` (double/string) **:** Smallest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"Upp"` (double/string) **:** Biggest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
    - `"BarionicDensities"` **:** Number densities per each particle at the point.
        - `"$PARTICLE_NAME"` **:** Properties, related to a specific particle. **Substitute** the entry name with actual (fermionic only) species among provided in "Particles" array one by one.
            - `"Column"` (uint, required*) **:** Column number with number density for given particle, counting from 0th.
            - `"ProvidedAs` (string) **:** The way the number density is provided. Choose from ["Density", "DensityFraction", "KFermi"]. "Density" mode expects actual number density, "DensityFraction" mode expects ratio of the number density to the total baryonic density, "KFermi" mode expects Fermi momentum of the particle. Defaults to "Density".
            ```{note}
            As meaningless as it is, in its current state "Units" are shared between particle entries and cannot be overriden. This leads to limiting all particles to the same units. Something to work on as soon.
            ```
        - `"Units"` (string/double) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev3", "Fm-3", "DimLess", "Gev", "Fm-1"], or specify an actual multiplier. Defaults to the "Units" entry of "BaryonicDensity".
        
        ```{warning}
        The dimensionality are different for different provision modes; it is user's responsibility to provide correct units.
        ```
    - `"EffectiveMasses"` **:** [Effective masses](https://en.wikipedia.org/wiki/Effective_mass_(solid-state_physics)) per each particle at the point.
        - `"$PARTICLE_NAME"` **:** Properties, related to a specific particle. **Substitute** the entry name with actual (fermionic only) species among provided in "Particles" array one by one.
            - `"Column"` (uint) **:** Column number with effective mass for given particle, counting from 0th.
            - `"ProvidedAs"` (string, required*) **:** The way the effective mass is provided. Choose from ["FermiEnergy", "EffectiveMass"]. "FermiEnergy" mode incurs the effective mass via relativistic formula with Fermi momentum (occasionally applicable for light particles) and "EffectiveMass" mode expects actual effective mass. If "FermiEnergy" mode is chosen, the "Column" entry is ignored.
            
            ```{note}
            As meaningless as it is, in its current state "Units" are shared between particle entries and cannot be overriden. This leads to limiting all particles to the same units. Something to work on as soon.
            ```
