# EoSSetup

## Overview

This entry contains all relevant information regarding how program is supposed to interpret the supplied EoS. In particular, one specifies here EoS datafile format, interpolation kind, column to quantity mapping, units, particle species etc. In what follows we discuss all the possible settings under this entry.

In order to refer to an example, see <span style="color:blue">_presupplied/Inputfile/RHMconfig.json_</span> and <span style="color:blue">_presupplied/EoS/APR_EOS_Acc_Fe_RHMstandard.dat_</span>.

## Description

- `"Datafile"` **:** Contains all relevant knowledge regarding EoS datafile structure.
    - `"Path"` (string, required) **:** valid path to the EoS datafile.
    - `"Rows"` (uint pair) **:** [first, last) row in the inputfile to consider (counting from 0th). Assumes whole file, if not specified. Last row may be specified as 0 to indicate the last row in the file.
    - `"Columns"` (uint pair) **:** [first, last) row in the inputfile to consider (counting from 0th). Deduces from the first row, if not specified. Last column may be specified as 0 to indicate the last column in the file.
    - `"Interpolation"` (string) **:** Interpolation kind to be used across all columns. Choose from ["Linear", "Cubic"], with "Linear" being default. 
    ```{warning}
    Using "Cubic" mode may lead to issues with positivity of some functions, so it is at the moment advised to instead populate the datafile with more points, yet to use "Linear" mode.
    ```
    ```{warning}
    Be advised that any uint entry (a.k.a. `size_t` in C++) will happily consume negative values(causing overflow), if supplied as such. Since it is usually not the intended behavior, though silently allowed, make sure you supply meaningful values. This applies to **all uint entries** in the inputfile, as well as within the instantiator.
    ```
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
        - `"Units"` (string/double, required) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevFm-3", "GCm-3"], or specify an actual multiplier.
        - `"Column"` (uint, required) **:** Column number with energy density, counting from 0th.
        - `"Low"` (double/string) **:** Smallest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"Upp"` (double/string) **:** Biggest accessible $\rho_E$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"CoreLimit"` (double/string) **:** First $\rho_E$ upon entering core. Deduced automatically if specified as "Deduce" or if left blank.
    - `"Pressure"` **:** Total pressure at the point.
        - `"Units"` (string/double, required) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev4", "MevFm-3", "DyneCm-2"], or specify an actual multiplier.
        - `"Column"` (uint, required) **:** Column number with pressure, counting from 0th.
        - `"Low"` (double/string) **:** Smallest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
        - `"Upp"` (double/string) **:** Biggest accessible $P$ within the datafile. Deduced automatically if specified as "Deduce" or if left blank.
    - `"NumberDensities"` **:** Number densities per each particle at the point.
        - `"$PARTICLE_NAME"` **:** Properties, related to a specific particle. **Substitute** the entry name with actual species among provided in "Particles" array one by one.
            - `"Column"` (uint, required*) **:** Column number with number density for given particle, counting from 0th.
            - `"ProvidedAs` (string) **:** The way the number density is provided. Choose from ["Density", "DensityFraction", "KFermi"]. "Density" mode expects actual number density, "DensityFraction" mode expects ratio of the number density to the total baryonic density, "KFermi" mode expects Fermi momentum of the particle. Defaults to "Density".
            - `"Units"` (string/double) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev3", "Fm-3", "DimLess", "Gev", "Fm-1"], or specify an actual multiplier. Defaults to the "Units" entry of "BaryonicDensity".
        
        ```{warning}
        The dimensionality are different for different provision modes; it is user's responsibility to provide correct units.
        ```
    - `"EffectiveMasses"` **:** [Effective masses](https://en.wikipedia.org/wiki/Effective_mass_(solid-state_physics)) per each particle at the point.
        - `"$PARTICLE_NAME"` **:** Properties, related to a specific particle. **Substitute** the entry name with actual (fermion) species among provided in "Particles" array one by one.
            - `"Column"` (uint) **:** Column number with effective mass for given particle, counting from 0th.
            - `"ProvidedAs"` (string, required*) **:** The way the effective mass is provided. Choose from ["FermiEnergy", "EffectiveMass"]. "FermiEnergy" mode incurs the effective mass via relativistic formula with Fermi momentum (occasionally applicable for light particles) and "EffectiveMass" mode expects actual effective mass. If "FermiEnergy" mode is chosen, the "Column" and "Units" entries is ignored.
            - `"Units"` (string/double, required*) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev", "MeV", "NucleonMass"], or specify an actual multiplier. If "FermiEnergy" mode is chosen, "Units" are disregarded.
            
    - `"IonVolumeFraction"` : Ion volume fraction in the crust at the point. Affects neutrino bremhstrahlung in the crust.
        - `"Column"` (uint) **:** Column number with ion volume fraction, counting from 0th.
        - `"ProvidedAs"` (string, required*) **:** The way the ion volume fraction is provided. Choose from ["IonVolumeFraction", "Absent", "ExcludedVolume"]. "IonVolumeFraction" mode expects actual ion volume fraction, "Absent" mode renders the ratio zero and "ExcludedVolume" mode performs calculation in the crust
        $\eta = \text{min}\left[\dfrac{4}{3}\pi (1.1 \text{fm})^3 \dfrac{\rho_E}{m_{\text{nucleon}}}, 1.0\right]$. Defaults to "Absent". "Column" entry is only used in "IonVolumeFraction" mode.
        ```{note}
        Units are dimensionless for all modes.
        ```
        ```{note}
        This quantity could use some citation and explanation. 
        ```
    - `"QuarkSuperconductingGap"` : Superconductive gap for quarks at the point. Noticeably affects all quark cooling channels as $\sim \exp{\left[-\frac{\Delta}{T}\right]}$, with $\Delta$ being the gap and $T$ being local temperature. 
        - `"Column"` (uint) **:** Column number with quark superconducting gap, counting from 0th. If not specified, the gap is assumed to be absent for all quarks.
        - `"Units"` (string/double) **:** Conversion factor to natural units (GeV powers). Choose from ["Gev", "MeV", "Fm-1"], or specify an actual multiplier. 
        ```{note}
        This quantity could use some variety in provision modes.
        ```
- `"Misc"` : Miscellaneous EoS settings, unfit under entries above.
    - `"CrustalEta"` (double, required*) **:** Light element share in the atmosphere {cite}`potekhin1997internal`.
    - `"ProtonSuperfluidity1S0"` (string) **:** Critical temperature model for protons in 1S0 superfluid state. Choose from ["AO", "CCDK", "A", "B", "C", "A2", "HadronToQGP"], with normal fluidity being default.
    - `"NeutronSuperfluidity1S0"` (string) **:** Critical temperature model for neutrons in 1S0 superfluid state. Choose from ["AO", "CCDK", "A", "B", "C", "A2", "HadronToQGP"], with normal fluidity being default.
    - `"NeutronSuperfluidity3P2"` (string) **:** Critical temperature model for neutrons in 3P2 superfluid state. Choose from ["AO", "CCDK", "A", "B", "C", "A2", "HadronToQGP"], with normal fluidity being default.
    ```{note}
    There is an interplay, if one enables both 1S0 and 3P2 superfluidity for neutrons. In this case, the critical temperature will be imposed as 3P2 in the core and 1S0 beyond the core. This could use some citation.
    ```

