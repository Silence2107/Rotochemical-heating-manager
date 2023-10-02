# What data does program need?

RHM may be a powerful tool, but in a nutshell it cannot resolve the physics without an user telling the properties of the matter. These properties are usually summarized in a so-called equation of state (EoS), and are supplied to RHM as global data. Most EoS are complicated enough to not have analytic formulas, but rather to be presented as a tabulation $P(n_b), \rho_E(n_b), ...$. An example of such tabulation (APR4, {cite}`Akmal_1998`) is to be found under <span style="color:blue">_presupplied/EoS/APR\_EOS\_Acc\_Fe\_RHMstandard.dat_</span>.

While all the libraries in RHM do not explicitly depend on global data, supplied to RHM, they are transferred there via RHM main programs. High-level global settings (datafile interpolators, variable with solver settings etc -- whatever the main programs use) are initialized in <span style="color:blue">_include/instantiator.hpp_</span>. Structure is sketched on {numref}`rhm-structure`.

```{figure} ../plots/RHM_structure.png
---
name: rhm-structure
scale: 30%
---
RHM dependency structure.
```

In what follows all possible ways to instantiate RHM are discussed and covered.
