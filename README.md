# Rotochemical-heating-manager

Toolset for simulating compact star matter with provided nuclear equations of state (EoS).

The toolset is capable of constructing
- Mass Radius diagram
- Cooling curve based on classical cooling theory
- Cooling curve with non-equilibrium (rotochemical) heating enabled

## Input

In order to setup the EoS, one must fill in the inputfile under include/inputfile.hpp . Currently we have no definitive guide on how to construct one, so at best you can replicate the one we provide for APR4 EoS.

## Compilation

Any program that uses the libraries provided here, may be compiled in release/debug configuration by running (assume top directory)
```
make release/debug app=path/to/cxxfile_no_ext
```
If successful, the binary will be put under `bin/path/to/cxxfile_no_ext.out`.
