# Rotochemical-heating-manager

Toolset for simulating compact star matter with provided nuclear equations of state (EoS).

The toolset is capable of constructing
- Mass Radius diagram
- Cooling curve based on classical cooling theory
- Cooling curve with non-equilibrium (rotochemical) heating enabled

## Input

In order to setup the EoS, one must supply the program with global data under include/instantiator.hpp . We provide examples for presupplied APR4 EoS, alongside providing users with [inputfile manual](https://silence2107.github.io/Rotochemical-heating-manager/data_input/what_data_program_need.html).

## Compilation

Any program that uses the libraries provided here, may be compiled in release/debug configuration by running (assume top directory)
```
make release/debug app=path/to/cxxfile_no_ext
```
If successful, the binary will be put under `bin/path/to/cxxfile_no_ext.out`.

### ROOT support

The toolset heavily relies on ROOT framework for graphics. Wherever possible, we decouple the dependence on ROOT, however recentest code is not guaranteed to compile without it. In case you do not have a system-wide installation, you can still run make on the scripts that do not depend on ROOT.
