# Rotochemical-heating-manager

Toolset for simulating compact star matter with provided nuclear equations of state (EoS). With little restrictions, the toolset is capable of constructing
- Mass Radius diagram
- Cooling curve based on classical cooling theory
- Cooling curve with non-equilibrium (rotochemical) heating enabled
- Various raw additional functionality, including building radial dependencies, high-level variables etc.

## License

RHM is under construction and, until release, is **read-only**.

## Getting started

Once ready to use, the best approach is to familiarize yourself with [documentation](https://silence2107.github.io/Rotochemical-heating-manager).

Here follows a short reference guide:

### Input

In order to setup the RHM to work with your EoS, one must supply the program with global data under include/instantiator.hpp . We provide examples for [presupplied APR4 EoS](https://github.com/Silence2107/Rotochemical-heating-manager/blob/main/presupplied/Inputfile/RHMconfig.json), alongside providing users with [inputfile manual](https://silence2107.github.io/Rotochemical-heating-manager/data_input/what_data_program_need.html).

### Compilation

Any (C++) RHM program that uses the libraries provided here, may be compiled in release/debug configuration by running (assume top directory)
```
make release/debug app=path/to/cxxfile_no_ext
```
If successful, the binary will be put under `bin/path/to/cxxfile_no_ext.out`. See [link](https://silence2107.github.io/Rotochemical-heating-manager/on_compilation/compilation_usage.html) for more details.

### Plots

RHM may output raw text tabulations, alongside (optional) [CERN ROOT](https://root.cern/) outputs. Refer to [link](https://github.com/Silence2107/Rotochemical-heating-manager/tree/main/project/Processing) for quick plotting scripts.
