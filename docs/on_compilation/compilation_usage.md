# Usage

Make is used to compile any program that uses RHM library. Example in [Getting Started](../getting_started/first_run.md) is rather generic and, in general, executes as follows
```bash
make release/debug app=path_to_cxx_no_extension
```

Release mode is default and it shall prevail in all cases except debugging. If executed as above, the executable lands into
```bash
bin/path_to_cxx_no_extension.out
```

To compile all the stable programs in release mode, type
```bash
make all -j$(nproc)
```

## Example

Imagine we would like to compile a main program under <span style="color:blue">_project/Cooling/cooling_curve.cxx_</span> and execute it. The following commands will do the job
```bash
> make release app=project/Cooling/cooling_curve RHM_HAS_ROOT=0
g++ project/Cooling/cooling_curve.cxx src/tov_solver.o src/cooling.o src/auxiliaries.o -o bin/project/Cooling/cooling_curve.out -Wall -Wextra  -DRHM_HAS_ROOT=0 -O3

> bin/project/Cooling/cooling_curve.out --help
usage: cooling_curve [OPTIONS]...

Evaluates surface temperature time dependency based on EoS

keyword arguments:
  --help, -h        : show this help
  --inputfile VALUE : json input file path (required)

Argparse powered by SiLeader

> bin/project/Cooling/cooling_curve.out --inputfile presupplied/APR4/RHMconfig.json
M = 1.4 [Ms]
t [years]           Te^inf [K]          L^inf_ph [erg/s]    L^inf_nu [erg/s]
1e-12               6.80188e+06         3.15924e+36         3.27258e+45
2.14988e-12         6.80188e+06         3.15924e+36         3.27258e+45
3.47211e-12         6.80188e+06         3.15924e+36         3.27258e+45
4.99251e-12         6.80188e+06         3.15924e+36         3.27258e+45
6.7408e-12          6.80188e+06         3.15924e+36         3.27257e+45
8.75112e-12         6.80188e+06         3.15924e+36         3.27257e+45
1.10627e-11         6.80188e+06         3.15924e+36         3.27257e+45
1.37208e-11         6.80188e+06         3.15924e+36         3.27257e+45
1.67773e-11         6.80188e+06         3.15924e+36         3.27257e+45
...
```

## Environment variables

<span style="color:blue">_Makefile_</span> may be further customized by setting enviromental variables. Among currently supported ones are
    
- **RHM_HAS_ROOT** - if overriden to 1/0, the stable main programs will/will not compile parts of the code that require CERN ROOT library. Default is resolved automatically.

## Console options

[Argparse](https://github.com/SiLeader/argparse) is responsible for parsing console options passed to the main programs. Related functionality is revealed, as would be anticipated, by running the program's `--help` option.