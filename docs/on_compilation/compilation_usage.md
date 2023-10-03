# Usage

Make is used to compile any program that uses RHM library. Example in [Getting Started](../getting_started/first_run.md) is rather generic and, in general, executes as follows
```bash
make release/debug/clean app=path_to_cxx_no_extension
```

Release mode is default and it shall prevail in all cases except debugging. If executed as above, the executable lands into
```bash
bin/path_to_cxx_no_extension.out
```

## Example

Imagine we would like to compile a main program under <span style="color:blue">_project/Cooling/cooling_curve.cxx_</span> and execute it. The following commands will do the job
```bash
make app=project/Cooling/cooling_curve
bin/project/Cooling/cooling_curve.out --help
bin/project/Cooling/cooling_curve.out --inputfile=whatever/input/file/path.json
```

## Environment variables

<span style="color:blue">_Makefile_</span> may be further customized by setting enviromental variables. Among currently supported ones are
    
- **RHM_HAS_ROOT** - if overriden to 1/0, the stable main programs will/will not execute parts of the code that require CERN ROOT library. Default is resolved automatically. Overriding this flag <code><it> will only affect what parts of the codes are compiled and run </it></code>, the ROOT libraries inclusion is still up to make to decide.

## Console options

[Argparse](https://github.com/SiLeader/argparse) is responsible for parsing console options passed to the main programs. Related functionality is revealed, as would be anticipated, by running the program's `--help` option.