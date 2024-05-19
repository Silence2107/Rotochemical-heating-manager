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
make release app=project/Cooling/cooling_curve
bin/project/Cooling/cooling_curve.out --help
bin/project/Cooling/cooling_curve.out --inputfile=whatever/input/file/path.json
```

## Environment variables

<span style="color:blue">_Makefile_</span> may be further customized by setting enviromental variables. Among currently supported ones are
    
- **RHM_HAS_ROOT** - if overriden to 1/0, the stable main programs will/will not compile parts of the code that require CERN ROOT library. Default is resolved automatically.

- **RHM_REQUIRES_INPUTFILE** - if overriden to 1/0, the stable main programs will/will not expect an input file to be passed to them. Inputfile is required by default. See [Manual instantiation](../data_input/instantiate_manually.md) for more details, but assume it an advanced way.

## Console options

[Argparse](https://github.com/SiLeader/argparse) is responsible for parsing console options passed to the main programs. Related functionality is revealed, as would be anticipated, by running the program's `--help` option.