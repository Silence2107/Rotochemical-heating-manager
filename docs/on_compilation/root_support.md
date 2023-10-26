# ROOT support

There is a seamless integration of [CERN ROOT](https://root.cern/) libraries, used for graphics/I/O inessential purposes only. Make automatically deduces, whether your system carries ROOT, and thereafter makes sure most programs will run regardless, though with different functionality.

For example, on a system with ROOT installed, cooling curve executable gives more extensive choice of console options:
```
> bin/project/Cooling/cooling_curve.out --help
usage: cooling_curve [OPTIONS]...

Evaluates surface temperature time dependency based on EoS

keyword arguments:
  --help, -h            : show this help
  --inputfile VALUE     : json input file path (required)
  --pdf_path VALUE      : pdf output file path (optional, default: Cooling.pdf)
  --rootfile_path VALUE : root output file path (optional, default: None)

Argparse powered by SiLeader
```

Therefore, even running the executable as is
```bash
bin/project/Cooling/cooling_curve.out --inputfile presupplied/Inputfile/RHMconfig.json
```
will not only produce the console tabulation, but also a pdf preview {numref}`apr-example-cooling`.
```{figure} ../plots/apr_example_cooling.jpg
---
name: apr-example-cooling
scale: 30%
---
Cooling curve from RHM with presupplied settings. $2M_{\odot}$, APR4 EoS
```

Another advantage is that ROOT allows for quick I/O with high-level graphic objects, like `TGraph`. If you provide <span style="color:blue">_project/Cooling/nonequilibrium\_profiles.cxx_</span> with `--rootfile_path` value, not only will it produce a preview {numref}`apr-example-profile`,
```{figure} ../plots/apr_example_profiles.jpg
---
name: apr-example-profile
scale: 30%
---
Cooling profile evolving with time from RHM with presupplied settings. $2M_{\odot}$, APR4 EoS
```
but this same preview will be saved to a ROOT file, which can be later opened and manipulated {numref}`apr-example-rootfile`.
```{figure} ../plots/apr_example_rootfile.jpg
---
name: apr-example-rootfile
scale: 30%
---
ROOT file preview for nonequilibrium profile with same settings.
```
Contents of ROOT files vary between programs and occasionally running options, but never exceed the maximum verbosity of the console output.

Be advised that some programs may not have ROOT-free executable available, however we make sure that key functionality is there regardless. In case you attempt to compile such a program without ROOT, `make` will explicitly notify you about it.