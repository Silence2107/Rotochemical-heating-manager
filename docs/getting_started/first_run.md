# First run

## Overview

RHM follows standard <span style="color:blue">_bin/include/src_</span> structure, with program's logic assembled under <span style="color:blue">_project/_</span> folder. <span style="color:blue">_3rd-party/_</span> folder contains several IO libraries (ready-to-use), and <span style="color:blue">_presupplied/_</span> contains presupplied data. <span style="color:blue">_Makefile_</span> is also supplied for compilation.

Since there are many different main programs under <span style="color:blue">_project/_</span>, some of which may run with settings provided by the user and some may not, making executables is yet user's responsibility.

## Test run

To perform a test run with presupplied data, run make on some typical script, e.g.

```bash
make release app=project/M-R_diagram/m_r_diagram
```

If successful, the corresponding binary will be put under <span style="color:blue">_bin/_</span>, following same path as it took to the app. In this case, executing binary would work like follows

```bash
bin/project/M-R_diagram/m_r_diagram.out --help
```

`--help` invokes manual message for all standardized RHM programs.

To finally see whether physics is in order on your machine, this binary (most of the time) must be supplemented with inputfile. Let's run it by providing the one under <span style="color:blue">_presupplied/APR4/RHMconfig.json_</span>:

```bash
bin/project/M-R_diagram/m_r_diagram.out --inputfile presupplied/APR4/RHMconfig.json
```

**Expected output :**
- $(M, R)$ pairs of order $(2 \text{M}_\odot, 10 \text{km})$ based on various center pressure fractions.
    - phenomenal work

**Unexpected output :**
- $(M, R)$ pairs of $(0 \text{M}_\odot, 0 \text{km})$ 
    - inputfile is not supplied
- "(..) Cannot open file (..)"
    - inputfile path is supplied, but is not recognized as valid
    - Check path's spelling against <span style="color:blue">_presupplied/APR4/RHMconfig.json_</span>
    - If correct, make sure there exists a valid file under <span style="color:blue">_presupplied/APR4/APR\_EOS\_Acc\_Fe\_RHMstandard.dat_</span>
    - If exists, make sure there's an entry in <span style="color:blue">_presupplied/APR4/RHMconfig.json_</span> `["EoSSetup"]["Datafile"]["Path"]` that leads to the file above. Specify absolute path if in doubts.
- "keyword argument (..) must have value."
    - program's key is supplied, but actual value to it is not
- Whatever else happened
    - Tell me I messed

## What's next?

- Investigate the <span style="color:blue">_project/_</span> folder for more programs, to see which you may find useful. Stable ones are inspectable via `--help` flag.
- Further learn about RHM main programs running settings
- Learn about how RHM processes input data
- Supply it with your own EoS!