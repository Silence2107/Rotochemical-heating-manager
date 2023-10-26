# Instantiate manually (advanced)

<span style="color:green">_This section is essentially useless unless you develop your own contribution to the code._</span>

Manual instantiating implies directly putting global data into an instantiator, without relying on the external inputfile. An example is given under <span style="color:blue">_presupplied/EoS/instantiator\_apr4.hpp_</span>.

**Pros:**
- Most comfortable when adding new features
- Grants higher flexibility over the program (this may include putting whatever logic into callables, supplementing EoS datafile-free, enabling EoS dependence on variables beyond $n_b$ etc)
- Without limitations may also seek supply from inputfile

**Cons:**
- Requires recompilation each time changes are made to the supplied data
- Requires knowledge regarding what variables are essential and how they must be supplied
- Requires language skills

In case you for any reason prefer to store all the data purely within the instantiator, there is no requirement neither for inputfile nor for `instantiator::instantiate_system` implementation. Such a program may be built with 
```bash
make app=any/app/path RHM_REQUIRES_INPUTFILE=0
```
and then ran normally without `--inputfile` key. 

To get the best of both worlds, however, you may prefer to rely on manual instantiation only _partially_, still relying on inputfile with most variables, while some high-level ones are explicitly hard-coded within the instantiator.

```{note}
We do not provide an explicit algorithm for manual instantiation as of now. This essentially dublicates the [Inputfile Structure](../data_input/instantiate_via_inputfile.md) section in a less front-end way, so our best advice is to compare the inputfile entries with instantiator variables.
```