# Instantiate manually (advanced)

Manual instantiating implies directly putting global data into an instantiator, without relying on the external inputfile. An example is given under <span style="color:blue">_presupplied/EoS/instantiator\_apr4.hpp_</span>.

**Pros:**
- Most comfortable when adding new features
- Grants higher flexibility over the program (this may include putting whatever logic into callables, supplementing EoS datafile-free, enabling EoS dependence on variables beyond $n_b$ etc)
- Without limitations may also seek supply from inputfile

**Cons:**
- Requires recompilation each time changes are made to the supplied data
- Requires knowledge regarding what variables are essential and how they must be supplied
- Requires language skills