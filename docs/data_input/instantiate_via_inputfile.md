# Instantiate via inputfile

The instantiator (C++ header) is uniquely provided under <span style="color:blue">_include/instantiator.hpp_</span>, and an example of inputfile lies under <span style="color:blue">_presupplied/APR4/RHMconfig.json_</span>.

The inputfile is formed in [JSON](https://github.com/nlohmann/json) format, namely it has nested dictionary structure. An entry in this dictionary is considered commented if there is a hash key in front of its name. Names follow camel-case convention, except for "EoS" abbreviation.

Using inputfile imposes several restrictions on your modelling. Some of the restrictions might be advisory, but it is advisable to not probe your luck:

- EoS is to be formed in columnar tabulated fashion, with all quantities having their separate column. Quantities, irrelevant in the context (e.g. s-quark density in crust) must be appropriately filled (often zero-padded). The rows must be sorted along $n_b$ column.
- EoS datafile must have barionic density $n_b$ as its only parameter characterising all other quantities.
- Interpolation kind is shared among all columns.
- Quantities, such as pressure, energy density, barionic density, per-particle density and few others may only be supplied as columns, if the program requires them for calculation.
- All callables' logic is limited. For example, equilibration condition is limited to "how much time has passed" and "how different crust temperature is from core temperature".

See <span style="color:blue">_presupplied/APR4/APR_EOS_Acc_Fe_RHMstandard.dat_</span> for a reference.

In what follows we discuss all the possible settings available in the predefined instantiator. We mark recipients of each setting with <span style="color:red">red</span> color. If a given setting is not mandatory to set, it is explicitly marked as such.