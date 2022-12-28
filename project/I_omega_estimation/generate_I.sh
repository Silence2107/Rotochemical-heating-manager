#!/bin/bash

filepath="../../data/ist_ns_I_omega_i_lin_interp.txt"
executable="../../bin/estimator_for_I_omegai.out"
touch $filepath
echo "M/Msol I_e I_n I_p" > $filepath

for val in {1..99}
do
    frac=$(python3 -c "print(float($val)/100.0)")
    echo "./$executable $frac"
    echo $(./$executable $frac) >> $filepath
done