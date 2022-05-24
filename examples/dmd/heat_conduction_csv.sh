#!/bin/bash

#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p pbatch
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

rm -rf parameters.txt hc_list hc_data run/hc_data

./parametric_heat_conduction -r 0.40 -csv -out hc_data/hc_par1 > hc_par1.log
./parametric_heat_conduction -r 0.45 -csv -out hc_data/hc_par2 > hc_par2.log
./parametric_heat_conduction -r 0.55 -csv -out hc_data/hc_par3 > hc_par3.log
./parametric_heat_conduction -r 0.60 -csv -out hc_data/hc_par4 > hc_par4.log
./parametric_heat_conduction -r 0.50 -csv -out hc_data/hc_par5 > hc_par5.log

mkdir hc_list
mv -f run/hc_data .

awk 'END{print NR}' hc_data/hc_par1/step0/sol.csv > hc_data/dim.csv

mv hc_data/hc_par1/snap_list.csv hc_list/hc_par1.csv
mv hc_data/hc_par2/snap_list.csv hc_list/hc_par2.csv
mv hc_data/hc_par3/snap_list.csv hc_list/hc_par3.csv
mv hc_data/hc_par4/snap_list.csv hc_list/hc_par4.csv
mv hc_data/hc_par5/snap_list.csv hc_list/hc_par5.csv

echo "hc_par1, 0.40, 0.01, 0, 0"  > hc_list/hc_train_parametric.csv
echo "hc_par2, 0.45, 0.01, 0, 0" >> hc_list/hc_train_parametric.csv
echo "hc_par3, 0.55, 0.01, 0, 0" >> hc_list/hc_train_parametric.csv
echo "hc_par4, 0.60, 0.01, 0, 0" >> hc_list/hc_train_parametric.csv
echo "hc_par5, 0.50, 0.01, 0, 0"  > hc_list/hc_train_local.csv
echo "hc_par5, 0.50, 0.01, 0, 0"  > hc_list/hc_test.csv
