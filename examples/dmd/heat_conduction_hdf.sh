#!/bin/bash

#SBATCH -N 1
#SBATCH -t 0:05:00
#SBATCH -p pbatch
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

rm -rf parameters.txt hc_list hc_data run/hc_data

./parametric_heat_conduction -r 0.40 -save -hdf -out hc_data/hc_par1 -no-vis > hc_par1.log
./parametric_heat_conduction -r 0.45 -save -hdf -out hc_data/hc_par2 -no-vis > hc_par2.log
./parametric_heat_conduction -r 0.55 -save -hdf -out hc_data/hc_par3 -no-vis > hc_par3.log
./parametric_heat_conduction -r 0.60 -save -hdf -out hc_data/hc_par4 -no-vis > hc_par4.log
./parametric_heat_conduction -r 0.50 -save -hdf -out hc_data/hc_par5 -no-vis > hc_par5.log

mkdir hc_list
mv -f run/hc_data .

echo "hc_par1,0.40,0.01,0,0"  > hc_list/hc_train_parametric.csv
echo "hc_par2,0.45,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par3,0.55,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par4,0.60,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par5,0.50,0.01,0,0"  > hc_list/hc_train_local.csv
echo "hc_par5,0.50,0.01,0,0"  > hc_list/hc_test.csv
