#!/bin/bash

#SBATCH -N 1
#SBATCH -t 0:05:00
#SBATCH -p pdebug
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

rm -rf parameters.txt dmd_list dmd_data run/dmd_data

./parametric_heat_conduction -r 0.40 -save -hdf -out dmd_data/sim0 -no-vis > hc_par0.log
./parametric_heat_conduction -r 0.45 -save -hdf -out dmd_data/sim1 -no-vis > hc_par1.log
./parametric_heat_conduction -r 0.55 -save -hdf -out dmd_data/sim2 -no-vis > hc_par2.log
./parametric_heat_conduction -r 0.60 -save -hdf -out dmd_data/sim3 -no-vis > hc_par3.log
./parametric_heat_conduction -r 0.50 -save -hdf -out dmd_data/sim4 -no-vis > hc_par4.log

mkdir dmd_list
mv -f run/dmd_data .

echo "0,0.40,0.01,0,0"  > dmd_list/dmd_train_parametric.csv
echo "1,0.45,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "2,0.55,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "3,0.60,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "4,0.50,0.01,0,0"  > dmd_list/dmd_train_local.csv
echo "4,0.50,0.01,0,0"  > dmd_list/dmd_test.csv
