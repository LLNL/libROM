#!/bin/bash

#SBATCH -N 1
#SBATCH -t 0:05:00
#SBATCH -p pdebug
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

rm -rf parameters.txt dmd_list dmd_data run/dmd_data

./parametric_heat_conduction -r 0.40 -save -csv -out dmd_data/dmd_par1 -no-vis > hc_par1.log
./parametric_heat_conduction -r 0.45 -save -csv -out dmd_data/dmd_par2 -no-vis > hc_par2.log
./parametric_heat_conduction -r 0.55 -save -csv -out dmd_data/dmd_par3 -no-vis > hc_par3.log
./parametric_heat_conduction -r 0.60 -save -csv -out dmd_data/dmd_par4 -no-vis > hc_par4.log
./parametric_heat_conduction -r 0.50 -save -csv -out dmd_data/dmd_par5 -no-vis > hc_par5.log

mkdir dmd_list
mv -f run/dmd_data .

awk 'END{print NR}' dmd_data/dmd_par1/step0/sol.csv > dmd_data/dim.csv

mv dmd_data/dmd_par1/snap_list.csv dmd_list/dmd_par1.csv
mv dmd_data/dmd_par2/snap_list.csv dmd_list/dmd_par2.csv
mv dmd_data/dmd_par3/snap_list.csv dmd_list/dmd_par3.csv
mv dmd_data/dmd_par4/snap_list.csv dmd_list/dmd_par4.csv
mv dmd_data/dmd_par5/snap_list.csv dmd_list/dmd_par5.csv

echo "dmd_par1,0.40,0.01,0,0"  > dmd_list/dmd_train_parametric.csv
echo "dmd_par2,0.45,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "dmd_par3,0.55,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "dmd_par4,0.60,0.01,0,0" >> dmd_list/dmd_train_parametric.csv
echo "dmd_par5,0.50,0.01,0,0"  > dmd_list/dmd_train_local.csv
echo "dmd_par5,0.50,0.01,0,0"  > dmd_list/dmd_test.csv
