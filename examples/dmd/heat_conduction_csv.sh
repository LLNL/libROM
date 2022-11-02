#!/bin/bash

#SBATCH -N 1
#SBATCH -t 0:05:00
#SBATCH -p pbatch
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

rm -rf parameters.txt hc_list hc_data run/hc_data

./parametric_heat_conduction -r 0.40 -csv -out hc_data/hc_par1 -no-vis > hc_par1.log
./parametric_heat_conduction -r 0.45 -csv -out hc_data/hc_par2 -no-vis > hc_par2.log
./parametric_heat_conduction -r 0.55 -csv -out hc_data/hc_par3 -no-vis > hc_par3.log
./parametric_heat_conduction -r 0.60 -csv -out hc_data/hc_par4 -no-vis > hc_par4.log
./parametric_heat_conduction -r 0.50 -csv -out hc_data/hc_par5 -no-vis > hc_par5.log

mkdir hc_list
mv -f run/hc_data .

awk 'END{print NR}' hc_data/hc_par1/step0/sol.csv > hc_data/dim.csv

mv hc_data/hc_par1/snap_list.csv hc_list/hc_par1.csv
mv hc_data/hc_par2/snap_list.csv hc_list/hc_par2.csv
mv hc_data/hc_par3/snap_list.csv hc_list/hc_par3.csv
mv hc_data/hc_par4/snap_list.csv hc_list/hc_par4.csv
mv hc_data/hc_par5/snap_list.csv hc_list/hc_par5.csv

mv hc_data/hc_par1/numsnap hc_list/numsnap1
mv hc_data/hc_par2/numsnap hc_list/numsnap2
mv hc_data/hc_par3/numsnap hc_list/numsnap3
mv hc_data/hc_par4/numsnap hc_list/numsnap4
mv hc_data/hc_par5/numsnap hc_list/numsnap5

echo "hc_par1,numsnap1,0.40,0.01,0,0"  > hc_list/hc_train_parametric.csv
echo "hc_par2,numsnap2,0.45,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par3,numsnap3,0.55,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par4,numsnap4,0.60,0.01,0,0" >> hc_list/hc_train_parametric.csv
echo "hc_par5,numsnap5,0.50,0.01,0,0"  > hc_list/hc_train_local.csv
echo "hc_par5,numsnap5,0.50,0.01,0,0"  > hc_list/hc_test.csv
