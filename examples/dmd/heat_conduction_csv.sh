#!/bin/bash

#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p pbatch
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate

case "$(uname -s)" in
    Linux*)
		  COMMAND="srun -p pdebug"
			MACHINE="Linux";;
    Darwin*)
		  COMMAND="mpirun -oversubscribe"
			MACHINE="Darwin";;
    *)
			echo "The CSV databse test can only run on Linux and MAC."
			exit 1
esac

rm -rf parameters.txt hc_list hc_data run/hc_data

$COMMAND -n 1 ./parametric_heat_conduction -r 0.40 -offline -rdim 16 -csv -out hc_data/hc_par1 > hc_par1.log
$COMMAND -n 1 ./parametric_heat_conduction -r 0.45 -offline -rdim 16 -csv -out hc_data/hc_par2 > hc_par2.log
$COMMAND -n 1 ./parametric_heat_conduction -r 0.55 -offline -rdim 16 -csv -out hc_data/hc_par3 > hc_par3.log
$COMMAND -n 1 ./parametric_heat_conduction -r 0.60 -offline -rdim 16 -csv -out hc_data/hc_par4 > hc_par4.log
$COMMAND -n 1 ./parametric_heat_conduction -r 0.50 -online  -predict -csv -out hc_data/hc_par5 > hc_par5.log

mv -f run/hc_data .
mkdir hc_list

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

$COMMAND -n 1 ./local_tw_csv -o local_csv_serial -rdim 16 -dtc 0.01 > local_csv_serial.log

$COMMAND -n 1 ./local_tw_csv -o local_csv_tw -rdim 16 -nwinsamp 25 -dtc 0.01 > local_csv_tw.log

$COMMAND -n 1 ./parametric_tw_csv -offline -predict -o parametric_csv_serial -rdim 16 -dtc 0.01  > parametric_csv_serial.log
$COMMAND -n 1 ./parametric_tw_csv -online  -predict -o parametric_csv_serial -rdim 16 -dtc 0.01 >> parametric_csv_serial.log

$COMMAND -n 1 ./parametric_tw_csv -offline -predict -o parametric_csv_tw -rdim 16 -nwinsamp 25 -dtc 0.01  > parametric_csv_tw.log
$COMMAND -n 1 ./parametric_tw_csv -online  -predict -o parametric_csv_tw -rdim 16 -nwinsamp 25 -dtc 0.01 >> parametric_csv_tw.log
