#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "rm -rf parameters.txt dmd_list dmd_data run/dmd_data"
    "./parametric_heat_conduction -r 0.40 -save -csv -out dmd_data/dmd_par1 -no-vis"
    "./parametric_heat_conduction -r 0.45 -save -csv -out dmd_data/dmd_par2 -no-vis"
    "./parametric_heat_conduction -r 0.55 -save -csv -out dmd_data/dmd_par3 -no-vis"
    "./parametric_heat_conduction -r 0.60 -save -csv -out dmd_data/dmd_par4 -no-vis"
    "./parametric_heat_conduction -r 0.50 -save -csv -out dmd_data/dmd_par5 -no-vis"
    "mkdir dmd_list"
    "mv -f run/dmd_data ."
    "awk 'END{print NR}' dmd_data/dmd_par1/step0/sol.csv > dmd_data/dim.csv"
    "mv dmd_data/dmd_par1/snap_list.csv dmd_list/dmd_par1.csv"
    "mv dmd_data/dmd_par2/snap_list.csv dmd_list/dmd_par2.csv"
    "mv dmd_data/dmd_par3/snap_list.csv dmd_list/dmd_par3.csv"
    "mv dmd_data/dmd_par4/snap_list.csv dmd_list/dmd_par4.csv"
    "mv dmd_data/dmd_par5/snap_list.csv dmd_list/dmd_par5.csv"
    "echo 'dmd_par1,0.40,0.01,0,0'  > dmd_list/dmd_train_parametric.csv"
    "echo 'dmd_par2,0.45,0.01,0,0' >> dmd_list/dmd_train_parametric.csv"
    "echo 'dmd_par3,0.55,0.01,0,0' >> dmd_list/dmd_train_parametric.csv"
    "echo 'dmd_par4,0.60,0.01,0,0' >> dmd_list/dmd_train_parametric.csv"
    "echo 'dmd_par5,0.50,0.01,0,0'  > dmd_list/dmd_train_local.csv"
    "echo 'dmd_par5,0.50,0.01,0,0'  > dmd_list/dmd_test.csv"
    "./parametric_tw_csv -o hc_local -train-set dmd_train_local -rdim 16 -dtc 0.01 -offline"
    "./parametric_tw_csv -o hc_local -train-set dmd_train_local -rdim 16 -dtc 0.01 -online"
    "./parametric_tw_csv -o hc_local_tw -train-set dmd_train_local -nwinsamp 25 -dtc 0.01 -offline"
    "./parametric_tw_csv -o hc_local_tw -train-set dmd_train_local -nwinsamp 25 -dtc 0.01 -online"
    "./parametric_tw_csv -o hc_parametric -train-set dmd_train_parametric -rdim 16 -dtc 0.01 -offline"
    "./parametric_tw_csv -o hc_parametric -train-set dmd_train_parametric -rdim 16 -dtc 0.01 -online"
    "./parametric_tw_csv -o hc_parametric_tw -train-set dmd_train_parametric -nwinsamp 25 -dtc 0.01 -offline"
    "./parametric_tw_csv -o hc_parametric_tw -train-set dmd_train_parametric -nwinsamp 25 -dtc 0.01 -online"
)
TYPE="DMD"
OFFSET=0

# custom list of output files to compare
FILES=(
    "./run/hc_local/window0/par0_eigenvalue.csv"
    "./run/hc_local/window0/par0_singular_value.csv"
    "./run/hc_local/dmd_par5_prediction_error.csv"
    "./run/hc_local_tw/window1/par0_eigenvalue.csv"
    "./run/hc_local_tw/window1/par0_singular_value.csv"
    "./run/hc_local_tw/dmd_par5_prediction_error.csv"
    "./run/hc_parametric/window0/par0_eigenvalue.csv"
    "./run/hc_parametric/window0/par0_singular_value.csv"
    "./run/hc_parametric/dmd_par5_prediction_error.csv"
    "./run/hc_parametric_tw/window1/par0_eigenvalue.csv"
    "./run/hc_parametric_tw/window1/par0_singular_value.csv"
    "./run/hc_parametric_tw/dmd_par5_prediction_error.csv"
)

run_tests
