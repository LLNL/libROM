/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Compile with: make parametric_tw_csv
//
// Generate CSV or HDF database on heat conduction with either
// heat_conduction_csv.sh or heat_conduction_hdf.sh (HDF is more efficient).
//
// =============================================================================
//
// Parametric serial DMD command (for HDF version, append -hdf):
//   parametric_tw_csv -o hc_parametric_serial -rdim 16 -dtc 0.01 -offline
//   parametric_tw_csv -o hc_parametric_serial -rdim 16 -dtc 0.01 -online
//
// Final-time prediction error (Last line in run/hc_parametric_serial/dmd_par5_prediction_error.csv):
//   0.0012598331433506
//
// Parametric time windowing DMD command (for HDF version, append -hdf):
//   parametric_tw_csv -o hc_parametric_tw -nwinsamp 25 -dtc 0.01 -offline
//   parametric_tw_csv -o hc_parametric_tw -nwinsamp 25 -dtc 0.01 -online
//
// Final-time prediction error (Last line in run/hc_parametric_tw/dmd_par5_prediction_error.csv):
//   0.0006507358659606
//
// =============================================================================
//
// Description: Parametric time windowing DMD on general CSV datasets.
//
// User specify file locations and names by -list LIST_DIR -train-set TRAIN_LIST -test-set TEST_LIST -data DATA_DIR -var VAR_NAME -o OUT_DIR
//
// File structure:
// 1. LIST_DIR/TRAIN_LIST.csv             -- each row specifies one training DATASET
// 2. LIST_DIR/TEST_LIST.csv              -- each row specifies one testing DATASET
// 3. LIST_DIR/DATASET.csv                -- each row specifies the suffix of one STATE in DATASET
// 4. DATA_DIR/dim.csv                    -- specifies the dimension of VAR_NAME
// 5. DATA_DIR/DATASET/tval.csv           -- specifies the time instances
// 6. DATA_DIR/DATASET/STATE/VAR_NAME.csv -- each row specifies one value of VAR_NAME of STATE
// 7. DATA_DIR/DATASET/TEMPORAL_IDX.csv   -- (optional) specifies the first and last temporal index in DATASET
// 8. DATA_DIR/SPATIAL_IDX.csv            -- (optional) each row specifies one spatial index of VAR_NAME
// 9. run/OUT_DIR/indicator_val.csv       -- (optional) each row specifies one indicator endpoint value

#include "mfem.hpp"
#include "algo/DMD.h"
#include "algo/AdaptiveDMD.h"
#include "algo/NonuniformDMD.h"
#include "algo/manifold_interp/SnapshotInterpolator.h"
#include "linalg/Vector.h"
#include "linalg/Matrix.h"
#include "utils/HDFDatabase.h"
#include "utils/CSVDatabase.h"
#include <cmath>
#include <iostream>
#include <limits>

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

using namespace std;
using namespace CAROM;

int main(int argc, char *argv[])
{
    int n_out = 25;
    int n_snap = 11;
    //input times
    vector<double> t{ 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15};
    //Test function from original PCHIP paper
    vector<double> y{10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85}; 
    //sin(x)
    vector<double> y2{0, 0.909297426825682, 0.141120008059867, -0.958924274663138, -0.279415498198926,
                      0.989358246623382, 0.412118485241757, -0.999990206550703, -0.536572918000435,
                      0.990607355694870, 0.650287840157117};

    //query points
    vector<double> tq{0, 0.6000, 1.2000, 1.8000, 2.4000, 3.0000, 3.6000, 4.2000, 4.8000, 5.4000, 6.0000, 6.6000,
                      7.2000, 7.8000, 8.4000, 9.0000, 9.6000, 10.2000, 10.8000, 11.4000, 12.0000, 12.6000, 
                      13.2000, 13.8000, 14.4000, 15.000};
    //pchip reference solution for original test function.
    vector<double> yq1{10.000000000000000, 10.000000000000000, 10.000000000000002, 10.000000000000004, 
                       10.000000000000000, 10.000000000000000, 10.000000000000002, 10.000000000000000, 
                       10.000000000000004, 10.000000000000000, 10.000000000000000, 10.000000000000002,
                       10.000000000000000, 10.000000000000004, 10.102641509433962, 10.500000000000000,
                       11.106230625292378, 12.213163262123807, 14.128630750038976, 27.078413223140512,
                       50.000000000000000, 53.832363636363638, 55.720727272727267, 58.433818181818161,
                       67.055999999999969, 85.000000000000000};
    //pchip reference solution for sin(x)
    vector<double> yq2{0, 0.569748887844739, 0.833039030477175, 0.906694689302138, 0.701592409322499, 0.141120008059867,
                       -0.288488198334352, -0.697095554949408, -0.939878053603591, -0.782971083698006, -0.279415498198926,
                       0.188293444380395, 0.669217685146469, 0.965688937709033, 0.846474732609843, 0.412118485241757,
                       -0.077580693288345, -0.623537711025344, -0.971758328554163, -0.890773577229575, -0.536572918000435,
                       -0.041614069121016, 0.560852411851254, 0.957953731078007, 0.938810668593539, 0.650287840157117};

    std::vector<Vector*> snapshots;
    std::vector<Vector*> out_snapshots;
    std::vector<Vector*> reference_snapshots;
    std::vector<Vector*> times;
    std::vector<Vector*> out_times;
    for(int i = 0; i < t.size(); ++i)
    {
        Vector* temp_t = new Vector(1, false);
        temp_t->item(0) = t[i];
        times.push_back(temp_t);
        Vector* temp_y = new Vector(2,false);
        temp_y->item(0) = y[i];
        temp_y->item(1) = y2[i];
        snapshots.push_back(temp_y);
    }

    for(int i = 0; i < tq.size(); ++i)
    {
        Vector* temp_t = new Vector(1, false);
        temp_t->item(0) = tq[i];
        out_times.push_back(temp_t);
        Vector* temp_y = new Vector(2,false);
        temp_y->item(0) = yq1[i];
        temp_y->item(1) = yq2[i];
        reference_snapshots.push_back(temp_y);
    }

    SnapshotInterpolator* interp = new SnapshotInterpolator();
    
    std::cout << "Beginning base interpolation function" << std::endl;
    out_snapshots = interp->interpolate(times,snapshots,out_times);
    std::cout << "Finished interpolation " << std::endl;
/*
   for(int i = 0; i < out_snapshots.size(); ++i)
    {
        std::cout << "Time " << i << " is " << out_times[i]->item(0) << std::endl;
        std::cout << "Reference at " << i << " is (" << reference_snapshots[i]->item(0) << 
                     "," << reference_snapshots[i]->item(1) << ")" << std::endl;
        std::cout << "Snapshot interpolation at " << i << " is (" << out_snapshots[i]->item(0) << 
                     "," << out_snapshots[i]->item(1) << ")" << std::endl;
    }
 */   
    for(int i = 0; i < out_snapshots.size(); ++i)
    {
        std::cout << "Error at time[" << i << "] = " << out_times[i]->item(0) << " is (" << 
                reference_snapshots[i]->item(0) - out_snapshots[i]->item(0) << "," 
                << reference_snapshots[i]->item(1) - out_snapshots[i]->item(1) << ") " <<  std::endl;
    }


    std::cout << "Beginning variant interpolation function" << std::endl;
    out_snapshots = interp->interpolate(times,snapshots,26,&out_times);
    std::cout << "Finished interpolation " << std::endl;

    for(int i = 0; i < out_snapshots.size(); ++i)
    {
        std::cout << "Error at time[" << i << "] = " << out_times[i]->item(0) << " is (" << 
                reference_snapshots[i]->item(0) - out_snapshots[i]->item(0) << "," 
                << reference_snapshots[i]->item(1) - out_snapshots[i]->item(1) << ") " <<  std::endl;
    }

}