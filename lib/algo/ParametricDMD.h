/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the ParametricDMD algorithm on the given snapshot matrix. The
//              implemented dynamic mode decomposition algorithm is derived from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              This algorithm also works in the case that the first sample does
//              not start from t = 0.0 by incorporating a time offset.

#ifndef included_ParametricDMD_h
#define included_ParametricDMD_h

#include "DMD.h"

#include <vector>
#include <cstring>

namespace CAROM {

class Vector;

/**
 * @brief Constructor.
 *
 * @param[in] parameter_points The parameter points.
 * @param[in] dmds The DMD objects associated with
 *                       each parameter point.
 * @param[in] desired_point The desired point to create a parametric DMD at.
 * @param[in] rbf       The RBF type ("G" == gaussian, "MQ" == multiquadric,
 *                      "IQ" == inverse quadratic, "IMQ" == inverse
 *                      multiquadric)
 * @param[in] interp_method  The interpolation method type ("LS" == linear solve,
 *                      "IDW" == inverse distance weighting, "LP" == lagrangian polynomials)
 * @param[in] epsilon   The RBF parameter that determines the width of
                        influence.
 */
DMD* getParametricDMD(std::vector<Vector*> parameter_points,
                      std::vector<DMD*> dmds,
                      Vector* desired_point,
                      std::string rbf = "G",
                      std::string interp_method = "LS",
                      double epsilon = 1.0);

/**
 * @brief Constructor.
 *
 * @param[in] parameter_points The parameter points.
 * @param[in] dmd_paths The paths to the saved DMD objects associated with
 *                       each parameter point.
 * @param[in] desired_point The desired point to create a parametric DMD at.
 * @param[in] rbf       The RBF type ("G" == gaussian, "MQ" == multiquadric,
 *                      "IQ" == inverse quadratic, "IMQ" == inverse
 *                      multiquadric)
 * @param[in] interp_method  The interpolation method type ("LS" == linear solve,
 *                      "IDW" == inverse distance weighting, "LP" == lagrangian polynomials)
 * @param[in] epsilon   The RBF parameter that determines the width of
                        influence.
 */
DMD* getParametricDMD(std::vector<Vector*> parameter_points,
                      std::vector<std::string> dmd_paths,
                      Vector* desired_point,
                      std::string rbf = "G",
                      std::string interp_method = "LS",
                      double epsilon = 1.0);
}

#endif
