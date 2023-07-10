/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMD algorithm on the given snapshot matrix
//              with uniform sampling time steps, using incremenatal SVD
//              to update the DMD matrix.
//              Instead of approximating the discrete dynamics, this algorithm
//              approximates the continuous dynamics linearly.
//              This algorithm also works in the case that the first sample does
//              not start from t = 0.0 by incorporating a time offset.

#ifndef included_IncrementalDMD_h
#define included_IncrementalDMD_h

#include "DMD.h"
#include "linalg/BasisGenerator.h"

namespace CAROM {

/**
 * Class IncrementalDMD implements the standard DMD algorithm on the
 * given snapshot matrix with uniform sampling time steps, using
 * incremental SVD for update.
 * Instead of linearly approximating the discrete dynamics
 * x(t+dt) = Ax(t) in the original DMD, this algorithm approximates
 * the continuous dynamics linearly by dx/dt = Ax.
 */
class IncrementalDMD : public DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim               The full-order state dimension.
     * @param[in] alt_output_basis  Whether to use the alternative basis for  
     *                              output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset      The state offset.
     * @param[in] derivative_offset The derivative offset.
     */
    IncrementalDMD(int dim,
		   double dt,
		   Options svd_options,
		   std::string svd_base_file_name,
                   bool alt_output_basis = false,
                   Vector* state_offset = NULL);

    /**
     * @brief Destroy the IncrementalDMD object
     */
    ~IncrementalDMD();

    /**
     * @brief Load the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to load the
     *                           database to.
     */
    void load(std::string base_file_name) override;

    /**
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the filename to save the
     *                           database to.
     */
    void save(std::string base_file_name) override;

    /**
     * @param[in] k               The number of modes to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
                                  of W is linearly independent with W0.
     */
    void train(int k,
               const Matrix* W0 = NULL,
	       double linearity_tol = 0.0) override;


protected:

    void updateDMD(const Matrix* f_snapshots);
    
    BasisGenerator* bg = NULL;

private:

};

}

#endif
