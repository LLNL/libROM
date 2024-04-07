/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMD algorithm on the given snapshot matrix
//              with uniform sampling time steps, using incremenatal SVD
//              to update the DMD matrix.

#ifndef included_IncrementalDMD_h
#define included_IncrementalDMD_h

#include "DMD.h"
#include "linalg/BasisGenerator.h"
#include "linalg/svd/IncrementalSVDBrand.h"

namespace CAROM {

/**
 * Class IncrementalDMD implements the standard DMD algorithm on the
 * given snapshot matrix with uniform sampling time steps, using
 * incremental SVD for update.
 */
class IncrementalDMD : public DMD
{
public:

    /**
     * @brief Constructor. Basic DMD with uniform time step size.
     *
     * @param[in] dim              The full-order state dimension.
     * @param[in] dt               The dt between samples.
     * @param[in] svd_options       Options for SVD.
     * @param[in] alt_output_basis Whether to use the alternative basis for
     *                             output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset     The state offset.
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

    /**
     * @brief Make a prediction for just one dt.
     */
    Vector* predict_dt(Vector* u) override;


protected:

    void updateDMD(const Matrix* f_snapshots);

    BasisGenerator* bg = NULL;
    IncrementalSVDBrand* svd = NULL;

private:

};

}

#endif
