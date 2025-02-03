/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "Hyperreduction.h"
#include "DEIM.h"
#include "QDEIM.h"
#include "GNAT.h"
#include "S_OPT.h"

#include "linalg/Matrix.h"
#include "utils/Utilities.h"

namespace CAROM {

Hyperreduction::Hyperreduction(const char* sampling_type)
{
    auto iter = SamplingTypeMap.find(sampling_type);
    CAROM_VERIFY(iter != std::end(SamplingTypeMap));

    samplingType = iter->second;
}

void Hyperreduction::ComputeSamples(const std::shared_ptr<Matrix> & f_basis,
                                    int num_f_basis_vectors_used,
                                    std::vector<int>& f_sampled_row,
                                    std::vector<int>& f_sampled_rows_per_proc,
                                    Matrix& f_basis_sampled_inv,
                                    int myid,
                                    int num_procs,
                                    const int num_samples_req,
                                    std::vector<int> *init_samples,
                                    bool qr_factorize)
{
    switch (samplingType)
    {
    case deim:
        CAROM_VERIFY(num_samples_req == f_basis->numColumns());
        DEIM(*f_basis,
             num_f_basis_vectors_used,
             f_sampled_row,
             f_sampled_rows_per_proc,
             f_basis_sampled_inv,
             myid, num_procs);
        return;
    case gnat:
        GNAT(*f_basis,
             num_f_basis_vectors_used,
             f_sampled_row,
             f_sampled_rows_per_proc,
             f_basis_sampled_inv,
             myid, num_procs,
             num_samples_req,
             init_samples);
        return;
    case qdeim:
        QDEIM(*f_basis,
              num_f_basis_vectors_used,
              f_sampled_row,
              f_sampled_rows_per_proc,
              f_basis_sampled_inv,
              myid, num_procs,
              num_samples_req);
        return;
    case sopt:
        S_OPT(*f_basis,
              num_f_basis_vectors_used,
              f_sampled_row,
              f_sampled_rows_per_proc,
              f_basis_sampled_inv,
              myid, num_procs,
              num_samples_req,
              init_samples,
              qr_factorize);
        return;
    default:
        CAROM_ERROR("Sampling type not supported");
    }
}

}
