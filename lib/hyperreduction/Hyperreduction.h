/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to hyperreduction algorithms.

#ifndef included_Hyperreduction_h
#define included_Hyperreduction_h


#include <string>
#include <unordered_map>

#include <vector>

namespace CAROM {

enum SamplingType
{
    deim,      // Default, DEIM
    gnat,      // GNAT
    qdeim,     // QDEIM
    sopt       // S-OPT
};

static std::unordered_map<std::string, SamplingType> SamplingTypeMap =
{
    {"deim", deim},
    {"gnat", gnat},
    {"qdeim", qdeim},
    {"sopt", sopt}
};

class Matrix;
class Vector;

class Hyperreduction
{
public:

    Hyperreduction(SamplingType stype) :
        samplingType(stype)
    { }

    Hyperreduction(const char* sampling_type);

    void SetSamplingType(SamplingType stype)
    {
        samplingType = stype;
    }

    void ComputeSamples(const Matrix* f_basis,
                        int num_f_basis_vectors_used,
                        std::vector<int>& f_sampled_row,
                        std::vector<int>& f_sampled_rows_per_proc,
                        Matrix& f_basis_sampled_inv,
                        int myid,
                        int num_procs,
                        const int num_samples_req = -1,
			bool precond=false,
			Vector* K=nullptr,
			std::vector<int> *init_samples=NULL,
                        bool qr_factorize = false,
			bool init_trunc = false);

private:
    SamplingType samplingType;
};

}
#endif
