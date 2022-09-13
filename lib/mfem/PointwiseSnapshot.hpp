/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifndef POINTWISE_SNAPSHOT
#define POINTWISE_SNAPSHOT

#include "mfem.hpp"
//#include "gslib.hpp"

using namespace mfem;
using namespace std;


namespace CAROM {

class PointwiseSnapshot
{
public:
  PointwiseSnapshot(const int sdim, int *dims_);

  void SetMesh(ParMesh *pmesh);
  void GetSnapshot(ParGridFunction const& f, mfem::Vector & s);

  ~PointwiseSnapshot();

private:
  FindPointsGSLIB *finder{nullptr};

  int npoints;
  int dims[3];
  const int spaceDim;

  Vector domainMin, domainMax;
  Vector xyz;
};

}

#endif // POINTWISE_SNAPSHOT
