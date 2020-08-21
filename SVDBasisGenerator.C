/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract wrapper class for an abstract SVD algorithm and
//              sampler.  This class provides interfaces to each so that an
//              application only needs to instantiate one concrete
//              implementation of this class to control all aspects of basis
//              vector generation.

#include "SVDBasisGenerator.h"
#include "BasisWriter.h"

namespace CAROM {

SVDBasisGenerator::SVDBasisGenerator(
   const std::string& basis_file_name,
   Database::formats file_format) :
   d_basis_writer(0),
   d_basis_reader(0)
{
   if (!basis_file_name.empty()) {
      d_basis_writer = new BasisWriter(this, basis_file_name, file_format);
   }
}

SVDBasisGenerator::~SVDBasisGenerator()
{
   if (d_basis_writer) {
      delete d_basis_writer;
   }
   if (d_basis_reader) {
     delete d_basis_reader;
   }
}

}
