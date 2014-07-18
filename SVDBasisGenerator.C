#include "SVDBasisGenerator.h"
#include "BasisWriter.h"

namespace CAROM {

SVDBasisGenerator::SVDBasisGenerator(
   const std::string& basis_file_name,
   Database::formats file_format) :
   d_basis_writer(0)
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
}

}
