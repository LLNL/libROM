#include "svd_basis_generator.h"
#include "basis_writer.h"

namespace CAROM {

svd_basis_generator::svd_basis_generator(
   const std::string& basis_file_name,
   Database::formats file_format) :
   d_basis_writer(0)
{
   if (!basis_file_name.empty()) {
      d_basis_writer = new BasisWriter(this, basis_file_name, file_format);
   }
}

svd_basis_generator::~svd_basis_generator()
{
   if (d_basis_writer) {
      delete d_basis_writer;
   }
}

}
