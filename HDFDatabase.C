/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: The concrete database implementation using HDF5.

#include "HDFDatabase.h"
#include "Utilities.h"

namespace CAROM {

const int HDFDatabase::KEY_DOUBLE_ARRAY = 0;
const int HDFDatabase::KEY_INT_ARRAY = 1;

HDFDatabase::HDFDatabase() :
   d_is_file(false),
   d_file_id(-1),
   d_group_id(-1)
{
}

HDFDatabase::~HDFDatabase()
{
   if (d_is_file) {
      close();
   }

   if (d_group_id != -1) {
      herr_t errf = H5Gclose(d_group_id);
      CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
      CAROM_NULL_USE(errf);
#endif
   }
}

bool
HDFDatabase::create(
   const std::string& file_name)
{
   CAROM_ASSERT(!file_name.empty());
   hid_t file_id = H5Fcreate(file_name.c_str(),
                             H5F_ACC_TRUNC,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
   bool result = file_id >= 0;
   CAROM_ASSERT(result);
   d_is_file = true;
   d_file_id = file_id;
   d_group_id = file_id;

   return result;
}

bool
HDFDatabase::open(
   const std::string& file_name)
{
   CAROM_ASSERT(!file_name.empty());
   hid_t file_id = H5Fopen(file_name.c_str(),
                           H5F_ACC_RDONLY,
                           H5P_DEFAULT);
   bool result = file_id >= 0;
   CAROM_ASSERT(result);
   d_is_file = true;
   d_file_id = file_id;
   d_group_id = file_id;

   return result;
}

bool
HDFDatabase::close()
{
   herr_t errf = 0;
   if (d_is_file) {
      errf = H5Fclose(d_file_id);
      CAROM_ASSERT(errf >= 0);

      if (d_group_id == d_file_id) {
         d_group_id = -1;
      }
      d_file_id = -1;
      d_is_file = false;
   }

   return errf >= 0;
}

void
HDFDatabase::putIntegerArray(
   const std::string& key,
   const int* const data,
   int nelements)
{
   CAROM_ASSERT(!key.empty());
   CAROM_ASSERT(data != 0);
   CAROM_ASSERT(nelements > 0);

   hsize_t dim[] = { static_cast<hsize_t>(nelements) };
   hid_t space = H5Screate_simple(1, dim, 0);
   CAROM_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t dataset = H5Dcreate(d_group_id,
                             key.c_str(),
                             H5T_STD_I32BE,
                             space,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
#else
   hid_t dataset = H5Dcreate(d_group_id,
                             key.c_str(),
                             H5T_STD_I32BE,
                             space,
                             H5P_DEFAULT);
#endif
   CAROM_ASSERT(dataset >= 0);

   herr_t errf = H5Dwrite(dataset,
                          H5T_NATIVE_INT,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          data);
   CAROM_ASSERT(errf >= 0);

   // Write attribute so we know what kind of data this is.
   writeAttribute(KEY_INT_ARRAY, dataset);

   errf = H5Sclose(space);
   CAROM_ASSERT(errf >= 0);

   errf = H5Dclose(dataset);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::putDoubleArray(
   const std::string& key,
   const double* const data,
   int nelements)
{
   CAROM_ASSERT(!key.empty());
   CAROM_ASSERT(data != 0);
   CAROM_ASSERT(nelements > 0);

   hsize_t dim[] = { static_cast<hsize_t>(nelements) };
   hid_t space = H5Screate_simple(1, dim, 0);
   CAROM_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t dataset = H5Dcreate(d_group_id,
                             key.c_str(),
                             H5T_IEEE_F64BE,
                             space,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
#else
   hid_t dataset = H5Dcreate(d_group_id,
                             key.c_str(),
                             H5T_IEEE_F64BE,
                             space,
                             H5P_DEFAULT);
#endif
   CAROM_ASSERT(dataset >= 0);

   herr_t errf = H5Dwrite(dataset,
                          H5T_NATIVE_DOUBLE,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          data);
   CAROM_ASSERT(errf >= 0);

   // Write attribute so we know what kind of data this is.
   writeAttribute(KEY_DOUBLE_ARRAY, dataset);

   errf = H5Sclose(space);
   CAROM_ASSERT(errf >= 0);

   errf = H5Dclose(dataset);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::getIntegerArray(
   const std::string& key,
   int* data,
   int nelements)
{
   CAROM_ASSERT(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(nelements);
#endif

   if (!isInteger(key)) {
      CAROM_ERROR("HDFDatabase::getIntegerArray() error in database\n"
         << "    Key = " << key << " is not an integer array." << std::endl);
   }

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
   CAROM_ASSERT(dset >= 0);

   hid_t dspace = H5Dget_space(dset);
   CAROM_ASSERT(dspace >= 0);

   hsize_t nsel = H5Sget_select_npoints(dspace);
   CAROM_ASSERT(static_cast<int>(nsel) == nelements);

   herr_t errf;
   if (nsel > 0) {
     errf = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
     CAROM_ASSERT(errf >= 0);
   }

   errf = H5Sclose(dspace);
   CAROM_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::getDoubleArray(
   const std::string& key,
   double* data,
   int nelements)
{
   CAROM_ASSERT(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(nelements);
#endif

   if (!isDouble(key)) {
      CAROM_ERROR("HDFDatabase::getDoubleArray() error in database\n"
         << "    Key = " << key << " is not a double array." << std::endl);
   }

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
   CAROM_ASSERT(dset >= 0);

   hid_t dspace = H5Dget_space(dset);
   CAROM_ASSERT(dspace >= 0);

   hsize_t nsel = H5Sget_select_npoints(dspace);
   CAROM_ASSERT(static_cast<int>(nsel) == nelements);

   herr_t errf;
   if (nsel > 0) {
     errf = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
     CAROM_ASSERT(errf >= 0);
   }

   errf = H5Sclose(dspace);
   CAROM_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif
}

bool
HDFDatabase::isInteger(
   const std::string& key)
{
   bool is_int = false;

   if (!key.empty()) {
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      hid_t this_set = H5Dopen(d_group_id, key.c_str());
#endif
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_INT_ARRAY) {
            is_int = true;
         }
         herr_t errf = H5Dclose(this_set);
         CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
         CAROM_NULL_USE(errf);
#endif
      }
   }

   return is_int;
}

bool
HDFDatabase::isDouble(
   const std::string& key)
{
   bool is_double = false;

   if (!key.empty()) {
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      hid_t this_set = H5Dopen(d_group_id, key.c_str());
#endif
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_DOUBLE_ARRAY) {
            is_double = true;
         }
         herr_t errf = H5Dclose(this_set);
         CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
         CAROM_NULL_USE(errf);
#endif
      }
   }

   return is_double;
}

void
HDFDatabase::writeAttribute(
   int type_key,
   hid_t dataset_id)
{
   hid_t attr_id = H5Screate(H5S_SCALAR);
   CAROM_ASSERT(attr_id >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t attr = H5Acreate(dataset_id,
                          "Type",
                          H5T_STD_I8BE,
                          attr_id,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#else
   hid_t attr = H5Acreate(dataset_id,
                          "Type",
                          H5T_STD_I8BE,
                          attr_id,
                          H5P_DEFAULT);
#endif
   CAROM_ASSERT(attr >= 0);

   herr_t errf = H5Awrite(attr, H5T_NATIVE_INT, &type_key);
   CAROM_ASSERT(errf >= 0);

   errf = H5Aclose(attr);
   CAROM_ASSERT(errf >= 0);

   errf = H5Sclose(attr_id);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif
}

int
HDFDatabase::readAttribute(
   hid_t dataset_id)
{
   hid_t attr = H5Aopen_name(dataset_id, "Type");
   CAROM_ASSERT(attr >= 0);

   int type_key;
   herr_t errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
   CAROM_ASSERT(errf >= 0);

   errf = H5Aclose(attr);
   CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
   CAROM_NULL_USE(errf);
#endif

   return type_key;
}

}
