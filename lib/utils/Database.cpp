/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract database class defines interface for databases.

#include "Database.h"
#include <iostream>
#include <fstream>

namespace CAROM {


Database::Database()
{
}

Database::~Database()
{
}

bool
Database::create(
    const std::string& file_name,
    const MPI_Comm comm)
{
    std::cout << "Creating file: " << file_name << std::endl;
    return true;
}

bool
Database::open(
    const std::string& file_name,
    const std::string& type,
    const MPI_Comm comm)
{
    std::cout << "Opening file: " << file_name << std::endl;
    return true;
}

}
