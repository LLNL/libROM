/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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

bool fileExists(const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
    // ifstream f will be closed upon the end of the function.
}

Database::Database()
{
}

Database::~Database()
{
}

bool
Database::create(const std::string& file_name)
{
    std::cout << "Creating file: " << file_name << std::endl;
}

bool
Database::open(
    const std::string& file_name,
    const std::string& type)
{
    std::cout << "Opening file: " << file_name << std::endl;
}

}
