/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Utility functions for error reporting, etc.

#include "Utilities.h"
#include "ParallelBuffer.h"

#include "mpi.h"

#include <iomanip>
#include <stdlib.h>
#include <sys/stat.h>

namespace CAROM {

void
Utilities::abort(
    const std::string& message,
    const std::string& filename,
    int line)
{
    ParallelBuffer perr_buffer;
    std::ostream perr(&perr_buffer);
    perr << "Program abort called in file ``" << filename
         << "'' at line " << line << std::endl;
    perr << "ERROR MESSAGE: " << std::endl << message.c_str() << std::endl;
    perr << std::flush;

    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if (size > 1) {
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else {
            exit(-1);
        }
    }
    else {
        exit(-1);
    }
}

std::string
Utilities::processorToString(
    int processorID)
{
    std::ostringstream os;
    if (processorID < 0) {
        os << '-' << std::setw(6) << std::setfill('0') << -processorID;
    } else {
        os << std::setw(7) << std::setfill('0') << processorID;
    }
    os << std::flush;

    return os.str();
}

bool
Utilities::file_exist(
    const std::string& filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

}
