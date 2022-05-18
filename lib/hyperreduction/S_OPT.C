/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the S_OPT algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "linalg/Matrix.h"
#include "linalg/scalapack_wrapper.h"
#include "mpi.h"
#include <cmath>
#include <vector>
#include <map>
#include <set>

#include "S_OPT.h"

using namespace std;

namespace CAROM {

void
S_OPT(const Matrix* f_basis,
      bool use_qr_basis,
      int num_f_basis_vectors_used,
      std::vector<int>& f_sampled_row,
      std::vector<int>& f_sampled_rows_per_proc,
      const int myid,
      const int num_procs,
      const int num_samples_req)
{
    int num_rows = f_basis->numRows();
    int num_cols = f_basis->numColumns();

    Matrix* Vo = NULL;

    if (use_qr_basis)
    {
        // Get QR factorization of basis
        int *row_offset = new int[num_procs + 1];
        row_offset[num_procs] = f_basis->numDistributedRows();
        row_offset[myid] = num_rows;

        CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
                                   1,
                                   MPI_INT,
                                   row_offset,
                                   1,
                                   MPI_INT,
                                   MPI_COMM_WORLD) == MPI_SUCCESS);
        for (int i = num_procs - 1; i >= 0; i--) {
            row_offset[i] = row_offset[i + 1] - row_offset[i];
        }

        std::cout << row_offset[0] << " " << row_offset[1] << std::endl;

        CAROM_VERIFY(row_offset[0] == 0);

        SLPK_Matrix slpk_f_basis;

        int nrow = num_procs;
        int ncol = 1;
        int blocksize = row_offset[num_procs] / num_procs;
        if (row_offset[num_procs] % num_procs != 0) blocksize += 1;
        initialize_matrix(&slpk_f_basis, f_basis->numDistributedRows(), num_cols,
                          nrow, ncol, blocksize, num_cols);
        for (int rank = 0; rank < num_procs; ++rank) {
            scatter_block(&slpk_f_basis, 1, row_offset[rank] + 1,
                          f_basis->getData(), row_offset[rank + 1] - row_offset[rank],
                          num_cols, rank);
        }
        delete f_basis;

        print_debug_info(&slpk_f_basis);

        QRManager QRmgr;
        qr_init(&QRmgr, &slpk_f_basis);
        qrfactorize(&QRmgr);

        print_debug_info(QRmgr.A);

        // Obtain Q
        qrcompute(&QRmgr);
        print_debug_info(QRmgr.A);
        Vo = new Matrix(row_offset[myid + 1] - row_offset[myid], num_cols, true);
        std::cout << Vo->numRows() << " " << Vo->numColumns() << std::endl;
        for (int rank = 0; rank < num_procs; ++rank) {
            gather_transposed_block(&Vo->item(0, 0), QRmgr.A,
                         row_offset[rank] + 1, 1,
                         row_offset[rank + 1] - row_offset[rank], num_cols,
                         rank);
        }
        // for (int i = 0; i < num_rows; i++)
        // {
        //     for (int j = 0; j < num_cols; j++)
        //     {
        //         std::cout << Vo->item(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }
        free(QRmgr.tau);
        free(QRmgr.ipiv);
        free_matrix_data(QRmgr.A);
    }
    else
    {
        Vo = new Matrix(*f_basis);
    }


    // Square Vo.
    Matrix* nVo = Vo->pointwise_square();

    std::vector<double> indices;
    double max_val = 0.0;
    int index = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if (std::abs(Vo->item(i, 0)) > max_val)
        {
            max_val = std::abs(Vo->item(i, 0));
            index = i;
        }
    }
    // Do A GATHER HERE TO FIND TRUE MAX
    indices.push_back(index);
    std::cout << indices[0] << std::endl;

    // if  inum == 0
    //     [~, i_idx] = max(abs(Vo(:,1))); % Get largest element in the absolute valued version of the first row, save its index. What if there are multiple? I guess just take the first.
    //     index = [];
    //     index = [index; i_idx]; % Save index
    //     inum = inum + 1;
    // end

    for (int i = 2; i < num_rows; i++)
    {
        if (i <= num_cols)
        {
            Matrix* A0 = new Matrix(i - 1, i - 1, true);
            Matrix* V1_last_col = new Matrix(i - 1, 1, true);
            for (int j = 0; j < indices.size(); j++)
            {
                for (int k = 0; k < i - 1; k++)
                {
                    A0->item(j, k) = Vo->item(indices[j], k);
                    // std::cout << Vo->item(indices[j], k) << std::endl;
                    // std::cout << "A" << std::endl;
                }
                V1_last_col->item(j, 0) = Vo->item(indices[j], i - 1);
                // std::cout << Vo->item(indices[j], i-2) << std::endl;
                // std::cout << Vo->item(indices[j], i-1) << std::endl;
                // std::cout << Vo->item(indices[j], i) << std::endl;
                // std::cout << Vo->item(indices[j], i+1) << std::endl;
                // std::cout << Vo->item(indices[j], i+2) << std::endl;
                // std::cout << "B" << std::endl;
            }
            // for (int j = 0; j < A0->numRows(); j++)
            // {
            //     for (int k = 0; k < A0->numColumns(); k++)
            //     {
            //         std::cout << A0->item(j, k) << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;
            // for (int j = 0; j < V1_last_col->numRows(); j++)
            // {
            //     for (int k = 0; k < V1_last_col->numColumns(); k++)
            //     {
            //         std::cout << V1_last_col->item(j, k) << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // abort();

            Matrix* atA0 = V1_last_col->transposeMult(A0);
            Matrix* tt = new Matrix(num_rows, i - 1, true);
            Matrix* tt1 = new Matrix(i - 1, num_rows, true);
            Matrix* V1_last_col_squared = V1_last_col->pointwise_square();

            double ata = 0.0;
            for (int j = 0; j < V1_last_col_squared->numRows(); j++)
            {
                for (int k = 0; k < V1_last_col_squared->numColumns(); k++)
                {
                    ata += V1_last_col_squared->item(j, k);
                }
            }

            Matrix* A0_T_mult_A0 = A0->transposeMult(A0);
            Matrix* rhs = new Matrix(i - 1, num_rows + (atA0->numRows() * atA0->numColumns()), true);
            for (int k = 0; k < rhs->numRows(); k++)
            {
                rhs->item(k, 0) = atA0->item(0, k);
            }
            for (int j = 1; j < rhs->numColumns(); j++)
            {
                for (int k = 0; k < rhs->numRows(); k++)
                {
                    rhs->item(k, j) = Vo->item(j - 1, k);
                }
            }



            // V1 = Vo(index, 1:i); % For all the saved rows, get their columns from 1 to i
            // A0 = V1(1:i-1,1:i-1); % Get upper left rectangle of the previous matrix that only has the saved rows and columns from 1 to i
            // atA0 = V1(1:i-1,i)'*A0; % Do this transpose and multiplication, store in atA0
            // tt = zeros(N_Bm,i-1);
            // tt1 = zeros(i-1,N_Bm);
            // ata = sum(V1(1:i-1,i).^2) ; % Square every term in the matrix then sum then up to a scalar
            // bbb = (A0'*A0)\[atA0; Vo(:,1:i-1)]'; % solve linear equation left * bbb = right
            // c = bbb(:,2:end);
            // b = 1 + sum(Vo(:,1:i-1).*c',2); % Do a sum, on axis = 2
            //
            // for zz = 1 : i-1
            //     tt(:,zz) = Vo(:,zz).*Vo(:,i);
            // end
            // g1 = repmat(atA0,N_Bm,1) + tt; % create array copied from fill value
            // g1 = g1';
            // g2 = bbb(:,1);
            // oneprc = 1 + sum(Vo(:,1:i-1).*c',2);
            // g3 = sum(c'.*g1',2)./oneprc;
            // for zz = 1 : i-1
            //     tt1(zz,:) = c(zz,:).*(Vo(:,i) - g3)';
            // end
            //
            // GG = repmat(g2,1,N_Bm) + tt1;
            //
            // A = ata + Vo(:,i).^2 - sum(g1'.*GG',2);
            // A(A<0) = 0; % set all negative values to zero
            // nV = sum(nVo(index,1:i),1);
            // noM = sum(log(repmat(nV,N_Bm,1) + nVo(:,1:i)),2);
            //
            // A = log(abs(A)) + log(b) - noM; % do natural logarithm (ln) of every value in matrix
        }
        else
        {

        }
    }

    delete Vo;
}

}
