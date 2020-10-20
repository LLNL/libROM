/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "Matrix.h"
#include "mpi.h"
#include <cmath>

#include <vector>
#include <set>
#include <algorithm>

/* Use Autotools-detected Fortran name-mangling scheme */
#define dgesdd CAROM_FC_GLOBAL(dgesdd, DGESDD)

extern "C" {
void dgesdd(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

namespace CAROM {

  /* Implement QDEIM algorithm from Zlatko Drmac, Serkan Gugercin, "A
     new selection operator for the discrete empirical interpolation
     method -- Improved a priori error bound and extensions", SIAM
     Journal on Scientific Computing, Volume 38, Issue 2, pages
     A631-A648, 2016.
  */
void
QDEIM(const Matrix* f_basis,
      int num_f_basis_vectors_used,
      int* f_sampled_row,
      int* f_sampled_rows_per_proc,
      Matrix& f_basis_sampled_inv,
      const int myid,
      const int num_procs,
      const int num_samples_req)
{
  // This algorithm determines the rows of f that should be sampled, the
  // processor that owns each sampled row, and fills f_basis_sampled_inv with the
  // sampled rows of the basis of the RHS.

  CAROM_VERIFY(num_f_basis_vectors_used == f_basis->numColumns());  // The QR implementation uses the entire matrix.
  CAROM_VERIFY(f_basis->numColumns() <= num_samples_req && num_samples_req <= f_basis->numRows());
  CAROM_VERIFY(num_samples_req == f_basis_sampled_inv.numRows() && f_basis->numColumns() == f_basis_sampled_inv.numColumns());

  // QR will determine (numCol) pivots, which will define the first (numCol) samples.
  const int numCol = f_basis->numColumns();
  const int num_samples_req_QR = numCol;
  
  std::vector<double> sampled_row_data;
  if (myid == 0) sampled_row_data.resize(num_samples_req * numCol);

  int *f_sampled_row_owner = (myid == 0 && f_basis->distributed()) ? new int[num_samples_req] : NULL;

  // QDEIM computes selection/interpolation indices by taking a
  // column-pivoted QR-decomposition of the transpose of its input
  // matrix.
  f_basis->qrcp_pivots_transpose(f_sampled_row,
				 f_sampled_row_owner,
				 num_samples_req_QR);

  if (f_basis->distributed())
    {
      // Gather the sampled rows to root process in f_basis_sampled_inv

      // On root, f_sampled_row contains all global pivots.
      // On non-root processes, it contains the local pivots.
      
      int * ns = (myid == 0) ? new int[num_procs] : NULL;
      int * disp = (myid == 0) ? new int[num_procs] : NULL;
      int * all_sampled_rows = (myid == 0) ? new int[num_samples_req] : NULL;
      if (myid == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    ns[r] = 0;

	  for (int i=0; i<num_samples_req_QR; ++i)
	    ns[f_sampled_row_owner[i]]++;

	  disp[0] = 0;
	  for (int r=1; r<num_procs; ++r)
	    disp[r] = disp[r-1] + ns[r-1];

	  CAROM_VERIFY(disp[num_procs-1] + ns[num_procs-1] == num_samples_req_QR);

	  for (int r=0; r<num_procs; ++r)
	    {
	      f_sampled_rows_per_proc[r] = ns[r];
	      ns[r] = 0;
	    }
	  
	  for (int i=0; i<num_samples_req_QR; ++i)
	    {
	      const int owner = f_sampled_row_owner[i];
	      all_sampled_rows[disp[owner] + ns[owner]] = f_sampled_row[i];
	      ns[owner]++;
	    }

	  // Reorder f_sampled_row and f_sampled_row_owner to match the order of f_basis_sampled_inv
	  for (int i=0; i<num_samples_req_QR; ++i)
	    f_sampled_row[i] = all_sampled_rows[i];

	  int os = 0;
	  for (int r=0; r<num_procs; ++r)
	    {
	      for (int i=0; i<f_sampled_rows_per_proc[r]; ++i)
		{
		  f_sampled_row_owner[os + i] = r;
		}

	      os += f_sampled_rows_per_proc[r];
	    }
	}

      MPI_Bcast(f_sampled_rows_per_proc, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
      
      // TODO: eliminate n and the following MPI_Scatter
      
      int count = 0;
      MPI_Scatter(ns, 1, MPI_INT, &count, 1, MPI_INT, 0, MPI_COMM_WORLD);

      int * my_sampled_rows = (count > 0) ? new int[count] : NULL;
      double * my_sampled_row_data = (count > 0) ? new double[count*num_f_basis_vectors_used] : NULL;  // TODO: change to numCol

      MPI_Scatterv(all_sampled_rows, ns, disp, MPI_INT, my_sampled_rows, count, MPI_INT, 0, MPI_COMM_WORLD);

      int *row_offset = new int[num_procs];
      row_offset[myid] = f_basis->numRows();

      CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, row_offset, 1,
				 MPI_INT, MPI_COMM_WORLD) == MPI_SUCCESS);

      int os = 0;
      for (int i=0; i<num_procs; ++i)
	{
	  os += row_offset[i];
	  row_offset[i] = os - row_offset[i];
	}
      
      for (int i=0; i<count; ++i)
	{
	  const int msr = my_sampled_rows[i];
	  const bool mycond = my_sampled_rows[i] >= row_offset[myid] && my_sampled_rows[i] < row_offset[myid] + f_basis->numRows();
	  CAROM_VERIFY(my_sampled_rows[i] >= row_offset[myid] && my_sampled_rows[i] < row_offset[myid] + f_basis->numRows());
	  const int row = my_sampled_rows[i] - row_offset[myid];
	  os = i*num_f_basis_vectors_used;
	  for (int j=0; j<num_f_basis_vectors_used; ++j)
	    my_sampled_row_data[os + j] = f_basis->item(row, j);
	}

      if (myid == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    {
	      ns[r] *= num_f_basis_vectors_used;
	      disp[r] *= num_f_basis_vectors_used;
	    }
	}

      // TODO: use numCols instead of num_f_basis_vectors_used
      MPI_Gatherv(my_sampled_row_data, count*num_f_basis_vectors_used, MPI_DOUBLE, sampled_row_data.data(), ns, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // At this point, only the first (numCol) entries of f_sampled_row and rows of f_basis_sampled_inv are set by QR.
      // Now set the remaining (num_samples_req - numCol) samples by GappyPOD+E.
      double *A = (myid == 0) ? new double[num_samples_req * numCol] : NULL;

      const int nf = f_basis->numRows();

      std::vector<int> rcnt(num_procs);
      std::vector<int> rdsp(num_procs);
      MPI_Gather(&nf, 1, MPI_INT, rcnt.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

      int nglobal = 0;
      rdsp[0] = 0;
      if (myid == 0)
	{
	  for (int i=0; i<num_procs; ++i)
	    nglobal += rcnt[i];

	  for (int i=1; i<num_procs; ++i)
	    rdsp[i] = rdsp[i-1] + rcnt[i-1];
	}

      std::set<int> globalSamples;
      if (myid == 0)
	{
	  int count = 0;
	  for (int i=0; i<num_procs; ++i)
	    {
	      for (int j=0; j<f_sampled_rows_per_proc[i]; ++j, ++count)
		{
		  //globalSamples.insert(f_sampled_row[count] + row_offset[i]);
		  globalSamples.insert(f_sampled_row[count]);
		}
	    }

	  CAROM_VERIFY(count == numCol);
	}
      
      std::vector<double> rg;
      if (myid == 0) rg.resize(nglobal);

      std::vector<int> isort;
      if (myid == 0) isort.resize(nglobal);

      int n = numCol;
      Matrix V(n, n, false);
      
      for (int s = numCol; s < num_samples_req; ++s)  // Determine sample s
	{
	  int m = s;

	  double g = 0.0;
	  if (myid == 0)
	    {
	      // Compute SVD of the first s rows of f_basis_sampled_inv locally on root
	      CAROM_VERIFY(!f_basis_sampled_inv.distributed());

	      //Matrix Ut(n, m, false);global

	      // Use lapack's dgesdd Fortran function to perform the svd. As this is
	      // Fortran, A and all the computed matrices are in column major order.

	      for (int i=0; i<m; ++i)
		for (int j=0; j<n; ++j)
		  {
		    //A[i + (j*m)] = f_basis_sampled_inv.item(i, j);
		    A[i + (j*m)] = sampled_row_data[(i*numCol) + j];
		    //Ut.item(j, i) = f_basis_sampled_inv.item(i, j);
		  }

	      std::vector<double> sigma(n);
	      char jobz = 'O';
	      int lda = m;
	      int ldu = 1;
	      int ldv = n;
	      int lwork = (3*n) + std::max(m, (5*n*n) + (4*n));
	      double* work = new double [lwork];
	      int iwork[8*n];
	      int info;
	      dgesdd(&jobz, &m, &n, A, &lda, sigma.data(), NULL, &ldu, &V.item(0, 0),
		     &ldv, work, &lwork, iwork, &info);

	      CAROM_VERIFY(info == 0);

	      delete [] work;

	      g = (sigma[n-2] * sigma[n-2]) - (sigma[n-1] * sigma[n-1]);

	      // Note that V stores the right singular vectors row-wise

	      // Set Ub = V' * U'
	      //Matrix *Ub = V.mult(Ut);  // size n x m

	      // Set Ubt = U * V
	      V.transpose();
	    }

	  MPI_Bcast(&g, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  // Broadcast the small n-by-n undistributed matrix V which is computed only on root.
	  MPI_Bcast(V.getData(), n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  Matrix *Ubt = f_basis->mult(V);  // distributed

	  CAROM_VERIFY(Ubt->distributed() && Ubt->numRows() == f_basis->numRows() && Ubt->numColumns() == n);

	  std::vector<double> r(nf);

	  for (int i=0; i<nf; ++i)
	    {
	      r[i] = g;
	      for (int j=0; j<n; ++j)
		r[i] += Ubt->item(i, j) * Ubt->item(i, j);  // column sums of Ub.^2, which are row sums of Ubt.^2

	      r[i] -= sqrt((r[i] * r[i]) - (4.0 * g * Ubt->item(i, n-1) * Ubt->item(i, n-1)));
	    }

	  MPI_Gatherv(r.data(), nf, MPI_DOUBLE, rg.data(), rcnt.data(), rdsp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  int owner = -1;
	  if (myid == 0)
	    {
	      for (int i=0; i<nglobal; ++i)
		{
		  isort[i] = i;
		}

	      //std::sort(isort.begin(), isort.end(), std::greater<double>());  // descending order
	      std::sort(isort.begin(), isort.end(), [&](const int& a, const int& b) {
		  return (rg[a] > rg[b]); }
		);  // descending order

	      // Choose sample s as the first entry in isort not already in f_sampled_row
	      f_sampled_row[s] = -1;
	      for (int i=0; i<nglobal; ++i)
		{
		  std::set<int>::iterator it = globalSamples.find(isort[i]);
		  if (it == globalSamples.end()) // not found
		    {
		      CAROM_VERIFY(i < s-1 && f_sampled_row[s] == -1);
		      f_sampled_row[s] = isort[i];
		      break;
		    }
		}

	      CAROM_VERIFY(f_sampled_row[s] >= 0);
	      for (int i=num_procs-1; i >= 0; --i)
		{
		  if (f_sampled_row[s] >= row_offset[i])
		    {
		      owner = i;
		      break;
		    }
		}

	      CAROM_VERIFY(owner >= 0);
	      f_sampled_rows_per_proc[owner]++;
	      f_sampled_row_owner[s] = owner;

	      for (int i=0; i < num_procs; ++i)
		ns[i] = -1;

	      ns[owner] = f_sampled_row[s];
	      globalSamples.insert(f_sampled_row[s]);
	    }

	  // Send one row of f_basis, corresponding to sample s, to root process for f_basis_sampled_inv
	  // First, scatter from root to tell the owning process the sample index.

	  int sample = -1;
	  MPI_Scatter(ns, 1, MPI_INT, &sample, 1, MPI_INT, 0, MPI_COMM_WORLD);

	  const int tagSendRecv = 111;
	  if (sample > -1)
	    {
	      CAROM_VERIFY(sample >= row_offset[myid]);
	      MPI_Send(f_basis->getData() + ((sample - row_offset[myid]) * numCol), numCol, MPI_DOUBLE, 0, tagSendRecv, MPI_COMM_WORLD);
	    }

	  if (myid == 0)
	    {
	      MPI_Status status;
	      MPI_Recv(sampled_row_data.data() + (s*numCol), numCol, MPI_DOUBLE, owner, tagSendRecv, MPI_COMM_WORLD, &status);
	    }

	  delete Ubt;
	}  // loop s over samples

      // Subtract row_offset to convert f_sampled_row from global to local indices
      // Also, reorder f_sampled_row by process
      /*
      os = row_offset[myid];
      for (int i=0; i<num_f_basis_vectors_used; ++i)
	{
	  f_sampled_row[i] -= os;
	}
      */

      if (myid == 0)
	{
	  ns[0] = 0;
	  disp[0] = 0;
	  for (int r=1; r<num_procs; ++r)
	    {
	      ns[r] = 0;
	      disp[r] = disp[r-1] + f_sampled_rows_per_proc[r-1];
	    }

	  for (int i=0; i<num_samples_req; ++i)
	    {
	      const int owner = f_sampled_row_owner[i];
	      f_sampled_row[i] -= row_offset[owner];

	      all_sampled_rows[disp[owner] + ns[owner]] = i; //f_sampled_row[i];
	      ns[owner]++;
	    }

	  std::vector<int> sortedRow(num_samples_req);

	  for (int r=0; r<num_procs; ++r)
	    {
	      CAROM_VERIFY(ns[r] == f_sampled_rows_per_proc[r]);

	      // Sort the local indices of f_sampled_row for process r.
	      for (int i=0; i<ns[r]; ++i)
		{
		  isort[i] = i;
		}

	      std::sort(isort.begin(), isort.begin() + ns[r], [&](const int& a, const int& b) {
		  //return (f_sampled_row[disp[r] + a] < f_sampled_row[disp[r] + b]); }
		  return (f_sampled_row[all_sampled_rows[disp[r] + a]] < f_sampled_row[all_sampled_rows[disp[r] + b]]); }
		);  // ascending order locally

	      for (int i=0; i<ns[r]; ++i)
		{
		  const int ig = all_sampled_rows[disp[r] + isort[i]];
		  const int s = (disp[r] + i);
		  // Put row ig of sampled_row_data into row s of f_basis_sampled_inv
		  // Put entry ig of f_sampled_row into entry s of sortedRow

		  for (int j=0; j<numCol; ++j)
		    {
		      f_basis_sampled_inv.item(s, j) = sampled_row_data[(ig*numCol) + j];
		    }

		  sortedRow[s] = f_sampled_row[ig];
		}
	    }

	  for (int i=0; i<num_samples_req; ++i)
	    f_sampled_row[i] = sortedRow[i]; // all_sampled_rows[i];
	}

      MPI_Bcast(f_sampled_row, num_samples_req, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(f_sampled_rows_per_proc, num_procs, MPI_INT, 0, MPI_COMM_WORLD);

      delete [] A;
      delete [] row_offset;
      delete [] ns;
      delete [] disp;
      delete [] my_sampled_rows;
      delete [] all_sampled_rows;
    }
  else
    {
      // With the known interpolation (sample) indices, copy over the
      // rows of the sampled basis
      for (int i = 0; i < num_samples_req_QR; i++) {
	for (int j = 0; j < num_f_basis_vectors_used; j++) {
	  f_basis_sampled_inv.item(i, j) = f_basis->item(f_sampled_row[i], j);
	}
      }

      f_sampled_rows_per_proc[0] = num_f_basis_vectors_used;
    }

  delete [] f_sampled_row_owner;

  if (!f_basis->distributed())
    {
      CAROM_VERIFY(false);  // TODO: implement GappyPOD+E for this case
    }
  
  // Now invert f_basis_sampled_inv, storing its transpose.
  if (myid == 0)  // Matrix is valid only on root process
    f_basis_sampled_inv.transposePseudoinverse();
} // end void QDEIM

} // end namespace CAROM
