/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the nonuniform DMD algorithm.

#include "mfem.hpp"
#include "IncrementalDMD.h"
#include "linalg/Matrix.h"
#include "utils/Utilities.h"

using namespace mfem;

namespace CAROM {

IncrementalDMD::IncrementalDMD(int dim,
			       double dt,
			       Options svd_options,
			       std::string svd_base_file_name,
                               bool alt_output_basis,
                               Vector* state_offset) :
    DMD(dim, dt, alt_output_basis, state_offset)
{
    bg = new BasisGenerator(svd_options,
		    	    true,
		    	    svd_base_file_name);
    svd = new IncrementalSVDBrand(svd_options,
		    		  svd_base_file_name);
}

IncrementalDMD::~IncrementalDMD()
{
    delete bg;
    delete svd;
}

void
IncrementalDMD::load(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    DMD::load(base_file_name);
}

void
IncrementalDMD::save(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());
    CAROM_VERIFY(d_trained);

    DMD::save(base_file_name);
}

Vector*
IncrementalDMD::predict_dt(Vector* u)
{
    IncrementalDMDInternal mats = svd->getAllMatrices();
    Matrix* Up_new;
    Matrix* U_new;
    if (mats.Up->numColumns() == mats.Uq->numRows()) {
	// Liearly dependent sample
	Up_new = mats.Up->mult(mats.Uq);
	U_new = mats.U;
    }
    else {
	// Linearly independent sample
	int r = mats.Up->numColumns();
	Up_new = new Matrix(r+1, r+1, false);
	for (int i = 0; i < r; i++) {
	    for (int j = 0; j < r; j++) {
		Up_new->item(i, j) = mats.Up->item(i, j);
	    }
	}
	Up_new->item(r, r) = 1;
	Up_new = Up_new->mult(mats.Uq);
	U_new = new Matrix(d_dim, r+1, true);
	for (int i = 0; i < d_dim; i++) {
	    for (int j = 0; j < r; j++) {
		U_new->item(i, j) = mats.U->item(i, j);
	    }
	    U_new->item(i, r) = mats.p->item(i);
	}
    }
    Vector* u_proj = Up_new->transposeMult(U_new->transposeMult(u));
    Vector* pred = U_new->mult(Up_new->mult(d_A_tilde->mult(u_proj)));
    
    delete u_proj;
    delete Up_new;
    delete U_new;

    return pred;
}

void
IncrementalDMD::train(int k, const Matrix* W0, double linearity_tol)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    updateDMD(f_snapshots);
    delete f_snapshots;
}

void
IncrementalDMD::updateDMD(const Matrix* f_snapshots)
{
    std::pair<Matrix*, Matrix*> f_snapshot_pair = computeDMDSnapshotPair(f_snapshots);
    Matrix* f_snapshots_in = f_snapshot_pair.first;
    Matrix* f_snapshots_out = f_snapshot_pair.second;
    //std::cout << f_snapshots->numRows() << " x " << f_snapshots->numColumns() << std::endl;
    //Matrix* f_snapshots_in = f_snapshots->getFirstNColumns(f_snapshots->numColumns()-1);
    //Matrix* f_snapshots_out = f_snapshots->getLastNColumns(f_snapshots->numColumns()-1);
    
    /* Incremental SVD
     * 
     * Does everything below all at once:
     * (1) Initialize SVD if not initialized
     * (2) Check linear dependence of new snapshot
     * (3) Rank-1 update
     *
    */
    double* u_in = f_snapshots_in->getColumn(f_snapshots_in->numColumns()-1)
    	    			 ->getData();

    StopWatch svd_timer, rest_timer;

    StopWatch timer1, timer2;

    svd_timer.Start();

    int num_samples_pre = svd->getNumSamples();
    timer2.Start();
    svd->takeSample(u_in, 0, false); // what if norm(u_in) < eps at init?
    delete[] u_in;
    timer2.Stop();
    int num_samples = svd->getNumSamples();
    if (num_samples > num_samples_pre)
    {
	if (num_samples == 1)
	{
	    Matrix* d_basis = new Matrix(*(svd->getSpatialBasis()));
	    Matrix* d_basis_right = new Matrix(*(svd->getTemporalBasis()));
	    Vector* d_sv = new Vector(*(svd->getSingularValues()));
	    Matrix* d_S_inv = new Matrix(1, 1, false);
	    d_S_inv->item(0, 0) = 1 / d_sv->item(0);

	    Matrix* br_Sinv = d_basis_right->mult(d_S_inv);
	    Matrix* f_br_Sinv = f_snapshots_out->mult(br_Sinv);
	    
	    d_A_tilde = d_basis->transposeMult(f_br_Sinv);

	    delete br_Sinv;
	    delete f_br_Sinv;

	    delete d_basis;
	    delete d_basis_right;
	    delete d_sv;
	    delete d_S_inv;
	}
	else
	{
	    if (d_rank == 0) {
	    	std::cout << "Added linearly independent sample" << std::endl;
	    }
	    IncrementalDMDInternal mats = svd->getAllMatrices();
	   
	    Matrix* d_A_tilde_tmp = new Matrix(num_samples, num_samples, false);
	    Vector* u_new = f_snapshots_out->getColumn(f_snapshots_out->numColumns()-1);

	    timer1.Start();
	    Vector* UTu = mats.U->transposeMult(u_new);
	    timer1.Stop();

	    Matrix* f_cropped = f_snapshots_out->getFirstNColumns(
			    			 f_snapshots_out->numColumns()-1);
	    Vector* fTp = f_cropped->transposeMult(mats.p);

	    Vector* UpTUTu = mats.Up->transposeMult(UTu);
	    Vector* WTfTp = mats.W->transposeMult(fTp);
	    for (int i = 0; i < num_samples-1; i++) {
		for (int j = 0; j < num_samples-1; j++) {
		    d_A_tilde_tmp->item(i, j) = d_A_tilde->item(i, j) * mats.s->item(j);
		}
		d_A_tilde_tmp->item(i, num_samples-1) = UpTUTu->item(i);
	    }
	    for (int j = 0; j < num_samples-1; j++) {
		d_A_tilde_tmp->item(num_samples-1, j) = WTfTp->item(j);
	    }
	    d_A_tilde_tmp->item(num_samples-1, num_samples-1) = mats.p->inner_product(u_new);
	   
	    Matrix* WSinv = mats.Wq->mult(mats.Sq_inv);
	    Matrix* AWSinv = d_A_tilde_tmp->mult(WSinv);

	    Matrix* d_A_tilde_new = mats.Uq->transposeMult(AWSinv);
	    
	    delete UTu;
	    delete f_cropped;
	    delete fTp;
	    delete UpTUTu;
	    delete WSinv;
	    delete AWSinv;

	    delete d_A_tilde;
	    d_A_tilde = d_A_tilde_new;

	    delete u_new;
	    delete d_A_tilde_tmp;
	}
    }
    else
    {
	if (d_rank == 0) {
	    std::cout << "Added linearly dependent sample" << std::endl;
	}
	IncrementalDMDInternal mats = svd->getAllMatrices();
	
	Matrix* d_A_tilde_tmp = new Matrix(num_samples, num_samples+1, false);
	Vector* u_new = f_snapshots_out->getColumn(f_snapshots_out->numColumns()-1);

	timer1.Start();
	Vector* UTu = new Vector(num_samples, false);
	mats.U->transposeMult(*u_new, UTu);
	timer1.Stop();

	Vector* UpTUTu = mats.Up->transposeMult(UTu);

	for (int i = 0; i < num_samples; i++) {
	    for (int j = 0; j < num_samples; j++) {
		d_A_tilde_tmp->item(i, j) = d_A_tilde->item(i, j) * mats.s->item(j);
	    }
	    d_A_tilde_tmp->item(i, num_samples) = UpTUTu->item(i);
	}

	Matrix* WSinv = mats.Wq->mult(mats.Sq_inv);
	Matrix* AWSinv = d_A_tilde_tmp->mult(WSinv);

	Matrix* d_A_tilde_new = mats.Uq->transposeMult(AWSinv);
	
	delete UTu;
	delete UpTUTu;
	delete WSinv;
	delete AWSinv;

	delete d_A_tilde;
	d_A_tilde = d_A_tilde_new;

	delete u_new;
	delete d_A_tilde_tmp;
    }

    if (d_rank == 0) {
    	std::cout << "Using " << num_samples << " basis vectors out of "
		  << f_snapshots_out->numColumns() << " snapshots" << std::endl;
    }

    svd_timer.Stop();
    if ( d_rank == 0 ) {
        std::cout << "Time svd:" << svd_timer.RealTime() << std::endl;
        std::cout << "Timer1:" << timer1.RealTime() << std::endl;
        std::cout << "Timer2:" << timer2.RealTime() << std::endl;
    }

    delete f_snapshots_in;
    delete f_snapshots_out;

    d_trained = true;

    return;

    /*
    rest_timer.Start();
    //bg->computeBasis();
    // Copy pointers from IncrementalSVD to DMD
    d_k = svd->getSingularValues()->dim(); // no. of singular values
    
    //d_U = bg->d_U;
    //Matrix* d_Up = bg->getUp();
    
    if (d_trained){ delete d_basis; } 
    Matrix* d_basis = new Matrix(*(svd->getSpatialBasis())); // left singular vectors
    Vector* d_sv = new Vector(*(svd->getSingularValues())); // singular values
    Matrix* d_basis_right = new Matrix(*(svd->getTemporalBasis())); // right singular vectors
    Matrix* d_S_inv = new Matrix(d_k, d_k, false);

    // Get inverse of singular values by multiplying by reciprocal.
    for (int i = 0; i < d_k; ++i)
    {
        d_S_inv->item(i, i) = 1 / d_sv->item(i);
    }

    // Calculate A_tilde = U_transpose * f_snapshots_out * V * inv(S)
    Matrix* d_basis_mult_f_snapshots_out = d_basis->transposeMult(f_snapshots_out);
    Matrix* d_basis_mult_f_snapshots_out_mult_d_basis_right =
        d_basis_mult_f_snapshots_out->mult(d_basis_right);
    Matrix* d_A_tilde_exact = d_basis_mult_f_snapshots_out_mult_d_basis_right->mult(d_S_inv);

    double A_normF = 0.0;
    for (int i = 0; i < d_A_tilde_exact->numRows(); i++) {
	for (int j = 0; j < d_A_tilde_exact->numColumns(); j++) {
	    double diff = d_A_tilde->item(i, j) - d_A_tilde_exact->item(i, j);
	    A_normF += diff * diff;
	}
    }
    std::cout << "Frobenius norm of A_tilde diff: " << sqrt(A_normF) << std::endl;

    delete d_A_tilde;
    d_A_tilde = new Matrix(*d_A_tilde_exact);

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(d_A_tilde);
    d_eigs = eigenpair.eigs;

    struct DMDInternal dmd_internal = {f_snapshots_in, f_snapshots_out, d_basis, d_basis_right, d_S_inv, &eigenpair};
    computePhi(dmd_internal);

    Vector* init = new Vector(f_snapshots_in->numRows(), true);
    for (int i = 0; i < init->dim(); i++)
    {
        init->item(i) = f_snapshots_in->item(i, 0);
    }

    // Calculate pinv(d_phi) * initial_condition.
    projectInitialCondition(init);

    d_trained = true;

    delete d_basis_right;
    delete d_S_inv;
    delete d_basis_mult_f_snapshots_out;
    delete d_basis_mult_f_snapshots_out_mult_d_basis_right;
    delete f_snapshots_in;
    delete f_snapshots_out;
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;
    delete init;

    rest_timer.Stop();
    std::cout << "Time rest:" << rest_timer.RealTime() << std::endl;

    */
}

}
