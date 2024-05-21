#ifndef included_SnapshotDMD_h
#define included_SnapshotDMD_h

#include "DMD.h"
#include <vector>

namespace CAROM {

class Vector;
class SnapshotDMD : public DMD
{
public:
    SnapshotDMD(int dim, double dt, bool alt_output_basis = false,
                Vector* state_offset = NULL) : DMD(dim,dt,alt_output_basis,state_offset) {}
    SnapshotDMD(std::string base_file_name) : DMD(base_file_name) {}
    ~SnapshotDMD();
    void interpolateToNSnapshots(int n, int run, int window);
    void interpolateToNSnapshots(int n);
    void train(int k, const Matrix* W0 = NULL, double linearity_tol = 0.0);
    void train(double energy_fraction, const Matrix* W0 = NULL,
               double linearity_tol = 0.0);
    std::vector<Vector*> getSnapshotVectors()
    {
        std::vector<Vector*> return_snapshots(d_snapshots);
        return return_snapshots;
    }
    double checkProjectionError(Vector* init);
protected:
    friend void getParametricDMD<SnapshotDMD>(SnapshotDMD*& parametric_dmd,
            std::vector<Vector*>& parameter_points,
            std::vector<SnapshotDMD*>& dmds,
            Vector* desired_point,
            std::string rbf,
            std::string interp_method,
            double closest_rbf_val,
            bool reorthogonalize_W);

    SnapshotDMD(std::vector<std::complex<double>> eigs, Matrix* phi_real,
                Matrix* phi_imaginary, int k,
                double dt, double t_offset,
                Vector* state_offset) :
        DMD(eigs, phi_real, phi_imaginary, k, dt, t_offset, state_offset) {}


private:
};
}
#endif