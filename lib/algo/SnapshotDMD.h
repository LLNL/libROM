#ifndef included_SnapshotDMD_h
#define included_SnapshotDMD_h

#include "DMD.h"
#include <vector>

namespace CAROM {

class Vector;
class SnapshotDMD : public DMD
{
public:

    /**
     * @brief Constructor. Basic DMD with uniform time step size.
     *        Inherited directly from base DMD class.
     *
     * @param[in] dim              The full-order state dimension.
     * @param[in] dt               The dt between samples.
     * @param[in] alt_output_basis Whether to use the alternative basis for
     *                             output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset     The state offset.
     */
    SnapshotDMD(int dim, double dt, bool alt_output_basis = false,
                std::shared_ptr<Vector> state_offset = nullptr) :
        DMD(dim, dt, alt_output_basis, state_offset) {}

    /**
     * @brief Constructor. DMD from saved models. Inherited directly
     *        from base DMD class.
     *
     * @param[in] base_file_name The base part of the filename of the
     *                           database to load when restarting from a save.
     */
    SnapshotDMD(std::string base_file_name) : DMD(base_file_name) {}

    /**
    * @brief Destroy the SnapshotDMD object
    */
    ~SnapshotDMD();

    /**
    * @brief Interpolate the current snapshots to n, new snapshots
    *        distributed uniformly over the currently sampled time
    *        domain.
    * @param[in] n                  The number of desired snapshots.
    */
    void interpolateToNSnapshots(int n);

    /**
     * @brief Train the DMD model with specified reduced dimension. If k is
     *        too large then new snapshots are computed using
     *        interpolateToNSnapshots(k+1).
     *
     * @param[in] k               The number of modes to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
     *                            of W is linearly independent with W0.
     */
    void train(int k, const Matrix* W0 = NULL, double linearity_tol = 0.0);

    /**
     * @brief Train the DMD model with energy fraction criterion.
     *
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
     *                            of W is linearly independent with W0.
     */
    void train(double energy_fraction, const Matrix* W0 = NULL,
               double linearity_tol = 0.0);

    /**
     * @brief Returns a copy of the current snapshot vector "d_snapshots"
     */
    std::vector<std::shared_ptr<Vector>> getSnapshotVectors()
    {
        return d_snapshots;
    }
protected:
    /**
     * @brief Obtain DMD model interpolant at desired parameter point by
     *        interpolation of DMD models from training parameter points.
     *
     * @param[in] parametric_dmd    The interpolant DMD model at the desired point.
     * @param[in] parameter_points  The training parameter points.
     * @param[in] dmds              The DMD objects associated with
     *                              each training parameter point.
     * @param[in] desired_point     The desired point at which to create a parametric DMD.
     * @param[in] rbf               The RBF type ("G" == gaussian,
     *                              "IQ" == inverse quadratic,
     *                              "IMQ" == inverse multiquadric)
     * @param[in] interp_method     The interpolation method type
     *                              ("LS" == linear solve,
     *                              "IDW" == inverse distance weighting,
     *                              "LP" == lagrangian polynomials)
     * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
     *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
     * @param[in] reorthogonalize_W Whether to reorthogonalize the interpolated W (basis) matrix.
     */
    friend void getParametricDMD<SnapshotDMD>(std::unique_ptr<SnapshotDMD>&
            parametric_dmd,
            const std::vector<Vector>& parameter_points,
            std::vector<SnapshotDMD*>& dmds,
            const Vector & desired_point,
            std::string rbf,
            std::string interp_method,
            double closest_rbf_val,
            bool reorthogonalize_W);

    /**
     * @brief Constructor.
     *
     * @param[in] eigs d_eigs
     * @param[in] phi_real d_phi_real
     * @param[in] phi_imaginary d_phi_imaginary
     * @param[in] k d_k
     * @param[in] dt d_dt
     * @param[in] t_offset d_t_offset
     * @param[in] state_offset d_state_offset
     * @param[in] derivative_offset d_derivative_offset
     */
    SnapshotDMD(std::vector<std::complex<double>> & eigs,
                std::shared_ptr<Matrix> & phi_real,
                std::shared_ptr<Matrix> & phi_imaginary, int k, double dt,
                double t_offset, std::shared_ptr<Vector> & state_offset) :
        DMD(eigs, phi_real, phi_imaginary, k, dt, t_offset, state_offset) {}

private:
};
}
#endif
