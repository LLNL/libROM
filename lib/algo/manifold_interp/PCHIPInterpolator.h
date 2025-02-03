#ifndef included_PCHIPInterpolator_h
#define included_PCHIPInterpolator_h

#include <vector>
#include <string>
#include <memory>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Implements PCHIP algorithm.  Based on "A METHOD FOR UCTING LOCAL MONOTONE
 * PIECEWISE CUBIC INTERPOLANTS", Fritsch and Butland (1984).  as well as "MONOTONE
 * PIECEWISE CUBIC INTERPOLATION," Fritsch and Carlson (1980)
 *
 */
class PCHIPInterpolator
{
public:


    PCHIPInterpolator()
    {}
    ~PCHIPInterpolator()
    {}

    /**
     * @brief Compute new snapshots interpolated from snapshot_ts to
     *        output_ts.
     *
     * @param[in] snapshot_ts       The parameter points.
     * @param[in] snapshots         The rotation matrices associated with
     *                              each parameter point.
     * @param[in] output_ts         Requested times for interpolated
     *                              snapshots
     * @param[out] output_snapshots snapshots at output_ts interpolated
     *                              from snapshot_ts
     */
    void interpolate(std::vector<Vector>& snapshot_ts,
                     std::vector<std::shared_ptr<Vector>>& snapshots,
                     std::vector<Vector>& output_ts,
                     std::vector<std::shared_ptr<Vector>>& output_snapshots);

    /**
     * @brief Compute new snapshots interpolated from snapshot_ts to
     *        output_ts.
     *
     * @param[in] snapshot_ts       The parameter points.
     * @param[in] snapshots         The rotation matrices associated with
     *                              each parameter point.
     * @param[in] n_out             Number of output snapshots requested
     * @param[out] output_ts        std::vector of CAROM::Vectors that are
                                    the times of the interpolated snapshots.
     * @param[out] output_snapshots snapshots at output_ts interpolated
     *                              from snapshot_ts
     */
    void interpolate(std::vector<Vector>& snapshot_ts,
                     std::vector<std::shared_ptr<Vector>>& snapshots,
                     int n_out,
                     std::vector<Vector>& output_ts,
                     std::vector<std::shared_ptr<Vector>>& output_snapshots);

private:

    double computeDerivative(double S1, double S2, double h1, double h2) const;
    double computeH1(double x, double xl, double xr) const;
    double computeH2(double x, double xl, double xr) const;
    double computeH3(double x, double xl, double xr) const;
    double computeH4(double x, double xl, double xr) const;
    double computePhi(double t) const;
    double computePsi(double t) const;
    int sign(double a) const;
};

}
#endif
