#ifndef included_SnapshotInterpolator_h
#define included_SnapshotInterpolator_h

#include "Interpolator.h"
#include <vector>
#include <string>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Implements PCHIP algorithm.  Based on "A METHOD FOR CONSTRUCTING LOCAL MONOTONE
 * PIECEWISE CUBIC INTERPOLANTS", Fritchs and Butland (1984).  as well as "MONOTONE 
 * PIECEWISE CUBIC INTERPOLATION," Fritsch and Carlson (1980)
 *
 */
class SnapshotInterpolator
{
public:
    SnapshotInterpolator();
    ~SnapshotInterpolator();

    std::vector<Vector*> interpolate(std::vector<Vector*> snapshot_ts, 
                                     std::vector<Vector*> snapshots,
                                     std::vector<Vector*> output_ts);

    std::vector<Vector*> interpolate(std::vector<Vector*> snapshot_ts, 
                                     std::vector<Vector*> snapshots,
                                     int n_out,
                                     std::vector<Vector*>* output_ts);

    std::vector<Vector*> interpolate();

    

private:

    double computeDerivative(double S1, double S2, double h1, double h2);
    double computeH1(double x, double xl, double xr);
    double computeH2(double x, double xl, double xr);
    double computeH3(double x, double xl, double xr);
    double computeH4(double x, double xl, double xr);
    double computePhi(double t);
    double computePsi(double t);
    int sign(double a);
};

}
#endif