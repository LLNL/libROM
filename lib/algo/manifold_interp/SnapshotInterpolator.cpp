/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

/**
 * Implements PCHIP algorithm.  Based on "A METHOD FOR CONSTRUCTING LOCAL MONOTONE
 * PIECEWISE CUBIC INTERPOLANTS", Fritsch and Butland (1984).  as well as "MONOTONE
 * PIECEWISE CUBIC INTERPOLATION," Fritsch and Carlson (1980)
 *
 */
#include "SnapshotInterpolator.h"

#include <cfloat>
#include <limits.h>
#include <cmath>
#include "linalg/Vector.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif
using namespace std;

namespace CAROM {

void SnapshotInterpolator::interpolate(std::vector<Vector*>& snapshot_ts,
        std::vector<Vector*>& snapshots,
        std::vector<Vector*>& output_ts,
        std::vector<Vector*>& output_snapshots)
{
    CAROM_VERIFY(snapshot_ts.size() == snapshots.size());
    CAROM_VERIFY(snapshot_ts.size() > 0);
    CAROM_VERIFY(output_ts.size() > 0);
    CAROM_VERIFY(output_snapshots.size() == 0);

    int n_out = output_ts.size();
    int n_snap = snapshots.size();
    int n_dim = snapshots[0]->dim();

    for(int i = 0; i < n_out; ++i)
    {
        Vector* temp_snapshot = new Vector(snapshots[0]->dim(),
                                          snapshots[0]->distributed());
        output_snapshots.push_back(temp_snapshot);
    }

    for(int i = 0; i < n_dim; ++i)
    {
        double h_temp,delta_temp, t;
        std::vector<double> h,d,delta,t_in,y_in;
        for(int j = 0; j < n_snap-1; ++j)
        {
            h_temp = snapshot_ts[j+1]->getData()[0] - snapshot_ts[j]->getData()[0];
            h.push_back(h_temp);
            delta_temp = (snapshots[j+1]->getData()[i] - snapshots[j]->getData()[i])/h_temp;
            delta.push_back(delta_temp);
            t_in.push_back(snapshot_ts[j]->getData()[0]);
            y_in.push_back(snapshots[j]->getData()[i]);
        }
        t_in.push_back(snapshot_ts[n_snap-1]->getData()[0]);
        y_in.push_back(snapshots[n_snap-1]->getData()[i]);

        double d_temp = ((2*h[0] + h[1])*delta[0] - h[0]*delta[1])/(h[0]+h[1]);
        if(sign(d_temp)!=sign(delta[0]))
        {
            d_temp = 0;
        }
        else if (sign(delta[0]) != sign(delta[1]) && abs(d_temp) > abs(3*delta[0]))
        {
            d_temp = 3*delta[0];
        }
        d.push_back(d_temp);
        int counter = 0;
        for(int j = 0; j < n_snap-2; ++j)
        {
            d_temp = computeDerivative(delta[j],delta[j+1],h[j],h[j+1]);
            d.push_back(d_temp);
            while(output_ts[counter]->getData()[0] <= t_in[j+1])
            {
                t = output_ts[counter]->getData()[0];
                output_snapshots[counter]->getData()[i] = y_in[j]*computeH1(t,t_in[j],
                        t_in[j+1]) +
                        y_in[j+1]*computeH2(t,t_in[j],t_in[j+1]) +
                        d[j]*computeH3(t,t_in[j],t_in[j+1]) +
                        d[j+1]*computeH4(t,t_in[j],t_in[j+1]);
                counter += 1;
            }
        }

        d_temp = ((2*h[n_snap-2]+h[n_snap-3])*delta[n_snap-2] - h[n_snap-2]*delta[n_snap
                  -3])/(h[n_snap-2]+h[n_snap-3]);
        if(sign(d_temp) != sign(delta[n_snap-2]))
        {
            d_temp = 0;
        }
        else if (sign(delta[n_snap-2]) != sign(delta[n_snap-3])
                 && abs(d_temp) > abs(3*delta[n_snap-2]))
        {
            d_temp = 3*delta[n_snap-2];
        }
        d.push_back(d_temp);
        
        while(counter < n_out && output_ts[counter]->getData()[0] - t_in[n_snap-1] <= FLT_EPSILON )
        {
            t = output_ts[counter]->getData()[0];
            output_snapshots[counter]->getData()[i] = y_in[n_snap-2]*computeH1(t,
                    t_in[n_snap-2],t_in[n_snap-1]) +
                    y_in[n_snap-1]*computeH2(t,t_in[n_snap-2],t_in[n_snap-1]) +
                    d[n_snap-2]*computeH3(t,t_in[n_snap-2],t_in[n_snap-1]) +
                    d[n_snap-1]*computeH4(t,t_in[n_snap-2],t_in[n_snap-1]);
            counter += 1;
        }
    }
}

void SnapshotInterpolator::interpolate(std::vector<Vector*>&
        snapshot_ts,
        std::vector<Vector*>& snapshots,
        int n_out,
        std::vector<Vector*>& output_ts,
        std::vector<Vector*>& output_snapshots)
{
    CAROM_VERIFY(snapshot_ts.size() == snapshots.size());
    CAROM_VERIFY(snapshot_ts.size() > 0);
    CAROM_VERIFY(n_out > 0);
    CAROM_VERIFY(output_ts.size() == 0 && output_snapshots.size() == 0);


    int n_snap = snapshots.size();
    int n_dim = snapshots[0]->dim();



    double t_min = snapshot_ts[0]->getData()[0];
    double t_max = snapshot_ts[n_snap-1]->getData()[0];
    double dt = (t_max-t_min)/(n_out-1);
    output_ts.clear();

    for(int i = 0; i < n_out; ++i)
    {
        Vector* temp_t = new Vector(1, false);
        temp_t->getData()[0] = t_min + i*dt;
        output_ts.push_back(temp_t);
    }
    std::cout << "Interpolating to " << n_out << " snapshots, from " << n_snap <<
              std::endl;
    interpolate(snapshot_ts,snapshots,output_ts, output_snapshots);
}

const double SnapshotInterpolator::computeDerivative(double S1, double S2, double h1,
        double h2)
{
    double d = 0.0;
    double alpha = (h1 + 2*h2)/(3*(h1+h2));
    if(S1*S2 > 0)
    {
        d = S1*S2/(alpha*S2 + (1-alpha)*S1);
    }
    else
    {
        d = 0.0;
    }
    return d;
}

const double SnapshotInterpolator::computeH1(double x, double xl, double xr)
{
    double h = xr - xl;
    return computePhi((xr-x)/h);
}

const double SnapshotInterpolator::computeH2(double x, double xl, double xr)
{
    double h = xr - xl;
    return computePhi((x-xl)/h);
}

const double SnapshotInterpolator::computeH3(double x, double xl, double xr)
{
    double h = xr-xl;
    return -h*computePsi((xr-x)/h);
}

const double SnapshotInterpolator::computeH4(double x, double xl, double xr)
{
    double h = xr-xl;
    return h*computePsi((x-xl)/h);
}

const double SnapshotInterpolator::computePhi(double t)
{
    return 3.*pow(t,2.) - 2*pow(t,3.);
}

const double SnapshotInterpolator::computePsi(double t)
{
    return pow(t,3.) - pow(t,2.);
}

const int SnapshotInterpolator::sign(double a)
{
    double TOL = 1e-15;
    if(abs(a) < TOL)return 0;
    else if(a > 0) return 1;
    else if (a < 0) return -1;
    return 0;
}

}