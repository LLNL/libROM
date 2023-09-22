// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// 3D flow over a cylinder benchmark example

#include "navier_solver.hpp"
#include <fstream>

using namespace mfem;
using namespace navier;

struct s_NavierContext
{
   int order = 3;
   double kin_vis = 0.01;
   double t_final = 10.0;
   double dt = 1e-2;
} ctx;

void vel(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);

   double U = 1;
   double L = 12;

   if (xi <= 1e-8 || yi <= 1e-8 || (L-yi) <= 1e-8)
   {
      //u(0) = 4.0*U*yi*(L-yi)/pow(L,2.0);
      u(0) = U;
   }
   else
   {
      u(0) = 0.0;
   }
   u(1) = 0.0;
}

double pres(const Vector &x, double t) {
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   int serial_refinements = 1;

   Mesh *mesh = new Mesh("square-circle.msh");

   for (int i = 0; i < serial_refinements; ++i)
   {
      mesh->UniformRefinement();
   }

   if (Mpi::Root())
   {
      std::cout << "Number of elements: " << mesh->GetNE() << std::endl;
   }

   auto *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   // Create the flow solver.
   NavierSolver flowsolver(pmesh, ctx.order, ctx.kin_vis);
   flowsolver.EnablePA(false); // Gmsh does not support partial assembly?

   // Set the initial condition.
   ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
   VectorFunctionCoefficient u_excoeff(pmesh->Dimension(), vel);
   u_ic->ProjectCoefficient(u_excoeff);

   bool read_vec = 1;
   bool write_vec = 0;

   ParGridFunction *u_gf = flowsolver.GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver.GetCurrentPressure();

   std::cout << "Rank: " << pmesh->GetMyRank()
	     << ", DOF(u): " << u_gf->Size()
	     << ", DOF(p): " << p_gf->Size() << std::endl;
   if (read_vec) {
      char fNameIn_u[7];
      char fNameIn_p[7];
      int rank = pmesh->GetMyRank();
      sprintf(fNameIn_u, "u_gf.%d", rank);
      sprintf(fNameIn_p, "p_gf.%d", rank);
      std::ifstream in(fNameIn_u);
      ParGridFunction u_gf_in(pmesh, in);
      in.close();
      in.open(fNameIn_p);
      ParGridFunction p_gf_in(pmesh, in);
      in.close();
      *u_gf = u_gf_in;
      *p_gf = p_gf_in;
   }

   // Add Dirichlet boundary conditions to velocity space restricted to
   // selected attributes on the mesh.
   Array<int> attr(pmesh->bdr_attributes.Max());
   Array<int> pattr(pmesh->bdr_attributes.Max());
   // Outlet is attribute 2.
   for (int i = 0; i < pmesh->bdr_attributes.Max(); i++) {
      attr[i] = 1;
   }
   attr[1] = 0;
   flowsolver.AddVelDirichletBC(vel, attr);
   pattr[1] = 1;
   flowsolver.AddPresDirichletBC(pres, pattr);

   double t = 0.0;
   double dt = ctx.dt;
   double t_final = ctx.t_final;
   bool last_step = false;
   double cfl;

   flowsolver.Setup(dt);

   ParGridFunction w_gf(*u_gf);
   flowsolver.ComputeCurl2D(*u_gf, w_gf);

   ParaViewDataCollection pvdc("2dfoc_fom", pmesh);
   pvdc.SetDataFormat(VTKFormat::BINARY32);
   pvdc.SetHighOrderOutput(true);
   pvdc.SetLevelsOfDetail(ctx.order);
   pvdc.SetCycle(0);
   pvdc.SetTime(t);
   pvdc.RegisterField("velocity", u_gf);
   pvdc.RegisterField("pressure", p_gf);
   pvdc.RegisterField("vorticity", &w_gf);
   pvdc.Save();

   StopWatch timer;

   for (int step = 0; !last_step; ++step)
   {
      timer.Clear();
      timer.Start();

      if (t + dt >= t_final - dt / 2)
      {
         last_step = true;
      }

      flowsolver.Step(t, dt, step);

      if (step % 10 == 0)
      {
	 flowsolver.ComputeCurl2D(*u_gf, w_gf);
         pvdc.SetCycle(step);
         pvdc.SetTime(t);
         pvdc.Save();
      }

      cfl = flowsolver.ComputeCFL(*u_gf, dt);

      if (Mpi::Root())
      {
         printf("%11s %11s %11s\n", "Time", "dt", "CFL");
         printf("%.5E %.5E %.5E\n", t, dt, cfl);
         fflush(stdout);
      }

      timer.Stop();
      if (Mpi::Root()) {
         std::cout << "Time per iter: " << timer.RealTime() << std::endl;
      }

   }

   flowsolver.PrintTimingData();

   if (write_vec) {
      char fNameOut_u[7];
      char fNameOut_p[7];
      sprintf(fNameOut_u, "u_gf.%d", pmesh->GetMyRank());
      sprintf(fNameOut_p, "p_gf.%d", pmesh->GetMyRank());
      std::ofstream of(fNameOut_u);
      u_gf->Save(of);
      of.close();
      of.open(fNameOut_p);
      p_gf->Save(of);
      of.close();
   }

   delete pmesh;

   return 0;
}
