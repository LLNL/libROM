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
//
// Navier double shear layer example
//
// Solve the double shear problem in the following configuration.
//
//       +-------------------+
//       |                   |
//       |      u0 = ua      |
//       |                   |
//  -------------------------------- y = 0.5
//       |                   |
//       |      u0 = ub      |
//       |                   |
//       +-------------------+
//
// The initial condition u0 is chosen to be a varying velocity in the y
// direction. It includes a perturbation at x = 0.5 which leads to an
// instability and the dynamics of the flow. The boundary conditions are fully
// periodic.

#include "navier_solver.hpp"
#include <fstream>
#include "algo/DMD.h"

using namespace mfem;
using namespace navier;

struct s_NavierContext
{
   int order = 6;
   double kinvis = 1.0 / 100000.0;
   double t_final = 1000 * 1e-3;
   double dt = 1e-3;
} ctx;

void vel_shear_ic(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);

   double rho = 30.0;
   double delta = 0.05;

   if (yi <= 0.5)
   {
      u(0) = tanh(rho * (yi - 0.25));
   }
   else
   {
      u(0) = tanh(rho * (0.75 - yi));
   }

   u(1) = delta * sin(2.0 * M_PI * xi);
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   int serial_refinements = 2;

   Mesh *mesh = new Mesh("../data/periodic-square.mesh");
   mesh->EnsureNodes();
   GridFunction *nodes = mesh->GetNodes();
   *nodes -= -1.0;
   *nodes /= 2.0;

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
   NavierSolver flowsolver(pmesh, ctx.order, ctx.kinvis);
   flowsolver.EnablePA(true);

   // Set the initial condition.
   ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
   VectorFunctionCoefficient u_excoeff(pmesh->Dimension(), vel_shear_ic);
   u_ic->ProjectCoefficient(u_excoeff);

   double t = 0.0;
   double dt = ctx.dt;
   double t_final = ctx.t_final;
   bool last_step = false;

   flowsolver.Setup(dt);

   ParGridFunction *u_gf = flowsolver.GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver.GetCurrentPressure();

   ParGridFunction w_gf(*u_gf);
   flowsolver.ComputeCurl2D(*u_gf, w_gf);

   ParaViewDataCollection pvdc("shear_output", pmesh);
   pvdc.SetDataFormat(VTKFormat::BINARY32);
   pvdc.SetHighOrderOutput(true);
   pvdc.SetLevelsOfDetail(ctx.order);
   pvdc.SetCycle(0);
   pvdc.SetTime(t);
   pvdc.RegisterField("velocity", u_gf);
   pvdc.RegisterField("pressure", p_gf);
   pvdc.RegisterField("vorticity", &w_gf);
   pvdc.Save();

   // Create the DMD object.
   Vector* u = u_gf->GetTrueDofs(); // Local
   Vector* p = p_gf->GetTrueDofs();
   int ndof_u = u->Size();
   int ndof_p = p->Size();
   int ndof = ndof_u + ndof_p;
   Vector* up = new Vector(ndof);
   for (int i = 0; i < ndof_u; i++) {
      up->Elem(i) = u->Elem(i);
   }
   for (int i = 0; i < ndof_p; i++) {
      up->Elem(ndof_u+i) = p->Elem(i);
   }
   CAROM::DMD* dmd_uvp = new CAROM::DMD(ndof, dt);
   dmd_uvp->takeSample(up->GetData(), t);

   // FOM solve
   for (int step = 0; !last_step; ++step)
   {
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

      if (Mpi::Root())
      {
         printf("%11s %11s\n", "Time", "dt");
         printf("%.5E %.5E\n", t, dt);
         fflush(stdout);
      }

      u = u_gf->GetTrueDofs();
      p = p_gf->GetTrueDofs();
      for (int i = 0; i < u->Size(); i++) {
         up->Elem(i) = u->Elem(i);
      }
      for (int i = 0; i < p->Size(); i++) {
         up->Elem(u->Size()+i) = p->Elem(i);
      }
      dmd_uvp->takeSample(up->GetData(), t);
   }

   flowsolver.PrintTimingData();

   // Train DMD
   dmd_uvp->train(0.999);

   t = 0.0;
   last_step = false;

   // Set the initial condition.
   u_ic->ProjectCoefficient(u_excoeff);

   ParGridFunction *u_dmd = flowsolver.GetCurrentVelocity();
   ParGridFunction *p_dmd = flowsolver.GetCurrentPressure();

   ParGridFunction w_dmd(*u_dmd);
   flowsolver.ComputeCurl2D(*u_dmd, w_dmd);

   ParaViewDataCollection pvdc_dmd("shear_output_dmd", pmesh);
   pvdc_dmd.SetDataFormat(VTKFormat::BINARY32);
   pvdc_dmd.SetHighOrderOutput(true);
   pvdc_dmd.SetLevelsOfDetail(ctx.order);
   pvdc_dmd.SetCycle(0);
   pvdc_dmd.SetTime(t);
   pvdc_dmd.RegisterField("velocity", u_dmd);
   pvdc_dmd.RegisterField("pressure", p_dmd);
   pvdc_dmd.RegisterField("vorticity", &w_dmd);
   pvdc_dmd.Save();

   // DMD prediction
   for (int step = 0; !last_step; ++step)
   {
      if (t + dt >= t_final - dt / 2)
      {
         last_step = true;
      }

      t += dt;
      
      CAROM::Vector* sol_dmd = dmd_uvp->predict(t);
      
      for (int i = 0; i < ndof_u; i++) {
	  u->Elem(i) = sol_dmd->item(i);
      }
      for (int i = 0; i < ndof_p; i++) {
	  p->Elem(i) = sol_dmd->item(ndof_u+i);
      }
      u_dmd->SetFromTrueDofs(*u);
      p_dmd->SetFromTrueDofs(*p);

      if (step % 10 == 0)
      {
         flowsolver.ComputeCurl2D(*u_dmd, w_dmd);
         pvdc_dmd.SetCycle(step);
         pvdc_dmd.SetTime(t);
         pvdc_dmd.Save();
      }

   }

   delete pmesh;

   return 0;
}
