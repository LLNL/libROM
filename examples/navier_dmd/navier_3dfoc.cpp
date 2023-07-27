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
#include "algo/DMD.h"
#include "algo/IncrementalDMD.h"

using namespace mfem;
using namespace navier;

struct s_NavierContext
{
   int order = 4;
   double kin_vis = 0.01;
   double t_final = 1.0;
   double dt = 1e-2;
} ctx;

void vel(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   //double zi = x(2);

   double U = 0.25;

   /*
   if (xi <= 1e-8)
   {
      u(0) = 16.0 * U * yi * zi * sin(M_PI * t / 8.0) * (0.41 - yi)
             * (0.41 - zi) / pow(0.41, 4.0);
   }
   else
   {
      u(0) = 0.0;
   }
   u(1) = 0.0;
   u(2) = 0.0;
   */
   u(0) = 0.;
   u(1) = 0.;
   if (std::fabs(xi+2.5) <= 1e-5) {
      u(0) = U*yi*(3-yi)/(1.5*1.5);
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   int serial_refinements = 2;

   Mesh *mesh = new Mesh("fluid-cht.mesh");

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
   flowsolver.EnablePA(true);

   // Set the initial condition.
   ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
   VectorFunctionCoefficient u_excoeff(pmesh->Dimension(), vel);
   u_ic->ProjectCoefficient(u_excoeff);

   // Add Dirichlet boundary conditions to velocity space restricted to
   // selected attributes on the mesh.
   Array<int> attr(pmesh->bdr_attributes.Max());
   // Inlet is attribute 1.
   attr[0] = 1;
   // Walls is attribute 3.
   attr[2] = 1;
   flowsolver.AddVelDirichletBC(vel, attr);

   double t = 0.0;
   double dt = ctx.dt;
   double t_final = ctx.t_final;
   bool last_step = false;

   flowsolver.Setup(dt);

   ParGridFunction *u_gf = flowsolver.GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver.GetCurrentPressure();

   ParaViewDataCollection pvdc("cht", pmesh);
   pvdc.SetDataFormat(VTKFormat::BINARY32);
   pvdc.SetHighOrderOutput(true);
   pvdc.SetLevelsOfDetail(ctx.order);
   pvdc.SetCycle(0);
   pvdc.SetTime(t);
   pvdc.RegisterField("velocity", u_gf);
   pvdc.RegisterField("pressure", p_gf);
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
   delete u;
   delete p;

   CAROM::DMD* dmd_up;
    if (true)
    {
	CAROM::Options svd_options(ndof, 5000, -1, true);
	svd_options.setIncrementalSVD(1e-6, dt, 1e-6, 10.0, true, false, true);
	svd_options.setMaxBasisDimension(5000);
	svd_options.setStateIO(false, false);
	svd_options.setDebugMode(false);
	std::string svd_base_file_name = "";
	dmd_up = new CAROM::IncrementalDMD(ndof, dt,
					  svd_options,
					  svd_base_file_name,
					  false);
    }
    else
    {
    	dmd_up = new CAROM::DMD(ndof, dt);
    }
    dmd_up->takeSample(up->GetData(), t);

    bool dmd_trained = false;
    Vector up_dmd(ndof);
    CAROM::Vector* up_dmd_carom = NULL;

   // FOM solve
   for (int step = 0; !last_step; ++step)
   {
      if (t + dt >= t_final - dt / 2)
      {
         last_step = true;
      }

      // ROM solve
	if (dmd_trained)
	{
	    if (myid == 0){
		std::cout << "DMD exists: make prediction" << std::endl;
	    }
	    
	    CAROM::Vector* up_pre = new CAROM::Vector(up->GetData(), ndof,
			    			     true, true); // previous solution
    	    if (true) {
		up_dmd_carom = dmd_up->predict_dt(up_pre);
	        up_dmd = up_dmd_carom->getData();
	    }
	    
	    if (true) // DMD prediction is accurate
	    {
		if (myid == 0){
		    std::cout << "DMD prediction accurate: use it" << std::endl;
		}

		up = &up_dmd; // use DMD solution as new snapshot
		t += dt;
	    }
	    else // FOM solve
	    {
		if (myid == 0){
		    std::cout << "DMD prediction not accurate: call FOM" << std::endl;
		}

      		flowsolver.Step(t, dt, step);
	    	u = u_gf->GetTrueDofs();
	        p = p_gf->GetTrueDofs();
      		for (int i = 0; i < u->Size(); i++) {
         	    up->Elem(i) = u->Elem(i);
      		}
      		for (int i = 0; i < p->Size(); i++) {
         	    up->Elem(u->Size()+i) = p->Elem(i);
      		}
      		delete u;
      		delete p;
	    }
	    
	    delete up_pre;

	}

	else {
      	    flowsolver.Step(t, dt, step);
	    u = u_gf->GetTrueDofs();
            p = p_gf->GetTrueDofs();
            for (int i = 0; i < u->Size(); i++) {
         	up->Elem(i) = u->Elem(i);
             }
      	     for (int i = 0; i < p->Size(); i++) {
         	up->Elem(u->Size()+i) = p->Elem(i);
      	      }
      	     delete u;
            delete p;
	}

      if (step % 10 == 0)
      {
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
      
      dmd_up->takeSample(up->GetData(), t);
      dmd_up->train(1);

      delete up_dmd_carom;

   }

   delete pmesh;
   delete up;
   delete dmd_up;

   return 0;

   // Train DMD
   //dmd_up->train(0.999);

   t = 0.0;
   last_step = false;

   // Set the initial condition.
   u_ic->ProjectCoefficient(u_excoeff);

   ParGridFunction *u_dmd = flowsolver.GetCurrentVelocity();
   ParGridFunction *p_dmd = flowsolver.GetCurrentPressure();

   ParGridFunction w_dmd(*u_dmd);
   flowsolver.ComputeCurl2D(*u_dmd, w_dmd);

   ParaViewDataCollection pvdc_dmd("cht_dmd", pmesh);
   pvdc_dmd.SetDataFormat(VTKFormat::BINARY32);
   pvdc_dmd.SetHighOrderOutput(true);
   pvdc_dmd.SetLevelsOfDetail(ctx.order);
   pvdc_dmd.SetCycle(0);
   pvdc_dmd.SetTime(t);
   pvdc_dmd.RegisterField("velocity", u_dmd);
   pvdc_dmd.RegisterField("pressure", p_dmd);
   pvdc_dmd.RegisterField("vorticity", &w_dmd);
   pvdc_dmd.Save();

   up_dmd_carom = new CAROM::Vector(up->GetData(), ndof, true, true);
   dmd_up->projectInitialCondition(up_dmd_carom, 0);

   // DMD prediction
   for (int step = 0; !last_step; ++step)
   {
      if (t + dt >= t_final - dt / 2)
      {
         last_step = true;
      }

      t += dt;
      
      CAROM::Vector* sol_dmd = dmd_up->predict(t);
      
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

   flowsolver.PrintTimingData();

   delete pmesh;

   return 0;
}
