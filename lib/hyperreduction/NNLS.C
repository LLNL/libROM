/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "NNLS.h"
#include <set>

using namespace std;

namespace CAROM {

void NNLS(const double tau, Matrix const& G, Vector const& w, Vector & x, int printLevel)
{
  const int n = G.numRows();
  const int m = G.numColumns();

  CAROM_VERIFY(tau > 0.0);
  CAROM_VERIFY(m == w.dim() && m == x.dim());

  Vector b(n, false);
  Vector r(n, false);
  Vector u(m, false);
  Vector z(m, false);

  G.mult(w, b);  // b = Gw

  r = b;
  x = 0.0;

  const double bnorm = b.norm();

  std::set<int> ids;

  int outer = 0;
  while (r.norm() >= tau * bnorm)
    {
      if (printLevel) cout << "NNLS Outer loop iter " << outer << ", r norm " << r.norm() << endl;
      outer++;

      G.transposeMult(r, u);

      // Set id to max_value_index(u) for entries not in ids, without absolute values.
      int id = -1;
      double maxu = 0.0;
      for (int i=0; i<m; ++i)
	{
	  auto search = ids.find(i);
	  if (search == ids.end()) // if i is not in ids
	    {
	      const double u_i = u(i);
	      if (id == -1 || u_i > maxu)
		{
		  maxu = u_i;
		  id = i;
		}
	    }
	}

      CAROM_VERIFY(id >= 0);
      ids.insert(id);

      int inner = 0;
      while (true)
	{
	  if (printLevel) cout << "NNLS inner loop iter " << inner
			       << ", size " << ids.size() << endl;
	  inner++;

	  // Extract the submatrix G_{ids}
	  Matrix Gsub(n, ids.size(), false);

	  int count = 0;
	  for (auto i : ids)
	    {
	      for (int j=0; j<n; ++j)
		{
		  Gsub(j, count) = G(j, i);
		}
	      count++;
	    }

	  CAROM_VERIFY(Gsub.numColumns() < n);

	  // Compute the pseudo-inverse of Gsub, storing its transpose in Gsub.
	  Gsub.transposePseudoinverse();

	  Vector t(Gsub.numColumns(), false);
	  Gsub.transposeMult(b, t);

	  z = 0.0;

	  double minz = t(0);
	  count = 0;
	  for (auto i : ids)
	    {
	      z(i) = t(count);

	      if (t(count) < minz)
		minz = t(count);

	      count++;
	    }

	  CAROM_VERIFY(count == t.dim());

	  if (minz > 0.0)
	    {
	      x = z;
	      break;
	    }

	  // Find max feasible step from x to z, which is at most 1.
	  double step = 1.0;
	  bool initStep = false;

	  for (auto i : ids)
	    {
	      // Solve for x + s * (z - x) = 0 => s = x / (x - z), assuming x != z.
	      if (x(i) != 0.0 && fabs(x(i) - z(i)) > 1.0e-8 && z(i) <= 0.0)
		{
		  const double s = x(i) / (x(i) - z(i));
		  if (initStep)
		    {
		      step = std::min(step, s);
		    }
		  else
		    {
		      initStep = true;
		      step = s;
		    }
		}
	    }

	  // Set x = x + step (z - x) = (1 - step) x + step z
	  x *= 1.0 - step;
	  z *= step;
	  x += z;

	  for (auto i : ids)
	    {
	      if (x(i) <= 0.0)
		{
		  ids.erase(i);
		}
	    }
	} // inner loop

      // Update r = b - Gx
      G.mult(x, r);
      r *= -1.0;
      r += b;
    }  // outer loop
}

}
