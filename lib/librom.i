%module librom

%include /usr/tce/packages/python/python-2.7.16/lib/python2.7/site-packages/mpi4py/include/mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
%include <std_string.i>
%include <std_vector.i>
// Instantiate templates used by example
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}
%include <std_shared_ptr.i>
%include <cpointer.i>

/* Create some functions for working with "int *" */
%pointer_functions(int, intp);

/* Create some functions for working with "double *" */
%pointer_functions(double, doublep);

%{
    #define SWIG_FILE_WITH_INIT
    #include "utils/CSVDatabase.h"
    #include "utils/HDFDatabase.h"
    #include "linalg/Vector.h"
    #include "linalg/Matrix.h"
    #include "linalg/Options.h"
    #include "linalg/BasisReader.h"
    #include "linalg/BasisGenerator.h"
    #include "algo/DMD.h"
    #include "algo/AdaptiveDMD.h"
    #include "algo/ParametricDMD.h"
    #include "algo/manifold_interp/VectorInterpolator.h"
    #include "algo/manifold_interp/MatrixInterpolator.h"
    #include "hyperreduction/DEIM.h"
    #include "hyperreduction/QDEIM.h"
    #include "hyperreduction/STSampling.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%ignore Database;
%include "utils/Database.h"
%include "utils/CSVDatabase.h"
%include "utils/HDFDatabase.h"
%include "linalg/Vector.h"
%ignore EigenPair;
%ignore ComplexEigenPair;
%ignore SerialSVD;
%ignore SerialSVDDecomposition;
%ignore SymmetricRightEigenSolve;
%ignore NonSymmetricRightEigenSolve;
%ignore SpaceTimeProduct;
%include "linalg/Matrix.h"
%include "linalg/Options.h"
%include "linalg/BasisReader.h"
%include "linalg/BasisGenerator.h"
%include "algo/DMD.h"
%include "algo/AdaptiveDMD.h"
%include "algo/ParametricDMD.h"
%ignore Interpolator;
%ignore convertClosestRBFToEpsilon;
%ignore obtainRBF;
%ignore obtainRBFToTrainingPoints;
%ignore rbfWeightedSum;
%ignore solveLinearSystem;
%include "algo/manifold_interp/Interpolator.h"
%ignore obtainInterpolatedVector;
%include "algo/manifold_interp/VectorInterpolator.h"
%include "algo/manifold_interp/MatrixInterpolator.h"
%ignore RowInfo;
%ignore RowInfoMax;
%include "hyperreduction/DEIM.h"
%include "hyperreduction/QDEIM.h"
%include "hyperreduction/STSampling.h"
