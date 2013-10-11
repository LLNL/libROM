// Returns the inner product of the 2 distributed vectors v1 and v2.
double
inner_product(
   const double* v1,
   const double* v2,
   int dim,
   int num_procs);

// Returns the norm of the distributed vector v.
double
norm(
   const double* v,
   int dim,
   int num_procs);

// Normalizes the distributed vector v and returns its norm.
double
normalize(
   double* v,
   int dim,
   int num_procs);
