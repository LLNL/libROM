// Multiplies 2 matrices which are completely local to the calling process.
// Returns a matrix which is completely local to that process.
double*
LocalMatLocalMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols);

// Multiplies the transpose of a matrix, a, with another matrix, b, both of
// which are completely local to the calling process.  Returns a matrix which
// is completely local to that process.
double*
LocalMatTransposeLocalMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols);

// Multiplies a matrix whose rows are distributed across multiple processes, a,
// with a matrix which is completely local to the calling process, b.  Returns
// the part of the distributed result local to the calling process.
double*
DistributedMatLocalMatMult(
   const double* local_a,
   int num_local_a_rows,
   int num_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols);

// Multiplies the transpose of a matrix whose rows are distributed across
// multiple processes, a, with another matrix whose rows are distributed across
// multiple processes, b.  Returns the part of the distributed result local to
// the calling process.
double*
DistributedMatTransposeDistributedMatMult(
   const double* local_a,
   int num_local_a_rows,
   int num_local_a_cols,
   const double* local_b,
   int num_local_b_rows,
   int num_local_b_cols,
   int num_procs);

// Multiplies the transpose of a matrix which is completely local to the
// calling process with another matrix whose rows are distributed across
// multiple processes.  Returns the part of the distributed result local to the
// calling process.
double*
LocalMatTransposeDistributedMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* local_b,
   int num_local_b_rows,
   int num_local_b_cols,
   int num_procs,
   int rank);
