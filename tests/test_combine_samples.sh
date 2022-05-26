# Tests the function combine_samples located in examples/misc.

# Subtract the mean from 2 snapshots, and compute the SVD over the mean-subtracted snapshot matrix.
./../build/examples/misc/combine_samples -m load_samples_data/sample1_snapshot load_samples_data/sample2_snapshot

# Read 2 basis files, and compute the SVD over the combined basis.
./../build/examples/misc/combine_samples -b load_samples_data/sample1_basis load_samples_data/sample2_basis

# Snapshot files can be written after data is collected using a static basis generator. 
# Given the generator:
#   std::unique_ptr<CAROM::BasisGenerator> static_basis_generator;
# And after collecting samples (see examples or tests),
# A snapshot file can be written:
#   static_basis_generator->writeSnapshot();
