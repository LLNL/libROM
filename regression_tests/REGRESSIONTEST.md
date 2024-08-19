# Regression Test Documentation

The script will git clone the libROM master branch into dependencies/libROM and run the regression tests on the user’s local machine, once with the user’s branch, and once with the master branch and compare the results. If something fails along the way (e.g. an error bound is exceeded), the test fails with an appropriate error message and the user can look into the run directories to examine the files more closely. Since the tests are run on a user’s local workstation, both the baseline and the new user branch are run in the same environment. Currently, the regression tests will catch errors between branches, but will not catch errors relating to differing results on different machines.

The usage instructions are output for any bad invocation. Currently, there are the include(-i) and basis comparison(-x) options. Include(-i) is used if a user only wants to run a single test. The offline(-x) option is used if the user wants to additionally compare the generated offline bases. The online solutions are compared by default. 

Note that the tests will compile libROM with MFEM if the library has not been built at the time of execution of the regression tests. 

How to run the tests on LC

1. The baseline branch placed in the 'dependencies' folder. It will be automatically brought in sync with origin/master. If there are changes to either the local or remote, the branch is recompiled.
2. Run the regression tests by invoking
   `sbatch regression_tests/run_regression_tests.sh`
   from the libROM directory
3. The slurm output file will be stored in sbatch.log in the directory you ran the previous command from. View this file for the general output of the regression tests.
4. Test logs are stored in regression_tests/results.
5. Generated example artifacts are moved to their respective folder in build/examples after each run. For example, the artifacts
produced by dg_advection are moved to dg_advection_out in build/examples/dmd.
6. To find any errors, scroll to the bottom of each test file.

How to run the tests on MAC

1. Follow step 1 above.
2. `./regression_tests/run_regression_tests.sh` (from the libROM directory)
3. Follow steps 4-5 from the instructions above.


Here are some example runs:

LC Run:
`sbatch regression_tests/run_regression_tests.sh`

To run all tests:
`sbatch run_regression_tests.sh`

To run only one test(e.g. dg_advection)
`sbatch run_regression_tests.sh -i dg_advection`

To run all tests including basis comparison
`sbatch run_regression_tests.sh -x`

To run only one test with basis comparison
`sbatch run_regression_tests.sh -xi dg_advection`

How to add a new regression test

1. Create the appropriate regression test file(with a .sh suffix) under regression_tests/tests
2. Place it under PROM or DMD depending on the type of test run.
3. Source the common.sh at the top of this file.
4. Specify the commands used to run the regression test in the CMDS array.
5. Add the type of test: TYPE="DMD" or TYPE="PROM"
6. Set OFFSET to the number of header lines in this test's generated solution file. These lines will be removed before comparing solutions.
7. (Optional) Specify a custom list of files (text-only) to compare in the FILES array. When FILES is defined, only those files will be compared, all other files will be ignored.
8. Call run_tests.

