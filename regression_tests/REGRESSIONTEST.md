# Regression Test Documentation

The script will git clone the libROM master branch into dependencies/libROM and run the regression tests on the user’s local machine, once with the user’s branch, and once with the master branch and compare the results. If something fails along the way (e.g. an error bound is exceeded), the test fails with an appropriate error message and the user can look into the run directories to examine the files more closely. Since the tests are run on a user’s local workstation, both the baseline and the new user branch are run in the same environment. Currently, the regression tests will catch errors between branches, but will not catch errors relating to differing results on different machines.

The usage instructions are output for any bad invocation. Currently, there are the include(-i) and basis comparison(-x) options. Include(-i) is used if a user only wants to run a single test. The offline(-x) option is used if the user wants to additionally compare the generated offline bases. The online solutions are compared by default. 

Note that the tests will compile libROM with MFEM if the library has not been built at the time of execution of the regression tests. 

How to run the tests on LC

1. It is recommended to check that your baseline branch is up-to-date with any recent libROM commits and do a git merge if necessary. In the event that your baseline in dependencies/libROM is not up to date, a pull will be performed to update it. 
2. sbatch regression_tests/run_regression_tests.sh (if in the base directory) or sbatch run_regression_tests.sh (if in the regression_tests directory). The baseline branch will be git cloned and rebuilt each time the script is called. The different tests will automatically run in parallel (can not be turned off). Look below for test options.
3. The slurm output file will be stored in sbatch.log in the directory you ran the previous command from. View this file for the general output of the regression tests. 
4. Test logs are stored in regression_tests/results.
5. Generated example artifacts are moved to their respective folder in build/examples after each run. 
6. To find the particular error and where it occurred, scroll to the bottom of each test file. 

How to run the tests on MAC

1. Follow step 1 above.
2. ./regression_tests/run_regression_tests.sh (if in the base directory) or ./run_regression_tests.sh (if in the regression_tests directory). Look below for test options.
3. Follow steps 4-5 from the instructions above.


Here are some example runs and results:

LC Run: sbatch -N 1 -t 1:00:00 -p pbatch -o sbatch.log  --open-mode truncate regression_tests/run_regression_tests.sh

./run_regression_tests.sh -> Run all tests.

./run_regression_tests.sh -i dg_advection -> Run only dg_advection.sh

