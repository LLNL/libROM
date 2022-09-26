# Regression Test Documentation

The script will git clone the master branch into dependencies/libROM and run the regression tests on the user’s local machine, once with the user’s branch, and once with the master branch and compare the results. The csv and test scripts in the user’s branch are used for both simulations. If something fails along the way (i.e. a file is missing, the error bound is exceeded), the test fails with an appropriate error message and the user can look into the run directories to examine the files more closely. Since the tests are run on a user’s local workstation, both the baseline and the new user branch are run in the same environment. Currently, the regression tests will catch errors between branches, but will not catch errors relating to differing results on different machines.

The usage instructions are output for any bad invocation. Currently, there is only the include, option. Include is used if a user only wants to run some tests.

How to run the tests on LC

1. Make sure your user branch is up-to-date with any recent libROM commits and do a git merge if necessary. Otherwise, your tests will most likely fail.
2. sbatch tests/run_regression_tests.sh (if in the base directory) or sbatch run_regression_tests.sh (if in the tests directory). The baseline branch will be git cloned and rebuilt each time the script is called. The different tests will automatically run in parallel (can not be turned off). Look below for test options.
3. The slurm output file will be stored in sbatch.log in the directory you ran the previous command from. View this file for the general output of the regression tests. 
4. Test commands/logs are stored in regression_tests/results.
5. To find the particular error and where it occurred, scroll to the bottom of each test file. 

How to run the tests on MAC

1. Follow step 1 above.
2. ./regression_tests/run_regression_tests.sh (if in the base directory) or ./run_regression_tests.sh (if in the regression_tests directory). Look below for test options.
3. Follow steps 4-5 from the instructions above.


Here are some example runs and results:

./run_regression_tests.sh -> Run all tests.

./run_regression_tests.sh -i dg_advection -> Run only dg_advection.sh

