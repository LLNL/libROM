#!/bin/bash
# If using sbatch, set the following flags
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p pbatch
#SBATCH -o sbatch.log
#SBATCH --open-mode truncate
# On Linx:
# sbatch -N 1 -t 1:00:00 -p pbatch -o sbatch.log --open-mode truncate ./regression_tests/run_regression_tests.sh
# On Mac:
# ./regression_tests/run_regression_tests.sh

#echo "PWD=$PWD"
if [[ -z ${GITHUB_WORKSPACE} ]]; then
# Set GITHUB_WORKSPACE variable
  if [[ -f "$PWD/run_regression_tests.sh" ]]; then
    GITHUB_WORKSPACE=$(cd .. && pwd)
  elif [[ -f "$PWD/regression_tests/run_regression_tests.sh" ]]; then
    GITHUB_WORKSPACE=$PWD
  else
    echo "Run tests from the libROM or libROM/regression_tests directory"
    exit
  fi
fi

case "$(uname -s)" in
    Linux*)
			MACHINE="Linux";;
    Darwin*)
			MACHINE="Darwin";;
    *)
			echo "The regression tests can only run on Linux and MAC."
			exit 1
esac
echo "MACHINE = $MACHINE"
export MACHINE
# echo "GITHUB_WORKSPACE = ${GITHUB_WORKSPACE}"

# Get options
while getopts ":i:e:" o;
do
	case "${o}" in
		i)
			i=${OPTARG}
      ;;
    *)
      echo "Usage:"
      if [[ $MACHINE = "Linux" ]]; then
          echo "To run all tests:"
          echo "sbatch -N 1 -t 1:00:00 -p pbatch -o sbatch.log --open-mode truncate ./regression_tests/run_regression_tests.sh"
          echo "To run one test (example)"
          echo "sbatch -N 1 -t 1:00:00 -p pbatch -o sbatch.log --open-mode truncate ./regression_tests/run_regression_tests.sh -i dg_advection"
      elif [[ $MACHINE = "Darwin" ]]; then
        echo "./regression_tests/run_regression_tests.sh"
        echo "Example: To run a single test:"
        echo "./regression_tests/run_regression_tests.sh -i dg_advection"
      fi
      echo "For more details, refer to REGRESSIONTEST.md"
			exit 1
      ;;
    esac
done
shift $((OPTIND-1))

echo "Setting up test suite"
echo "For detailed logs of the regression tests, please check regression_tests/results."

export GITHUB_WORKSPACE
export TMPDIR=/tmp
BASELINE_DIR=${GITHUB_WORKSPACE}/dependencies
TESTS_DIR=${GITHUB_WORKSPACE}/regression_tests/tests
BUILD_DIR=${BASELINE_DIR}/build
EXAMPLES_DMD_LOCAL=${GITHUB_WORKSPACE}/build/examples/dmd
EXAMPLES_PROM_LOCAL=${GITHUB_WORKSPACE}/build/examples/prom
EXAMPLES_DMD_BASELINE=${BASELINE_DIR}/libROM/build/examples/dmd
EXAMPLES_PROM_BASELINE=${BASELINE_DIR}/libROM/build/examples/prom
DIR=$GITHUB_WORKSPACE
scriptName="Unknown"
#clean up
cd ${EXAMPLES_DMD_LOCAL} && rm -rf ./*/
cd ${EXAMPLES_PROM_LOCAL} && rm -rf ./*/
cd ${EXAMPLES_DMD_BASELINE} && rm -rf ./*/
cd ${EXAMPLES_PROM_BASELINE} && rm -rf ./*/
#echo "My current dir = $DIR"
# Compile current branch if it isn't already compiled
# (assuming the branch isn't compiled if mfem doesn't exist)
if [ ! -d "$GITHUB_WORKSPACE/dependencies/mfem" ]; then
   cd ${GITHUB_WORKSPACE}
   echo "Compile libROM from the current branch"
   ./scripts/compile.sh -m
      echo "Compile libROM of current branch - done"
    if [[ "$?" -ne 0 ]]; then
       echo "Compilation failed for the current branch"
       exit 1
    fi
else
   echo "${GITHUB_WORKSPACE}/dependencies/libROM already exists"
fi

# Compile master branch if it isn't already compiled
# (assuming that the master branch isn't compiled if libROM isn't found in the dependencies)

if [ ! -d $BASELINE_DIR/libROM ]; then # Clone master branch to baseline directory
   #echo "Creating $BASELINE_DIR"
   mkdir -p $BASELINE_DIR
   cd ${BASELINE_DIR}
   echo "Clone libROM master into baseline"
   git clone https://github.com/LLNL/libROM.git
   cd libROM/scripts
   echo "Compile libROM master"
   ./compile.sh -m
   echo "Compile libROM master - done"
    if [[ "$?" -ne 0 ]]; then
       echo "Compilation failed for libROM master"
       exit 1
    fi
else
   echo "${BASELINE_DIR}/libROM master already exists"
   cd
fi
cd $DIR
RESULTS_DIR=$DIR/regression_tests/results
if [ ! -d $RESULTS_DIR ]; then
    echo "Creating $RESULTS_DIR"
    mkdir -p $RESULTS_DIR
else
  echo "Removing old files from $RESULTS_DIR"
	rm -rf $RESULTS_DIR/*
fi

# Get the number of processors
NUM_PROCESSORS=$(getconf _NPROCESSORS_ONLN)
re='^[0-9]+$'
if ! [[ $NUM_PROCESSORS =~ $re ]] ; then
   echo "Error: $NUM_PROCESSORS is not a number"
   exit 1
fi
# echo "Number of processors = $NUM_PROCESSORS"
totalTests=0
testNum=0
testNumPass=0
testNumFail=0
cd $TESTS_DIR
type_of_tests_to_execute=(*)
if [[ -z $i ]]; then
    echo "Running all regression tests except non_linear_elasticity"
else 
    echo "Running only $i"
fi
for type_of_test in ${type_of_tests_to_execute[@]}; do
  cd $type_of_test
  all_tests=(*)
  #echo "Tests to execute = ${all_tests}"
  #echo "i = $i"
  for test in ${all_tests[@]}; do 
      scriptName=$(basename $test ".sh")
      # Run a specific test by specifying the test (without the .sh suffix)
      if [[ -n $i && ! "$i" == "$scriptName" ]]; then
         continue
      fi

      if [[ "$scriptName" == "nonlinear_elasticity" ]] ; then
         #echo "Skipping $scriptName"
         continue
      fi
          
      simulationLogFile="${RESULTS_DIR}/${scriptName}.log"
      touch $simulationLogFile
      testNum=$((testNum+1))
      ./$test "$NUM_PROCESSORS" >> $simulationLogFile 2>&1
      if [[ $? -ne 0 || "${PIPESTATUS[0]}" -ne 0 ]];  
          then
            testNumFail=$((testNumFail+1))
            echo "$testNum. $test: FAIL"   
          else
            testNumPass=$((testNumPass+1))
            echo "$testNum. $test: PASS" 
      fi
  done
  cd ..
done
totalTests=$testNum
echo "${testNumPass} passed, ${testNumFail} failed out of ${totalTests} tests"
if [[ $testNumFail -ne 0 ]]; then
	exit 1
fi
unset GITHUB_WORKSPACE