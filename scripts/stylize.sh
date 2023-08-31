#!/bin/bash
if [ "$1" == "" ] || [ $# -gt 2 ]; then
    echo "Usage: ./stylize.sh [-f] /path/to/astyle/executable"
    echo "Checks if stylization is required."
    echo "  -f: enforce stylization. By default, only the check is executed."
    exit 1
fi
ASTYLE_BIN=$1
enforce=false

# parse the flags.
while getopts f: flag
do
    case "${flag}" in
        f)
            enforce=true
            ASTYLE_BIN=$2
            ;;
    *)
        echo "Unknown option."
        exit 1
      ;;
    esac
done

# astyle version check
ASTYLE_VER="Artistic Style Version 3.1"
astyle_version="$($ASTYLE_BIN --version)"
if [ "$astyle_version" != "$ASTYLE_VER" ]; then
    printf "%s\n" "Invalid astyle version: '$astyle_version'"\
            "Please use: '$ASTYLE_VER'"
    exit 1
fi

# astyle dry-run if not enforced.
if [ $enforce != true ]; then
    ASTYLE_BIN="$ASTYLE_BIN --dry-run"
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="$(dirname "$DIR")"

FILELIST=("lib/*.cpp"
          "lib/*.h"
          "lib/*.hpp"
          "lib/*.c"
          "lib/*h.in"
          "regression_tests/*.cpp"
          "unit_tests/*.cpp"
          "examples/*.cpp"
          "examples/*.hpp")

ASTYLE_COMMAND="$ASTYLE_BIN --recursive --project=../.astylerc"

result=false
for files in ${FILELIST[@]}
do
    echo $files
    if $ASTYLE_COMMAND "$DIR/$files" | grep "Formatted"; then
        result=true
    fi
done

if [ $enforce != true ] && [ $result = true ]; then
    echo "Files need stylization!"
    echo "Please run stylization before merging the pull-request."
    echo "    1. Install $ASTYLE_VER"
    echo "    2. cd scripts && ./stylize.sh -f /path/to/astyle"
    exit 1
fi
