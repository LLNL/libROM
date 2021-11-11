if [ "$1" == "" ] || [ $# -gt 1 ]; then
    echo "Usage: ./stylize.sh /path/to/astyle/executable"
    exit 0
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="$(dirname "$DIR")"

$1 $DIR/*.C $DIR/*.h $DIR/*.hpp $DIR/*.c $DIR/*h.in $DIR/tests/*.C $DIR/examples/dmd/*.cpp $DIR/examples/dmd/*.hpp $DIR/examples/prom/*.cpp
