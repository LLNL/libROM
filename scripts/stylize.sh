if [ "$1" == "" ] || [ $# -gt 1 ]; then
    echo "Usage: ./stylize.sh /path/to/astyle/executable"
    exit 0
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="$(dirname "$DIR")"

$1 --recursive "$DIR/lib/*.C" 
$1 --recursive "$DIR/lib/*.h" 
$1 --recursive "$DIR/lib/*.hpp" 
$1 --recursive "$DIR/lib/*.c" 
$1 --recursive "$DIR/*h.in" 
$1 --recursive "$DIR/tests/*.C"
$1 --recursive "$DIR/examples/*.cpp"
$1 --recursive "$DIR/examples/*.hpp" 
