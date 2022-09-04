
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    double offlineTime = stod(argv[1]);
    double onlineTime = stod(argv[2]);
    double speedupTol = stod(argv[3]);
    if (offlineTime / onlineTime < speedupTol) {
        cerr << "offlineTime / onlineTime = " << offlineTime / onlineTime << endl;
        cerr << "speedupTol: " << speedupTol << " was not surpassed." << endl;
        abort();
    }

    return 0;
}