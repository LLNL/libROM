#ifndef included_SnapshotDMD_h
#define included_SnapshotDMD_h

#include "DMD.h"
#include <vector>

namespace CAROM {

class Vector;
class SnapshotDMD : public DMD
{
public:

private:
    std::vector<Vector*> d_interp_snapshots;
    
    void interpolateSnapshots();
};
}
#endif