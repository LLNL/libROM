#include "SLEPcManager.h"

namespace CAROM {

SLEPcManager* SLEPcManager::s_manager_instance = 0;

SLEPcManager::SLEPcManager() :
   d_instance_ct(0)
{
}

SLEPcManager::~SLEPcManager()
{
}

}
