#include "slepc.h"

namespace CAROM {

// A singleton class to manage instances of entities which make use of the
// SLEPc library.  This class initializes SLEPc the first time that an entity
// making use of SLEPc is created and finalizes SLEPc the last time that one
// of these entities is destroyed.  It allows initialization and finalization
// of SLEPc to be hidden from applications but to occur at the necessary times.
class SLEPcManager
{
   public:
      // Return the pointer to the single instance of the SLEPcManager.
      // All access to the SLEPcManager is through getManager().
      static SLEPcManager*
      getManager()
      {
         if (!s_manager_instance) {
            s_manager_instance = new SLEPcManager;
         }
         return s_manager_instance;
      }

      // Registers an instance of an entity making use of the SLEPc
      // library.  If this is the first instance, SLEPc is initialized.
      void
      registerSLEPcInstance(int* argc, char*** argv)
      {
         if (d_instance_ct == 0) {
            SlepcInitialize(argc, argv, PETSC_NULL, PETSC_NULL);
         }
         ++d_instance_ct;
      }

      // Unregisters an instance of an entity making use of the SLEPc
      // library.  If this is the last instance, SLEPc is finalized.
      void
      unRegisterSLEPcInstance()
      {
         --d_instance_ct;
         if (d_instance_ct == 0) {
            SlepcFinalize();
         }
      }

   private:
      // Constructor is private.  This is consistent with the definition of a
      // singleton class.
      SLEPcManager();

      // Destructor is private.  This is consistent with the definition of a
      // singleton class.
      ~SLEPcManager();

      // The singleton instance.
      static SLEPcManager* s_manager_instance;

      // The number of entities needing the SLEPc library.
      int d_instance_ct;
};

}
