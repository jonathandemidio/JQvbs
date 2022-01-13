#ifndef FILEOPS_H
#define FILEOPS_H

#include "RKK_determ.h"
#include "hamil.h"
#include "mcstep.h"
#include <iomanip>

// Class declarations which define file
// names and operators for reading/writing
//
class FILEOPS
{
  public:
   char fname[14];  //configuration file
   char dname[14];  //data file

   FILEOPS(const int); 
   void WRITECONFIG(const PARAMS p,const SSEDATAS& sdatas, const MCSTEP& mcs);
   void READCONFIG(PARAMS& p,SSEDATAS& sdatas, MATRIXELEMS& mel ,MCSTEP& mcs, const LATTICE& latt);
   void TDprint(const PARAMS,const MCSTEP&, const LATTICE&);

};

#endif
