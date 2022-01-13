#ifndef HAMIL_H 
#define HAMIL_H 

#include "RKK_determ.h" 
#include "simparam.h"
#include "lattice.h"
#include "MersenneTwister.h"

// Class declarations which define the
// Matrix elements and vertices according to 
// the Hamiltonian
//

class MATRIXELEMS
{
 public:

    int SWRE[8];
    int SWREbnd[4];

    MATRIXELEMS(const PARAMS& param, const LATTICE& latt); //constructor;
    double plx(int s1, int s2, int s3, int s4); //
    double bnd(int s1, int s2);
    
    

}; 

class SSEDATAS
{
  public:
    long nn;  //number of non-zero operators = nnp + nnb
    long Ma;  //length of operator string (redundant)
    

    vector<int> oprtr; //Operator element, 1=diag bond, 4=diag plaq, 0=null, -1 = OD bond, -2=ODplaq -3=ODplaq -4=ODplaq (see flipop function)
    vector<int> loc; // bond or plaquette number as the case maybe
    vector<int> Spin;  //Sz basis state (base 0)!!!!
    vector <int> nnpos; //stores the imaginary time positions of the non-null operators (to easily find them for loop updates)

    SSEDATAS(MTRand& ran,const LATTICE& latt,const PARAMS& param);

    int flipop(const int type, const int bnd);
    void INCREASEM(const int newM);
    void print(const LATTICE& latt);


}; //SSEDATAS 





#endif
