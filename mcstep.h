#ifndef MCSTEP_H 
#define MCSTEP_H 

//#include "RKK_determ.h"
#include "lattice.h"
#include "hamil.h"

// Class declarations which define the
// Monte Carlo Steps
//
class MCSTEP
{
  private:
    vector<int> OPLinkList;
    vector<int> LinkVXType;
    vector<int> First;
    vector<int> Last;
    vector<int> bop;
    vector<int> bnd1;
    vector<int> bnd2;
    vector<int> Bndry_Op; //0=middle, 1=start, 2=end (for t-winding #)



    long LLsize;
    double aa;  //factor that M>n
    int isEQ; //flag to determine if equilibriating
    long countVx;  //the length of each loop in #vertex
    double ll_uneq2;
    int numOperType;
    void ADJNL(const long tmpvtx);  //function to adjuct the number of loops
    double NLV; //number of vertices to cover each MCS loop update (see ADJNL)
    bool SW;



  public:
	// number of counted loops in SW algorithm OR estimated number of loops in the Wolff algorithm
	double numLoops;
	int m; // number of free spins
    // number of loops to do in each MCS for Wolff algorithm
    int numWOLFF;
    double eqlnumWOLFF; //compute running average to get best estimate of numWOLFF after equilibration.
    //some thermodynamic estimators
    double nnT;  //for total energy


    double nxe;
    double nxo;
    double nxenxe;
    double nxonxo;
    double nxenxo;

    double magsq;
    
    MCSTEP(const LATTICE& latt, const PARAMS& p);//constructor
    
    void DIAGUPDATE(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, MATRIXELEMS& mel, const PARAMS& p);
    void LINKOPERATOR(const LATTICE&,SSEDATAS&);
    void OPERATORLOOP(MTRand& ran, const LATTICE& latt,SSEDATAS& sdatas, const MATRIXELEMS& mel, const PARAMS& p);
    void UPDATEALPHA(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const PARAMS& p);
    void MEASURE(MTRand& ran, const PARAMS& p,const LATTICE& latt, SSEDATAS& sdatas);

    //some measurement functions
    void MEASURE_clear(const PARAMS& p, const LATTICE& latt);

};
#endif
