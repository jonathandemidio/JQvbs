#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <iomanip>

//Class to read in the simulation parameters from a file
class PARAMS
{
  public:
  int L1_;
  int L2_; 
  double J_; // J for projection operator
  double Q_; // Q: product of two J's
  double delta_; // field coupling to VBS order parameter: breaks lattice translation symmetry
  double deltamin_;
  double deltamax_;
  int Ndelta_; //I'll be simulating many deltas on one node
  double Beta_;  //the inverse temperature
  int suN_; // N of SU(N) symmetry of model
  int EQL_; //the number of equilibriation MC steps
  int MCS_; //the number of production MC steps
  int nBin_; // # of Bins (for statistics)
  long SEED_; //Seed

  double JbN_; // J/suN
  double deltabN_;// delta/suN => bond anisotropy
  double QbNN_;// Q/(suN*suN)

  PARAMS(int my_rank){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");    
    pfin >> L1_;
    pfin >> L2_;
    pfin >> J_;
    pfin >> Q_;
    pfin >> deltamin_;
    pfin >> deltamax_;
    pfin >> Ndelta_;
    pfin >> Beta_;
    pfin >> suN_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin.close();


    vector<double> dvec;
    dvec.resize(Ndelta_,0.0);
    double tmpd=deltamin_;
    double fac;
    if(Ndelta_>1)
        fac = pow(deltamax_/deltamin_, 1.0/(1.0*Ndelta_ - 1.0));
    else
        fac = 1.0;

    for(int d=0; d<Ndelta_; d++){
        dvec[d]=tmpd;
        tmpd=tmpd*fac;
    }

    delta_=dvec[my_rank%Ndelta_];


    JbN_=J_/suN_;
    deltabN_=delta_/suN_;
    QbNN_=Q_/(suN_*suN_);

    // print out parameter file to make sure its read correctly

    cout <<"SSE simulations of SU(N) square lattice JQ simulations in a field" <<endl;  
    cout <<"(Lx, Ly) = ("<<L1_<<", "<<L2_<<")"<<endl;
    cout <<"J = "<<J_<<endl;
    cout <<"Q = "<<Q_<<endl;
    cout <<"J/Q = "<<J_/Q_<<endl;
    cout <<"deltamin = "<<deltamin_<<" deltamax = "<<deltamax_<<" Ndelta = "<<Ndelta_<<endl;
    cout <<"my delta = "<<delta_<<endl;
    cout <<"Beta = "<<Beta_<<", SU(N), N = "<<suN_<<endl;
    cout <<"EQL steps = "<<EQL_<<endl;
    cout <<"number of steps in a bin = "<<MCS_<<", number of bins = "<<nBin_<<endl;
    cout <<"SEED = "<<SEED_<<endl;

    if(suN_>2)
      {cout <<"Not ready to simulate 2<N<5, exiting..."<<endl;exit(1);}
    
    
  }//constructor

}; //PARAMS

#endif
