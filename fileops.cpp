
#include "RKK_determ.h"
#include "fileops.h"

FILEOPS::FILEOPS(const int N)
{ 
  int Ttemp=N;

  //need to establish conf file name, as well as data file name.
  sprintf(dname,"%d.data",Ttemp);
  sprintf(fname,"%d.conf",Ttemp);


}//constructor


void FILEOPS::WRITECONFIG(const PARAMS p,const SSEDATAS& sdatas, const MCSTEP& mcs)
{ 
  int i,j;
  ofstream cfout;

  cfout.open(fname);

  cfout<<mcs.numWOLFF<<"\n";
  cfout <<sdatas.Ma<<"\n"; 
  
  for (i=0 ; i<sdatas.Spin.size(); i++) 
    cfout<<sdatas.Spin[i]<<" ";
  cfout<<"\n";
  
  cfout<<"-77\n";  /*end of array tag*/

  for (i=0; i < sdatas.oprtr.size(); i++)  
    cfout<<sdatas.oprtr[i]<<" ";
  cfout<<"\n";

  cfout<<"-88\n";  /*end of array tag*/

  for (i=0; i < sdatas.loc.size(); i++)  
    cfout<<sdatas.loc[i]<<" ";
  cfout<<"\n";

  cfout<<"-99\n";  /*end of array tag*/

  cfout.close();
  
} /*end WRITECONFIG*/



void FILEOPS::READCONFIG(PARAMS& p,SSEDATAS& sdatas, MATRIXELEMS& mel ,MCSTEP& mcs, const LATTICE& latt)
{
  int opi, endOFf;
  int Ncount, Ncount2;

  Ncount=Ncount2=0;


  ifstream cfin;
  
  cfin.open(fname);


  cfin >> mcs.numWOLFF;
  if (mcs.numWOLFF< 1)
    cout<<"CHECK CONFIG FILE: Nl"<<"\n";

  cfin >> sdatas.Ma;
  if (sdatas.Ma< 1)
    cout<<"CHECK CONFIG FILE: Nl"<<"\n";

  ////////////////////////////////
  sdatas.Spin.clear();
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -77)    /*EOF tag*/
      endOFf = 0;
    else sdatas.Spin.push_back(opi);
  }//end while
  if (sdatas.Spin.size() != latt.Nsite) cout<<"CONFIG Spin FUCKUP \n";
  ////////////////////////////////

  ////////////////////////////////
  sdatas.oprtr.clear();  //empty operator string
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -88)    /*EOF tag*/
      endOFf = 0;
    else {
      sdatas.oprtr.push_back(opi);
      if (opi != 0) Ncount ++;
    } 
  }  /*end while*/
  if (sdatas.Ma != sdatas.oprtr.size())  cout<<"CONFIG FILE ERROR oprtr \n"; 
  ////////////////////////////////

  ////////////////////////////////
  sdatas.loc.clear();  //empty operator string
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -99)    /*EOF tag*/
      endOFf = 0;
    else {
      sdatas.loc.push_back(opi);
      if (opi != -1) Ncount2++;
    } 
  }  /*end while*/
  ////////////////////////////////



  if (Ncount != Ncount2)  {cout<<Ncount<<" "<<Ncount2<<"CONFIG FILE ERROR N \n"; exit(1);}

  sdatas.nn = Ncount;
  sdatas.nnpos.resize(sdatas.Ma,0);//this will be created in linkoperator, no worries

  cfin.close();

}//READCONFIG

//--------------------------------------------------------------------------------
void FILEOPS::TDprint(const PARAMS p, const MCSTEP& mcs, const LATTICE& latt)
//Print thermodynamic estimators to file.
{
  ofstream dout;
  dout.open(dname,ios::app);

  //calculate a few measurements directly, the rest will be processed later
  double energy;
  energy=-mcs.nnT/(1.0*p.MCS_*p.Beta_*latt.Nsite);

  //double Ovbs; --> calculate this in post processing
  //Ovbs = (1.0/(0.5*latt.Nsite*p.Beta_*p.MCS_))*( mcs.nxe/(p.J_+p.delta_) - mcs.nxo/(p.J_-p.delta_) );

  double magsq;
  magsq=mcs.magsq/(1.0*p.MCS_*latt.Nsite*latt.Nsite);

  dout<<setprecision(8)<<energy<<" ";
  dout<<setprecision(8)<<magsq<<" ";
  //dout<<setprecision(8)<<Ovbs<<" "; //defined as <P01> - <P12>  : a difference of two bonds --> will be redefined.
  dout<<setprecision(8)<<mcs.nxe/(1.0*p.MCS_)<<" "; //number of x-even operators
  dout<<setprecision(8)<<mcs.nxo/(1.0*p.MCS_)<<" "; //number of x-odd operators
  dout<<setprecision(8)<<mcs.nxenxe/(1.0*p.MCS_)<<" "; //number of x-even operators squared
  dout<<setprecision(8)<<mcs.nxonxo/(1.0*p.MCS_)<<" "; //number of x-odd operators squared
  dout<<setprecision(8)<<mcs.nxenxo/(1.0*p.MCS_)<<" "; //the cross term
 

  // done 
  dout<<endl;
  dout.flush();
  dout.close();

}//TDprint
