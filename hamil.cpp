//functions for classes found in hamil.h

#include "hamil.h"

SSEDATAS::SSEDATAS(MTRand& ran,const LATTICE& latt,const PARAMS& param)
{//constructor
 

  nn = 0;  /*Operator string empty*/
  Ma = 20; /*Starting upper limit for OpString array*/
  
  nnpos.resize(Ma,0); //time locations of all non-nulls (only ever look at first nn elements)

  for(int i=0 ; i<latt.Nsite; i++)  //random starting spin configuration
    Spin.push_back(ran.randInt(param.suN_-1)); // randomly assign spins


  for (int i=0 ; i<Ma; i++)  /*initialize operator string to 0*/
    {
      oprtr.push_back(0);// blanks
      loc.push_back(-1);
    }

}//constructor

void SSEDATAS::INCREASEM(const int newM){
//******* Increases the expansion M
  int numFill;

  numFill = newM - Ma;  /*# of fill in operators needed*/
  
  
  for(int i=0; i<numFill; i++)
    {
      oprtr.push_back(0);  //pushes blank at end
      loc.push_back(-1);  //pushes to end
    }
  
  Ma = newM;
  nnpos.resize(Ma,0);

  if (Ma != oprtr.size() ) {cout<<"ERROR_INCREASEM"<<endl;exit(1);}
  if (Ma != loc.size() ) {cout<<"ERROR_INCREASEM"<<endl; exit(1);}
  
}//INCREASEM


void SSEDATAS::print(const LATTICE& latt)
{

    for(int i=0;i<Ma;i++)
    {
        cout<<oprtr[i]<<" "<<loc[i]<<" ";
        cout <<endl;
    }
    cout <<endl;
    for(int i=0;i<latt.Nsite;i++)
        cout<<Spin[i]<<" ";
 
    cout << "|"<<nn<<" "<<Ma<<endl<<endl;

  
}


int SSEDATAS::flipop(const int type, const int bnd)
{

    //this will flip my operator type depending on which bond I enter on.

    // for J terms the flip always exchanges -1 <==> 1

    // for Q terms it depends on what bnd I enter in on.

    //labelling scheme

    //  bnd =0,1
    //    ____
    //   |  1 |
    //   ------
    //    ____      
    //   |  0 |
    //   ------
    // if 0 and 1 are both diag, type=4 (positive)
    // if 0 is diag but 1 is off-diag type = -2 
    // if 0 is off-diag but 1 is diag type = -3
    // if both are off-diag type = -4


    if(abs(type)==1)
        return -1*type;
    else{
        if(bnd==0){
            if(type==4)
                return -3;
            else if(type==-3)
                return 4;
            else if(type==-4)
                return -2;
            else
                return -4;
        }
        else{
            //bnd == 1
            if(type==4)
                return -2;
            else if(type==-2)
                return 4;
            else if(type==-4)
                return -3;
            else
                return -4;
        }
    }


}


MATRIXELEMS::MATRIXELEMS(const PARAMS& p, const LATTICE& latt)
{ //constructor
  

    SWRE[0]=1;  SWRE[1]=0;
    SWRE[2]=3;  SWRE[3]=2;
    SWRE[4]=5;  SWRE[5]=4;
    SWRE[6]=7;  SWRE[7]=6;
    
    SWREbnd[0]=1;  SWREbnd[1]=0;
    SWREbnd[2]=3;  SWREbnd[3]=2;


}


double MATRIXELEMS::plx(int s1, int s2, int s3, int s4)
{

  if((s1==s2) && (s3==s4))
    return 1.0;
  else
    return 0.0;
}

double MATRIXELEMS::bnd(int s1, int s2)
{

 if(s1==s2)
    return 1.0;
 else
   return 0.0;
}
