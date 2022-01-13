//functions for classes found in mcstep.h
#include "mcstep.h"

MCSTEP::MCSTEP(const LATTICE& latt, const PARAMS& p){

  isEQ = 1;
  eqlnumWOLFF=0.0;
  numWOLFF = 5;
  NLV = 2.0;
  aa = 1.2;
  


}//constructor

void MCSTEP::DIAGUPDATE(MTRand& ran,const LATTICE& latt, SSEDATAS& sdatas, MATRIXELEMS& mel,const PARAMS& p)
{
  
    int time;
    long newM;
    int ranplx,ranbnd; /*Random plaqs for Heat-Bath attempt - called b in A.S. ref.*/
    double wplx,wbnd,wnull,wden;
    double dice;



    /**********************************************/
    /** GO THRU TIME SLICES ONCE ALLOWING CHANGES BETWEEN: NULL JD AND QDD **/
    /**********************************************/
    for (time=0; time < sdatas.Ma; time++) 
    {
        if (sdatas.oprtr[time] >=0) // this means its "fully" diagonal including null operator which is 0
	    {    
	        if(sdatas.oprtr[time]==1)
	            ranbnd = sdatas.loc[time];
	        else
	            ranbnd = ran.randInt(latt.Nbond-1);	 

            wbnd=p.Beta_*latt.Nbond*latt.Jcplg[ranbnd]*mel.bnd(sdatas.Spin[latt.Bnd[ranbnd][0]],sdatas.Spin[latt.Bnd[ranbnd][1]]);
        
	        if(sdatas.oprtr[time]==4)
	            ranplx=sdatas.loc[time];
	        else
	            ranplx = ran.randInt(latt.Nplaq-1);

	        wplx=p.Beta_*latt.Nplaq*p.QbNN_*mel.plx(sdatas.Spin[latt.Plq[ranplx][0]],sdatas.Spin[latt.Plq[ranplx][1]],sdatas.Spin[latt.Plq[ranplx][2]],sdatas.Spin[latt.Plq[ranplx][3]]);

	  
	        if(sdatas.oprtr[time]==0)
	            wnull=(double)(sdatas.Ma-sdatas.nn);
	        else
	            wnull=(double)(sdatas.Ma-sdatas.nn+1);
	  
	        wden=wbnd+wplx+wnull;
	        dice = ran.randDblExc(); 
	  
	        if(dice > ((wbnd+wplx)/wden)){// INSERT NULL

	            if(sdatas.oprtr[time] > 0){ // non-null J INTERACTION
	                sdatas.nn--; //reduce nn by 1
                }

	            sdatas.oprtr[time]=0; // insert blank
	            sdatas.loc[time]=-1;
	        }
	        else if(dice > (wplx/wden)){// INSERT BND
                if(sdatas.oprtr[time] == 0)
	                sdatas.nn++; //update # of non-zero operators in string
                
	            sdatas.oprtr[time] = 1; //JD=1
	            sdatas.loc[time]=ranbnd;
	            // dont need to update bndvtx[time][2 & 3 & 6 & 7] because its a bond
	        }
	        else{ //INSERT PLX
	            if(sdatas.oprtr[time] == 0)
	                sdatas.nn++; //update # of non-zero operators in string
                
	            sdatas.oprtr[time] = 4; // Plaquette diagonal =4
	            sdatas.loc[time]=ranplx;
	        }
	    }// END DIAG OPERATORS
        ///////// THEN  IF OFF_DIAGONAL PROPOGATE STATE
        else if(sdatas.oprtr[time]==-1) { // JO
	        sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]]+1)%2;
	        sdatas.Spin[latt.Bnd[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Bnd[sdatas.loc[time]][1]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -2){ // Plaq Off-diagonal
            sdatas.Spin[latt.Plq[sdatas.loc[time]][2]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][3]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][3]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -3){
            sdatas.Spin[latt.Plq[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][1]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -4){
            sdatas.Spin[latt.Plq[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][1]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][2]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][3]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][3]]+1)%2;
        }
        else
	        {cout <<"TIME TRAVEL :: SHOULD NOT END UP HERE!!!"<<endl;exit(1);}
     
      
    }//END TIME LOOP 

  


    newM = (int)( ((double)sdatas.nn) * aa);
  
    if(isEQ==1)
    {
        if (newM  > sdatas.Ma)      /*Adjust the length of the OpString array only during EQ*/     
	    sdatas.INCREASEM(newM);
    }
    else if ( ((double)sdatas.nn)>(0.95*(double)sdatas.Ma)) // if op string gets too long during simiulation (after EQ), quit!
        {cout <<"INCREASE EQUILIBRIATION TIME ... EXITING ... SHOULDNT CHANGE M after Equilibiration"<<endl;exit(1);}  
  
}//DIAGUPDATE

void MCSTEP::LINKOPERATOR(const LATTICE& latt, SSEDATAS& sdatas)
{

  //initialize OPLinkList to -1: note because of bonds, this array is NOT fully packed!
  // some legs have -1 entries: this mean they dont exist in the "real" linked list
  // when starting a loop, should check that entered vertex is != -1
  if (isEQ == 1){
    OPLinkList.clear();
    OPLinkList.resize(8*sdatas.Ma,-1);// has to be initialized at each call
  }
  else{
    OPLinkList.assign(8*sdatas.Ma,-1);// seems to be faster to initialize this way
  }

  //has to be initialized at each call to -1
  First.assign(latt.Nsite,-1);  
  Last.assign(latt.Nsite,-1);   

  LLsize=0;  //linked list size (sdatas.oprtr[time] excluding -1 elements)
  int countnn=0;
  for(int time=0;time<sdatas.Ma;time++)
    if(sdatas.oprtr[time]!=0) {

     sdatas.nnpos[countnn]=time;
     countnn++;

     /*FIRST RECORD IN THE LINKED LIST*/
      if(abs(sdatas.oprtr[time])==1){//bnd operator
        for(int site=0;site<2;site++){
          if(First[latt.Bnd[sdatas.loc[time]][site]]==-1){
            First[latt.Bnd[sdatas.loc[time]][site]]=8*time+site;
            Last[latt.Bnd[sdatas.loc[time]][site]]=8*time+site+4;/*+4 to move up in time for a bnd*/
          }
          else{/*enter the link into the linked list*/
            OPLinkList[Last[latt.Bnd[sdatas.loc[time]][site]]]=8*time+site;		  
            OPLinkList[8*time+site]=Last[latt.Bnd[sdatas.loc[time]][site]];		  
            Last[latt.Bnd[sdatas.loc[time]][site]]=8*time+site+4;
            LLsize+=2;
          }
        }//end for site 0..2  
      }//end bnd operator
      else{  // its a plaquette operator
       //cout<<"SHOULD NOT GO HERE FOR Q=0!!!"<<endl;
        for(int site=0;site<4;site++){
          if(First[latt.Plq[sdatas.loc[time]][site]]==-1){
            First[latt.Plq[sdatas.loc[time]][site]]=8*time+site;
            Last[latt.Plq[sdatas.loc[time]][site]]=8*time+site+4;/*+4 to move up in time for a plaq*/
          }
          else{/*enter into the linked list*/
            OPLinkList[Last[latt.Plq[sdatas.loc[time]][site]]]=8*time+site;
            OPLinkList[8*time+site]=Last[latt.Plq[sdatas.loc[time]][site]];		  
            Last[latt.Plq[sdatas.loc[time]][site]]=8*time+site+4;
            LLsize+=2;
          }
        }//end for site 0..4
      }//end plaquette operator
  
    }// end if its a nullop, nothing to do

  
  /*join Last and First*/
  for(int spin=0;spin<latt.Nsite;spin++){
    if (First[spin] != -1){ 
      LLsize+=2;
      OPLinkList.at(Last[spin])=First[spin];
      OPLinkList.at(First[spin])=Last[spin];
    }
  }

  
} //LINKOPERATOR

void MCSTEP::OPERATORLOOP(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const MATRIXELEMS& mel, const PARAMS& param) 
{  

    int Nlegmax;
    int startleg;
    int time,current,lo,li;
    int loop;
    bool exitOL;
    int ranpos;
    int site;

    Nlegmax=8*sdatas.Ma; // max number of legs 	
    countVx = 0;


    for(int countW=0; countW<numWOLFF; countW++){  //number of Wolff clusters	
        if(sdatas.nn){// avoids crashing
	
            ranpos = ran.randInt(sdatas.nn-1);
            if(abs(sdatas.oprtr[sdatas.nnpos[ranpos]])==1)
                startleg = 8*sdatas.nnpos[ranpos] + 4*ran.randInt(1)+ran.randInt(1);
            else
                startleg = 8*sdatas.nnpos[ranpos] + ran.randInt(7);


		    current=startleg;// start at startleg
		    exitOL = false;
		    loop = 0;
            


		    while(loop<Nlegmax && exitOL==false)
		    {
			    loop++;
			    time=current/8; // index of interaction from 0 to Ma-1     
			    li=current%8;
	

                lo=mel.SWRE[li];
                sdatas.oprtr[time] = sdatas.flipop(sdatas.oprtr[time], (li%4)/2 );
                current=8*time+lo;
                countVx += 2;



                if( (1-2*(lo/4))*(OPLinkList[current]-current) > 0){
                    //this should be the condition for crossing the boundary
                    //flip spin state
                    if(abs(sdatas.oprtr[time])==1)
                        site= latt.Bnd[sdatas.loc[time]][lo%4];
                    else
                        site=latt.Plq[sdatas.loc[time]][lo%4];
                    sdatas.Spin[site] = (sdatas.Spin[site]+1)%2; 
                }


			    current=OPLinkList[8*time+lo];// connect thru the LinkList

			    if(current==startleg) {exitOL = true; break;} // NEED THIS!!
		    }
	  
		    if(loop==Nlegmax) 
		        cout<<"LOOP SIZES ARE GETTING TOO BIG!!!"<<endl;
		  
	    }// if nn!=0
    } //loop over wolff loops  

    //constantly adjust loop length
    if (isEQ == 1){
        ADJNL(countVx);
        eqlnumWOLFF+=1.0*numWOLFF;//get running average during equilibration
    }
  

}//OPERATORLOOP



void MCSTEP::UPDATEALPHA(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const PARAMS& param)
{

    for(int site=0; site<latt.Nsite; site++){
        if(First[site]==-1)
            sdatas.Spin[site]=ran.randInt(param.suN_-1); //random between 0 and suN-1
    }// finish updating Spin vector


}//UPDATEALPHA


void MCSTEP::ADJNL(const long tempvtx)
//adjusts the number of simulation loops
{


   if (tempvtx < LLsize)
     numWOLFF++;  //incease the number of loops                                                                               
   else if(numWOLFF> 1) /*decrease # loops*/
   numWOLFF--;
                                                                                       
}//ADJNL


void MCSTEP::MEASURE_clear(const PARAMS& p, const LATTICE& latt){
  isEQ = 0;  //equilibriation is over
  nnT =0.0;

  nxe = 0.0;
  nxo = 0.0;
  nxenxe = 0.0;
  nxonxo = 0.0;
  nxenxo = 0.0;

  magsq=0.0;

}

//---------------------------------------------------------------------------

void MCSTEP::MEASURE(MTRand& ran, const PARAMS& p, const LATTICE& latt, SSEDATAS& sdatas)
{


    nnT += 1.0*sdatas.nn;


    double tmpnxe=0.0;
    double tmpnxo=0.0;
    
    double mag0 = 0.0; //magnetization at time slice 0
    double dmag=0.0;    //change in magnetization
    double sumdmag=0.0; //for measuring the squared magnetization

    for(int site=0; site<latt.Nsite; site++)
        mag0 += (1.0*sdatas.Spin[site]-0.5);



    for(int time=0; time<sdatas.Ma; time++){
        //measure even x-operators and odd x-operators for VBS order parameter
        if(abs(sdatas.oprtr[time])==1 && sdatas.loc[time] < latt.Nsite){
            if(sdatas.loc[time]%2==0)
                tmpnxe+=1.0;
            else
                tmpnxo+=1.0;
        }
    
        // measure magnetization thoughout entire imaginary time history
        if(sdatas.oprtr[time]==-1) { // JO
            dmag -= 2.0*(2*sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]]-1);
            sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]]+1)%2;
            sdatas.Spin[latt.Bnd[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Bnd[sdatas.loc[time]][1]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -2){ // Plaq Off-diagonal
            dmag -= 2.0*(2*sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]-1);
            sdatas.Spin[latt.Plq[sdatas.loc[time]][2]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][3]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][3]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -3){
            dmag -= 2.0*(2*sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]-1);
            sdatas.Spin[latt.Plq[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][1]]+1)%2;
        }
        else if (sdatas.oprtr[time] == -4){
            dmag -= 2.0*(2*sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]-1);
            dmag -= 2.0*(2*sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]-1);
            sdatas.Spin[latt.Plq[sdatas.loc[time]][0]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][0]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][1]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][1]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][2]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][2]]+1)%2;
            sdatas.Spin[latt.Plq[sdatas.loc[time]][3]] = (sdatas.Spin[latt.Plq[sdatas.loc[time]][3]]+1)%2;
        }

        
        sumdmag += ((mag0+dmag)*(mag0+dmag))/(1.0*sdatas.Ma);

    }//end loop over time


    magsq += sumdmag;

    nxe += tmpnxe;
    nxo += tmpnxo;
    nxenxe += tmpnxe*tmpnxe;
    nxonxo += tmpnxo*tmpnxo;
    nxenxo += tmpnxe*tmpnxo;


}// end MEASUREMENT!
