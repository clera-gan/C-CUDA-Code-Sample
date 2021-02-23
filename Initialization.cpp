#include "dempacking.h"
//#include "Initialization.h"

void SetConstants()
{
	//printf("set constant. rank=%4d\n",rank);
         gg=9.81;
         pi=atan(1.0)*4.0;
        realtime=freqpreflow;
        real_dt=dt*sqrt(diam/9.81);
        fac=(pi*denp*9.81*diam*diam*diam)/6.0;  
}

//    ******************************************************************
//    *	  initial host golbal variables                              *
//    ******************************************************************

void InitialGlobalVar()
{
		//printf("set InitialGlobalVar. rank=%4d\n",rank);
       timet=0.0;
        itime=0;
       it=0;
	   itstop=0;
//----------- particle number related -------------
       nparticledat=0;
	m_numParticles=0;
	m_numParticles_total=0;
       m_oldnumParticles=0;
       NParticlesrank=0;
      idelete=0;
       m_numParIn=0;

       m_hprehead=0;
       m_holdprehead=0;

	   allocatenum =0;

       allocatenei =0;
       oldallocatenum=0;
       m_maxIlist=0;
       oldm_maxIlist=0;
       m_totalcontacts=0;

//  MPI zones variables ----------------------------------
      Nlb=0;
      Nlh=0;
      Nrh=0;
      Nrb=0; 
      Nlb1=0;
      Nlh1=0;
      Nrh1=0;
      Nrb1=0; 
// introduce virtual process
  if(rank>0)
  left=rank-1;
  else
  left=MPI_PROC_NULL;

  if(rank<numprocs-1)
  right=rank+1;
  else
  right=MPI_PROC_NULL;

     m_hGDiv.x=DivX;
     m_hGDiv.y=DivY;
     m_hGDiv.z=DivZ;
//------------------------------------
   
	   m_hiall_total=1;
       m_hiall=1;
// feed ------------------------------------------
       FeedId=NumOfFeed;
	   tnewpar=0.002/real_dt; //feed frequency
	   Feedstart=0;

        m_hAddcount=0;
      averkeng=0.0;
      totalCN=0;
     m_hout=0;
     m_htotalout=0;

      ibatch=0;
       NallocateH=0;
       NallocateM=0;
       NallocateL=0;
	   Nallocate =0;
       totalNallocate=0;
	 // MPI start-----------------------------------------------------------------
// Energy dissipition -----------------
	   fmat=0.169;
       Wmin=249.8;
       rmin=0.0015/diam;
	   Emin=1e-20;
	   Emax=1e+1;
	   gap=0.5;
	   Nsect= (uint) ((log10(Emax)-log10(Emin))/gap+1);
	   hminflag=0;
	   itstarthmin=0;
}

void InitialCPUParticleArray(uint nn)
{
       memset(m_hpos,0,3*nn*sizeof(double));
       memset(m_hangv,0,3*nn*sizeof(double));
       memset(m_hpdis,0,nn*3*sizeof(uint));

       memset(m_hrad,0,nn*sizeof(double));
       memset(m_hinert,0,nn*sizeof(double));
       memset(m_hrmass,0,nn*sizeof(double));
       memset(m_hmatId,0,nn*sizeof(uint));

	   memset(m_holdsjgi,0,nn*sizeof(uint));
	   memset(m_holddisptw,0,nn*3*sizeof(double));
	   memset(m_holdfricpw,0,nn*sizeof(double));

	   memset(m_hforcei,0,nn*3*sizeof(double));
	   memset(m_hkeng,0,nn*sizeof(double));
	   memset(m_hcontacti,0,nn*sizeof(uint));

	   memset(m_hParticleHash,0xffffffff,nn*2*sizeof(uint));
	   memset(m_holdIP,0xffffffff,nn*sizeof(uint));

  // add energy dissipition variables
	   memset(m_hcontactEnd,0,nn*3*sizeof(uint));
	   memset(m_hEngdisp,0,nn*3*sizeof(double));
	   memset(m_hEngdispw,0,nn*sizeof(double));
	   memset(m_hEngdispVarw,0,nn*4*sizeof(double));
}

void InitialCPUNeigListArray(uint nn)
{
	   memset(m_holdLn,0,nn*sizeof(uint));

	   memset(m_holddispt,0,nn*3*sizeof(double));
	   memset(m_holdfricp,0,nn*sizeof(double));
  // add energy dissipition variables
	  memset(m_hEngdispVarpair,0,nn*4*sizeof(double));	   
}


void InitialGPUForceData()
{
   cpuErrchk(cudaMemset(m_dforcei,0, m_numParticles*3*sizeof(double)));
   cpuErrchk(cudaMemset(m_dtorquei,0, m_numParticles*3*sizeof(double)));
   cpuErrchk(cudaMemset(m_dcontacti,0, m_numParticles*sizeof(uint)));
   cpuErrchk(cudaMemset(m_dkeng,0,m_numParticles*sizeof(double)));
   //printf("InitialGPUForceData rank=%8d it=%8d\n",rank,it);
}
