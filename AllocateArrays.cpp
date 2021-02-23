#include "dempacking.h"

void AllocateCPUParticleArray(uint nn)
{
//----- particle property -----------------------------------------------
	  m_hpos=(double *)malloc(nn*3*sizeof(double));
	  m_hpdis=(double *)malloc(nn*3*sizeof(double));
	  m_hangv=(double *)malloc(nn*3*sizeof(double));

 	  m_hrad=(double *)malloc(nn*sizeof(double));
	  m_hinert=(double *)malloc(nn*sizeof(double));
	  m_hrmass=(double *)malloc(nn*sizeof(double));
	  m_hmatId=(uint *)malloc(nn*sizeof(uint));

	   m_holdsjgi=(uint *)malloc(nn*sizeof(uint));
       if(m_holdsjgi==NULL)exit(1);
	   m_holddisptw=(double *)malloc(nn*3*sizeof(double));
       if(m_holddisptw==NULL)exit(1);
	   m_holdfricpw=(double *)malloc(nn*sizeof(double));
       if(m_holdfricpw==NULL)exit(1);

        m_hforcei=(double *)malloc(nn*3*sizeof(double));
        if(m_hforcei==NULL)exit(1);
        m_hkeng=(double *)malloc(nn*sizeof(double));
        if(m_hkeng==NULL)exit(1);
	   m_hcontacti=(uint *)malloc(nn*sizeof(uint));
       if(m_hcontacti==NULL)exit(1);

	   m_hParticleHash=(uint *)malloc(nn*2*sizeof(uint));
       if(m_hParticleHash==NULL)exit(1);

       m_holdIP=(uint *)malloc(nn*sizeof(uint));
       if(m_holdIP==NULL)exit(1);

// add energy dissipition variables
	   m_hcontactEnd=(uint *)malloc(nn*3*sizeof(uint));
       if(m_hcontactEnd==NULL)exit(1);
	   m_hEngdisp=(double *)malloc(nn*3*sizeof(double));
       if(m_hEngdisp==NULL)exit(1);
	   m_hEngdispw=(double *)malloc(nn*sizeof(double));
       if(m_hEngdispw==NULL)exit(1);
	   m_hEngdispVarw=(double *)malloc(nn*4*sizeof(double));
       if(m_hEngdispVarw==NULL)exit(1);
}

void AllocateCPUNeigListArray(uint nn)
{
	   m_holdLn=(uint *)malloc(nn*sizeof(uint));
       if(m_holdLn==NULL)exit(1);
 	   m_holddispt=(double *)malloc(nn*3*sizeof(double));
       if(m_holddispt==NULL)exit(1);
  	   m_holdfricp=(double *)malloc(nn*sizeof(double));
       if(m_holdfricp==NULL)exit(1);

	   m_hforcepair=(double *)malloc(nn*3*sizeof(double));
       if(m_hforcepair==NULL)exit(1);
	   m_hpospairi=(double *)malloc(nn*3*sizeof(double));
       if(m_hpospairi==NULL)exit(1);
	   m_hpospairj=(double *)malloc(nn*3*sizeof(double));
       if(m_hpospairj==NULL)exit(1);

// add energy dissipition variables
	   m_holdcontactEndpair=(uint *)malloc(nn*sizeof(uint));
       if(m_holdcontactEndpair==NULL)exit(1);

	   m_hEngdisppair=(double *)malloc(nn*sizeof(double));
       if(m_hEngdisppair==NULL)exit(1);
	   m_hEngdispVarpair=(double *)malloc(nn*4*sizeof(double));
       if(m_hEngdispVarpair==NULL)exit(1);
}


void AllocateMPIParticleBuffer(uint nn)
{
   m_hParticleHashsbuf=(uint *)malloc(nn*2*sizeof(uint));
   if(m_hParticleHashsbuf==NULL)exit(1);
   m_hParticleHashrbuf=(uint *)malloc(nn*2*sizeof(uint));
   if(m_hParticleHashrbuf==NULL)exit(1);

   m_hpossbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hpossbuf==NULL)exit(1);
   m_hposrbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hposrbuf==NULL)exit(1);
   m_hpdissbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hpdissbuf==NULL)exit(1);
   m_hpdisrbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hpdisrbuf==NULL)exit(1);
   m_hangvsbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hangvsbuf==NULL)exit(1);
   m_hangvrbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hangvrbuf==NULL)exit(1);

   m_hradsbuf=(double *)malloc(nn*sizeof(double));
   if(m_hradsbuf==NULL)exit(1);
   m_hradrbuf=(double *)malloc(nn*sizeof(double));
   if(m_hradrbuf==NULL)exit(1);
   m_hrmasssbuf=(double *)malloc(nn*sizeof(double));
   if(m_hrmasssbuf==NULL)exit(1);
   m_hrmassrbuf=(double *)malloc(nn*sizeof(double));
   if(m_hrmassrbuf==NULL)exit(1);
   m_hinertsbuf=(double *)malloc(nn*sizeof(double));
   if(m_hinertsbuf==NULL)exit(1);
   m_hinertrbuf=(double *)malloc(nn*sizeof(double));
   if(m_hinertrbuf==NULL)exit(1);

      m_hmatIdsbuf=(uint *)malloc(nn*sizeof(uint));
   if(m_hmatIdsbuf==NULL)exit(1);
   m_hmatIdrbuf=(uint *)malloc(nn*sizeof(uint));
   if(m_hmatIdrbuf==NULL)exit(1);

   m_hdisptwsbuf=(double *)malloc(allocatebuf*3*sizeof(double));
   if(m_hdisptwsbuf==NULL)exit(1);
   m_hdisptwrbuf=(double *)malloc(allocatebuf*3*sizeof(double));
   if(m_hdisptwrbuf==NULL)exit(1);
   m_hfricpwsbuf=(double *)malloc(allocatebuf*sizeof(double));
   if(m_hfricpwsbuf==NULL)exit(1);
   m_hfricpwrbuf=(double *)malloc(allocatebuf*sizeof(double));
   if(m_hfricpwrbuf==NULL)exit(1);

// add energy dissipition variables
     m_hcontactEndsbuf=(uint *)malloc(allocatebuf*3*sizeof(uint));
   if(m_hcontactEndsbuf==NULL)exit(1);
     m_hcontactEndrbuf=(uint *)malloc(allocatebuf*3*sizeof(uint));
   if(m_hcontactEndrbuf==NULL)exit(1);

     m_hEngdispsbuf=(double *)malloc(allocatebuf*3*sizeof(double));
   if(m_hEngdispsbuf==NULL)exit(1);
        m_hEngdisprbuf=(double *)malloc(allocatebuf*3*sizeof(double));
   if(m_hEngdisprbuf==NULL)exit(1);

     m_hEngdispwsbuf=(double *)malloc(allocatebuf*sizeof(double));
   if(m_hEngdispwsbuf==NULL)exit(1);
        m_hEngdispwrbuf=(double *)malloc(allocatebuf*sizeof(double));
   if(m_hEngdispwrbuf==NULL)exit(1);

     m_hEngdispVarwsbuf=(double *)malloc(allocatebuf*4*sizeof(double));
   if(m_hEngdispVarwsbuf==NULL)exit(1);
     m_hEngdispVarwrbuf=(double *)malloc(allocatebuf*4*sizeof(double));
   if(m_hEngdispVarwrbuf==NULL)exit(1);
}
/*
void AllocateMPINeigListBuffer(uint nn)
{
   m_hLnsbuf=(uint *)malloc(nn*sizeof(uint));
   if(m_hLnsbuf==NULL)exit(1);
   m_hLnrbuf=(uint *)malloc(nn*sizeof(uint));
   if(m_hLnsbuf==NULL)exit(1);

   m_hdisptsbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hdisptsbuf==NULL)exit(1);
   m_hdisptrbuf=(double *)malloc(nn*3*sizeof(double));
   if(m_hdisptrbuf==NULL)exit(1);

   m_hfricpsbuf=(double *)malloc(nn*sizeof(double));
   if(m_hfricpsbuf==NULL)exit(1);
   m_hfricprbuf=(double *)malloc(nn*sizeof(double));
   if(m_hfricprbuf==NULL)exit(1);

   // add energy dissipition variables
   m_hEngdispVarpairsbuf=(double *)malloc(nn*4*sizeof(double));
   if(m_hEngdispVarpairsbuf==NULL)exit(1);
      m_hEngdispVarpairrbuf=(double *)malloc(nn*4*sizeof(double));
   if(m_hEngdispVarpairrbuf==NULL)exit(1);
}*/

void AllocateParticleDeviceArrays(uint nn, uint nnei)
{
   cpuErrchk(cudaMalloc((void**)&m_dpos, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dpdis, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dangv, nn*3*sizeof(double)));

   cpuErrchk(cudaMalloc((void**)&m_drad, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_drmass, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dinert, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dmatId, nn*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_doldpos, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldpdis, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldangv, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldrad, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldrmass, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldinert, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldmatId, nn*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_dParticleHash[0], nn*2*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dParticleHash[1], nn*2*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_dAnei, nnei*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dnjgi, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dnjli, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dsjgi, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_doldsjgi, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dCN, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dsumCN, nn*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_ddisptw, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dfricpw, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dolddisptw, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldfricpw, nn*sizeof(double)));

   cpuErrchk(cudaMalloc((void**)&m_dforcei, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dtorquei, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dcontacti, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dqi, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dkeng, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldIP, nn*sizeof(uint))); 

// add energy dissipition variables
   cpuErrchk(cudaMalloc((void**)&m_dcontactEnd, nn*3*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dEngdisp, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dEngdispw, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dEngdispVarw, nn*4*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dBriProb, nn*sizeof(double)));

   cpuErrchk(cudaMalloc((void**)&m_doldcontactEnd, nn*3*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_doldEngdisp, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldEngdispw, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldEngdispVarw, nn*4*sizeof(double)));

// initialization Particle Device arrays --------------------
   cpuErrchk(cudaMemset(m_dsjgi,0,nn*sizeof(uint)));
   cpuErrchk(cudaMemset(m_ddisptw,0,nn*3*sizeof(double)));
   cpuErrchk(cudaMemset(m_dfricpw,0,nn*sizeof(double)));
   cpuErrchk(cudaMemset(m_dParticleHash[0],0xffffffff,nn*2*sizeof(uint)));
   cpuErrchk(cudaMemset(m_dParticleHash[1],0xffffffff,nn*2*sizeof(uint)));

// add energy dissipition variables initialization
  cpuErrchk(cudaMemset(m_dcontactEnd,0,nn*3*sizeof(uint)));
  cpuErrchk(cudaMemset(m_dEngdisp,0,nn*3*sizeof(double)));
  cpuErrchk(cudaMemset(m_dEngdispw,0,nn*sizeof(double)));
  cpuErrchk(cudaMemset(m_dEngdispVarw,0,nn*4*sizeof(double)));
  cpuErrchk(cudaMemset(m_dBriProb,0,nn*sizeof(double)));
}

void AllocateNeigListDeviceArrays(uint nn)
{
   cpuErrchk(cudaMalloc((void**)&m_dLp, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dLn, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_doldLn, nn*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_ddispt, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dfricp, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dolddispt, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldfricp, nn*sizeof(double)));

   cpuErrchk(cudaMalloc((void**)&m_dforcepair, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dtorqpairi, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dtorqpairj, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dcontactpair, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dpospairi, nn*3*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dpospairj, nn*3*sizeof(double)));

// add energy dissipition variables
   cpuErrchk(cudaMalloc((void**)&m_doldEngdisppair, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_doldcontactEndpair, nn*sizeof(uint)));

   cpuErrchk(cudaMalloc((void**)&m_dcontactEndpair, nn*sizeof(uint)));
   cpuErrchk(cudaMalloc((void**)&m_dEngdisppair, nn*sizeof(double)));
   cpuErrchk(cudaMalloc((void**)&m_dEngdispVarpair, nn*4*sizeof(double)));

   cpuErrchk(cudaMemset(m_doldEngdisppair,0,nn*sizeof(double)));
   cpuErrchk(cudaMemset(m_doldcontactEndpair,0,nn*sizeof(uint)));

  cpuErrchk(cudaMemset(m_dcontactEndpair,0,nn*sizeof(uint)));
   cpuErrchk(cudaMemset(m_dEngdisppair,0,nn*sizeof(double)));
   cpuErrchk(cudaMemset(m_dEngdispVarpair,0,nn*4*sizeof(double)));
}

void AllocateParticleInCellDeviceArray()
{
   if(NallocateH>0)
   cudaMalloc((void**)&NumParticleInHCell, m_numgridsH[FeedId]*sizeof(uint));
   if(NallocateM>0)
   cudaMalloc((void**)&NumParticleInMCell, m_numgridsM[FeedId]*sizeof(uint));
   if(NallocateL>0)
   cudaMalloc((void**)&NumParticleInLCell, m_numgridsL[FeedId]*sizeof(uint));
}


void AllocateCollisionArray()
{
	cudaMalloc((void**)&m_dcontactSizepair, 6*Nsect*sizeof(uint));
	cpuErrchk(cudaMemset(m_dcontactSizepair,0,6*Nsect*sizeof(uint)));

   m_hcontactSizepair=(uint *)malloc(6*Nsect*sizeof(uint));
   memset(m_hcontactSizepair,0,6*Nsect*sizeof(uint));
	//AllocateContactSizePair(6,Nsect);
}