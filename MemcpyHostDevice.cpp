#include "dempacking.h"

void CopyParticleArrayDeviceToHost(uint prehead, uint nn)
{
  cpuErrchk(cudaMemcpy(m_hpos, m_dpos+prehead*3,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost)); 
  cpuErrchk(cudaMemcpy(m_hpdis, m_dpdis+prehead*3,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hangv, m_dangv+prehead*3,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost));

   cpuErrchk(cudaMemcpy(m_hrad,m_drad+prehead,nn*sizeof(double), //check here m_oldNumParticle
	                        cudaMemcpyDeviceToHost));
   cpuErrchk(cudaMemcpy(m_hrmass,m_drmass+prehead,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));
   cpuErrchk(cudaMemcpy(m_hinert,m_dinert+prehead,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));

   cpuErrchk(cudaMemcpy(m_hmatId,m_dmatId+prehead,nn*sizeof(uint),
	                        cudaMemcpyDeviceToHost));

  cpuErrchk(cudaMemcpy(m_holdsjgi,m_dsjgi,nn*sizeof(uint),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_holddisptw, m_ddisptw+prehead*3,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_holdfricpw, m_dfricpw+prehead,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));

// add energy dissipition variables
  cpuErrchk(cudaMemcpy(m_hcontactEnd,m_dcontactEnd+prehead*3,nn*3*sizeof(uint),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hEngdisp,m_dEngdisp+prehead*3,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hEngdispw,m_dEngdispw+prehead,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hEngdispVarw,m_dEngdispVarw+prehead*4,nn*4*sizeof(double),
	                        cudaMemcpyDeviceToHost));
}

void CopyNeigListArrayDeviceToHost(uint nn)
{
  cpuErrchk(cudaMemcpy(m_holdLn, m_doldLn,nn*sizeof(uint),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_holddispt, m_dolddispt,nn*3*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_holdfricp, m_doldfricp,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));

  // add energy dissipition variables
  cpuErrchk(cudaMemcpy(m_holdcontactEndpair, m_doldcontactEndpair,nn*sizeof(uint),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hEngdisppair, m_doldEngdisppair,nn*sizeof(double),
	                        cudaMemcpyDeviceToHost));
  cpuErrchk(cudaMemcpy(m_hEngdispVarpair, m_dEngdispVarpair,nn*4*sizeof(double),
	                        cudaMemcpyDeviceToHost));
}

void CopyParticlePropertyHostToDevice(uint nn)
{
   cpuErrchk(cudaMemcpy(m_dpos, m_hpos,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dpdis, m_hpdis,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dangv, m_hangv,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_drad, m_hrad,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_drmass, m_hrmass,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dinert, m_hinert,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dmatId, m_hmatId,nn*sizeof(uint),
	                        cudaMemcpyHostToDevice));

  cpuErrchk(cudaMemcpy(m_dsjgi,m_holdsjgi,nn*sizeof(uint),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_ddisptw, m_holddisptw,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dfricpw, m_holdfricpw,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));

 // add energy dissipition variables
 cpuErrchk(cudaMemcpy(m_dcontactEnd, m_hcontactEnd,nn*3*sizeof(uint),
	                        cudaMemcpyHostToDevice));
 cpuErrchk(cudaMemcpy(m_dEngdisp, m_hEngdisp,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
 cpuErrchk(cudaMemcpy(m_dEngdispw, m_hEngdispw,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dEngdispVarw, m_hEngdispVarw,nn*4*sizeof(double),
	                        cudaMemcpyHostToDevice));
}

void CopyNeigListArrayHostToDevice(uint nn)
{
    cpuErrchk(cudaMemcpy(m_doldLn, m_holdLn,nn*sizeof(uint),
	                        cudaMemcpyHostToDevice));
    cpuErrchk(cudaMemcpy(m_dolddispt, m_holddispt,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
    cpuErrchk(cudaMemcpy(m_doldfricp, m_holdfricp,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));

// add energy dissipition variables
	cpuErrchk(cudaMemcpy( m_doldcontactEndpair, m_holdcontactEndpair,nn*sizeof(uint),
	                        cudaMemcpyHostToDevice));
	cpuErrchk(cudaMemcpy( m_doldEngdisppair, m_hEngdisppair,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
    cpuErrchk(cudaMemcpy(m_dEngdispVarpair, m_hEngdispVarpair,nn*4*sizeof(double),
	                        cudaMemcpyHostToDevice));
}

void CopyParticleBatchHostToDevice(uint nn,uint dprehead,uint hprehead)
{
   cpuErrchk(cudaMemcpy(m_dpdis+dprehead*3,m_hpdisbatch[FeedId]+hprehead*3,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dangv+dprehead*3,m_hangvbatch[FeedId]+hprehead*3,nn*3*sizeof(double),
	                        cudaMemcpyHostToDevice));

   cpuErrchk(cudaMemcpy(m_drad+dprehead,m_hradbatch[FeedId]+hprehead,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_drmass+dprehead,m_hrmassbatch[FeedId]+hprehead,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dinert+dprehead,m_hinertbatch[FeedId]+hprehead,nn*sizeof(double),
	                        cudaMemcpyHostToDevice));
   cpuErrchk(cudaMemcpy(m_dmatId+dprehead,m_hmatIdbatch[FeedId]+hprehead,nn*sizeof(uint),
	                        cudaMemcpyHostToDevice));
}


void CopyDataHostToDevice()
{
	CopyParticlePropertyHostToDevice(m_oldnumParticles);
	CopyNeigListArrayHostToDevice(m_totalcontacts);

   if(m_oldnumParticles==0)cudaMemset(NumParticleInHCell, 0, m_numgridsH[FeedId]*sizeof(uint));
   if(m_oldnumParticles==0 && NallocateL==1)
     cudaMemset(NumParticleInLCell, 0, m_numgridsL[FeedId]*sizeof(uint));
  // printf("read data right, rank=%4d\n",rank);
}	

void CopyDataDeviceToHost()
{
  CopyParticleArrayDeviceToHost(m_hprehead, m_numParticles);
  CopyNeigListArrayDeviceToHost(m_maxIlist);
}
