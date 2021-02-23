#include "dempacking.h"
//#include "Particle.h"
//#include "FreeArrays.h"

void SetParNumInDomain(uint *nn, uint *m_maxIlistntot)
{
	   if(DivX==1)
       *nn= (int) ((m_gridSize[0]+4.0)/(m_gridSize[0]-2.0)*1.2*ntot_total/numprocs); // need to modify when using loading balance
       else if(DivY==1)
       *nn= (int) ((m_gridSize[1]+4.0)/(m_gridSize[1]-2.0)*1.2*ntot_total/numprocs);
       else
       *nn= (int) ((m_gridSize[2]+4.0)/(m_gridSize[2]-2.0)*1.2*ntot_total/numprocs);

	   if(numprocs==1)*nn=ntot_total;
	   *m_maxIlistntot=*nn*10;
}

 void AllocateGPUParticleArrays()
{
	if(m_numParticles>0)
   {
    allocatenum=m_numParticles;
   }
   else
   {
    allocatenum=1;
   }
 allocatenei=allocatenum*maxCnPerParticle;
  m_oldnumParticles=m_numParticles;
   if(m_totalcontacts>0)
     m_maxIlist=m_totalcontacts;
   else
     m_maxIlist=allocatenum*10;
   m_hprehead=0;
   m_hmaxhz_total=m_hGlim.poszmin;

   //note here m_numParticles is replaced by allocatenum
   AllocateParticleDeviceArrays(allocatenum,allocatenei);
   AllocateNeigListDeviceArrays(m_maxIlist);
   AllocateParticleInCellDeviceArray();
}

void CheckIfEnoughAllocateNum()
{
 if((m_hPidrh>allocatenum && m_numParticles>0 && m_hPidrh !=0xffffffff) || (allocatenum < m_numParIn+Nlh1+Nlb1+Nrb1+Nrh1))
   { 
   // printf("m_hPidrh>allocatenum, new rank=%4d,it=%8d,m_numParticles=%8d,m_hprehead=%8d,m_hPidrh=%8d,allocatenum=%8d,m_numParIn=%8d\n",
    //   rank,it,m_numParticles,m_hprehead,m_hPidrh,allocatenum,m_numParIn); 

	  if(ntot<allocatenum)
	   {
         ntot=(int) (1.2*allocatenum);
		 FreeCPUParticleArrays();
         AllocateCPUParticleArray(ntot);
         InitialCPUParticleArray(ntot);
		 }

	  if(m_maxIlistntot<m_maxIlist)
	   {
         m_maxIlistntot=(int) (1.2*m_maxIlist);
		 FreeCPUNeigListArrays();
         AllocateCPUNeigListArray(m_maxIlistntot);
         InitialCPUNeigListArray(m_maxIlistntot);
		 }

     cpuErrchk(cudaMemcpy(m_hParticleHash,m_dParticleHash[0],allocatenum*2*sizeof(uint),cudaMemcpyDeviceToHost));
	 CopyParticleArrayDeviceToHost(0, allocatenum);
	 CopyNeigListArrayDeviceToHost(m_totalcontacts);

     FreeParticleDeviceArray();
	 FreeNeigListDeviceArrays();

     oldallocatenum=allocatenum;
     oldm_maxIlist=m_maxIlist;
     allocatenum=(int)(1.2*m_hPidrh);
     if(allocatenum < m_numParIn+Nlh1+Nlb1+Nrb1+Nrh1)
     allocatenum=m_numParIn+Nlh1+Nlb1+Nrb1+Nrh1;

     if(allocatenum<m_hPidrh) allocatenum=m_hPidrh;

     if(layer>1 && m_numParIn !=0)
      AverIlist=m_totalcontacts/m_numParIn+3;
     else
      AverIlist=6;
     if(allocatenum>ntot_total)allocatenum=ntot_total;
     m_maxIlist=(allocatenum)*AverIlist;

   // printf("re allocatenum, new rank=%4d,it=%8d,m_numParticles=%8d,m_hprehead=%8d,m_hPidrh=%8d,allocatenum=%8d\n",
    //    rank,it,m_numParticles,m_hprehead,m_hPidrh,allocatenum);

    /*if(m_numParticles >0 && (allocatenum> 10*m_numParticles || AverIlist>30))
    printf("wrong re allocatenum,rank=%4d,it=%8d,m_numParticles=%8d,m_hprehead=%8d,m_hPidrh=%8d,allocatenum=%8d,m_totalcontacts=%8d,AverIlist=%8d,m_maxIlist=%8d\n",
        rank,it,m_numParticles,m_hprehead,m_hPidrh,allocatenum,m_totalcontacts,AverIlist,m_maxIlist);*/

    if(m_maxIlist>15*m_numParticles)m_maxIlist=15*m_numParticles;
    if(m_maxIlist<m_totalcontacts)m_maxIlist=m_totalcontacts;

   allocatenei=allocatenum*maxCnPerParticle;

// allocate memery for new array with new size 

  // AllocateParticleDeviceArrays(allocatenum);
   AllocateParticleDeviceArrays(allocatenum, allocatenei);

   AllocateNeigListDeviceArrays(m_maxIlist);

   cpuErrchk(cudaMemset(m_dParticleHash[0],0xffffffff,allocatenum*2*sizeof(uint)));
   cpuErrchk(cudaMemset(m_dParticleHash[1],0xffffffff,allocatenum*2*sizeof(uint)));

   if(m_numParIn>0 && m_numParIn<UINT_MAX)
   cpuErrchk(cudaMemcpy(m_dParticleHash[0],m_hParticleHash,m_numParIn*2*sizeof(uint),cudaMemcpyHostToDevice));

   CopyParticlePropertyHostToDevice(oldallocatenum);
   CopyNeigListArrayHostToDevice(m_totalcontacts);

   if(allocatenum>0)
   cpuErrchk(cudaMemset(m_dqi,0,allocatenum*3*sizeof(double)));

   if(allocatenum<m_hPidrh)
   {printf("wrong,allocatenum<m_hPidrh,it=%8d,rank=%8d,m_numParticles=%8d,%8d,%8d\n",
   it,rank,m_numParticles,allocatenum,m_hPidrh);
    WriteParticleForRerun();
   exit(1);}

   }
   }

void ReOrderParticleArrays()
{
  //printf("rank=%4d,Pidlh=%8d,Pidlb=%8d,Pidrb=%8d,Pidrbs=%8d,Pidrh=%8d,Nlh=%8d,Nlb=%8d,Nrb=%8d,Nrh=%8d,Nlh1=%8d,Nlb1=%8d,Nrb1=%8d,Nrh1=%8d,m_numParIn=%8d,allocatenum=%8d\n",
  //        rank,m_hPidlh,m_hPidlb,m_hPidrb,m_hPidrbs,m_hPidrh,Nlh,Nlb,Nrb,Nrh,Nlh1,Nlb1,Nrb1,Nrh1,m_numParIn,allocatenum);
  //need to check allocation of arrays
   resortParticleId(m_numParIn,
	     maxCnPerParticle,
	     m_hprehead,
            idelete,
            Nlh,Nrh,
            Nlh1,Nlb1,Nrb1,Nrh1,
            m_doldIP,
            m_dParticleHash[0], m_dParticleHash[1],
	     m_dpos, m_doldpos,
	     m_dpdis, m_doldpdis,
	     m_dangv, m_doldangv,
	     m_drad, m_doldrad,
	     m_drmass, m_doldrmass,
	     m_dinert, m_doldinert,
	     m_dmatId,m_doldmatId,
            //m_dsjgi,   m_doldsjgi, 
            m_ddisptw, m_dolddisptw,
            m_dfricpw, m_doldfricpw,
			m_dcontactEnd, m_doldcontactEnd,
		    m_dEngdisp, m_doldEngdisp,
		    m_dEngdispw, m_doldEngdispw,
		    m_dEngdispVarw, m_doldEngdispVarw);

    m_oldnumParticles=m_numParticles;
    m_holdprehead=m_hprehead;

   if(m_numParIn !=m_numParticles && m_numParticles>0) 
   {
   //printf("m_numParIn !=m_numParticles,it=%8d,rank=%8d,m_numParIn=%8d,m_numParticles=%8d\n",it,rank,m_numParIn,m_numParticles);
    m_numParticles=m_numParIn;
   }

   if(m_numParticles> allocatenum)
   printf("wrong allocate sjgi,it=%8d,rank=%8d,m_numParIn=%8d,m_numParticles=%8d,allocatenum\n",
          it,rank,m_numParIn,m_numParticles,allocatenum);

//new index point to the beginning index of lh,lb,rb,rh, etc. for arraies after reordering
if(m_hPidrh>0 && m_hPidrh != 0xffffffff && m_numParticles>0) // swap old and new pointers address
{
if(allocatenum<m_hPidrh)
 printf("reorder wrong,allocatenum<m_hPidrh,it=%8d,rank=%8d,m_numParIn=%8d,m_numParticles=%8d,%8d,%8d\n",it,rank,m_numParIn,m_numParticles,allocatenum,m_hPidrh);

    std::swap(m_dParticleHash[0],m_dParticleHash[1]);
    std::swap(m_dpos,m_doldpos);
    std::swap(m_dpdis,m_doldpdis);
    std::swap(m_dangv,m_doldangv);

    std::swap(m_drad,m_doldrad);
    std::swap(m_drmass,m_doldrmass);
    std::swap(m_dinert,m_doldinert);

    std::swap(m_dmatId,m_doldmatId);

   // std::swap(m_dsjgi,m_doldsjgi);
    std::swap(m_ddisptw,m_dolddisptw);
    std::swap(m_dfricpw,m_doldfricpw);

	std::swap(m_dcontactEnd, m_doldcontactEnd);
	std::swap(m_dEngdisp, m_doldEngdisp);
	std::swap(m_dEngdispw, m_doldEngdispw);
	std::swap(m_dEngdispVarw, m_doldEngdispVarw);
}

}

void ReAllocateAneiArrays()
{
  CopyNeigListArrayDeviceToHost(m_maxIlist);
  FreeNeigListDeviceArrays();

  m_maxIlist=m_totalcontacts;
// allocate new space 
   AllocateNeigListDeviceArrays(m_maxIlist);
   CopyNeigListArrayHostToDevice(m_maxIlist);
}
