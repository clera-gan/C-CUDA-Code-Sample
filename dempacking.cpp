
/******************************************************************
mainprog.cpp:
This is a part of the DEM code programmed by Dr. Jieqing Gan for the ellipsoidal particles.
Copyright (C) 2018 Monash University. All rights reserved.
********************************************************************/
/*
#include "Boundary.h"
#include "Feed.h"
#include "FreeArrays.h"
#include "GPUSet.h"
#include "Initialization.h"
#include "Materials.h"
#include "Movement.h"
#include "mpiFunctions.h"
#include "Particle.h"
#include "ReadData.h"
#include "WriteData.h"
*/
#include "dempacking.h"
#include "dempacking.cuh"

bool myfuncomp (double i,double j) { return (i>j); }


int main(int argc, char* argv[])
{
   MPI_Init(&argc,&argv);
   StartupMPI();
   readinputparameters();
//---------------------  
   readmaterialproperties(); //get diam
   SetMatMaxRad();

   SetConstants();

   readmeshinfo();
   SetGeometryBoundary();

   readfeedinfo();
   SetFeedArrays();
//---------------------  
   readmovementinfo();

//------------------------------------------------------
      SetParNumInDomain(&ntot,&m_maxIlistntot);
      AllocateCPUParticleArray(ntot);
	  AllocateCPUNeigListArray(m_maxIlistntot);

       InitialCPUParticleArray(ntot);
       InitialCPUNeigListArray(m_maxIlistntot);

      AllocateMaterialArray();
       InitialMaterialArray();
//------------------------------------------------------
      InitialGlobalVar();

	  allocatebuf=ntot;
	  m_maxIlistbuf=allocatebuf*10;
	  AllocateMPIParticleBuffer(allocatebuf); //newly added
	  //AllocateMPINeigListBuffer(m_maxIlistbuf);
//------------------------------------------------------
     // SetGeometryBoundary();
//------------------------------------------------------
// GPU start 
     InitialCUDA();

     AllocateGPUArray();
	 AllocateCollisionArray();
     CopyConstantsToCUDA();

//--- treat mesh data, once for all-----------------------
     TreatMesh();
 
     ReadRestartData();

     //printf("rank,it=%8d,nrestart=%4d\n",rank,it,nrestart);
   
     if(nrestart==0)
     {
     ReadParticleData();

     AllocateGPUParticleArrays();
     CopyDataHostToDevice();
	 m_hiall_total=1;
      } 
//------------------------------------------------------
//    main program starts here		
//------------------------------------------------------
    start= clock();

    while( it< tstop) 
 	{  
     it=it+1;
     timet=dt*it;
     tsd=(it)*real_dt;

   //printf("rank=%8d it=%8d\n",rank,it);

// poured packing
     FeedId=IsTimeToFeed();
 
     if(FeedId<NumOfFeed)
	 {
 // add particles to system
        itime++;
   //printf("add particles rank=%8d it=%8d\n",rank,it);
        //tnewpar=m_hfeedpro[FeedId].freq/real_dt;
       //tnewpar=0.002/real_dt;

       if(it==1|| itime>tnewpar)
	 {
         itime=0;
         m_oldnumParticles=m_numParticles;

 	   if(m_oldnumParticles>0)
	   {
         m_oldnumParticles=m_numParticles;	
		 oldm_maxIlist=m_maxIlist;	

		 if(ntot<m_oldnumParticles)
		 {
         ntot=(int) (1.2*m_oldnumParticles);

		 FreeCPUParticleArrays();
         AllocateCPUParticleArray(ntot);
         InitialCPUParticleArray(ntot);
		 }

		 if(m_maxIlistntot<oldm_maxIlist)
		 {
         m_maxIlistntot=(int) (1.2*oldm_maxIlist);

		 FreeCPUNeigListArrays();
         AllocateCPUNeigListArray(ntot);
         InitialCPUNeigListArray(ntot);
		 }

         CopyDataDeviceToHost();

          FreeParticleDeviceArray();
		  FreeNeigListDeviceArrays();
		  FreeParticleInCellDeviceArray();
	 }
  // printf("AddParticlesBatch rank=%8d it=%8d\n",rank,it);

          AddParticlesBatch();
          SetDataForParticles();

       m_hiall=1;
	   m_hprehead=0;
        } //if(it==1|| itime>tnewpar)

	 mpiReduceParNum();  
    }  // if(FeedId<=NumOfFeed)
        
   //if(it%1000==0)m_holdmaxhz=m_hmaxhz;

// set force,torque,coordination number, and m_diall to zero for each time step
    InitialGPUForceData();

   if(m_hiall_total>0) 
   {
   cpuErrchk(cudaMemset(m_dqi,0,m_numParticles*3*sizeof(double)));

// GPU code part: generate neighbour list ------------------------------------
// calculate hash
    calcHash(m_hprehead,m_dpos,m_drad,//modify
               m_dParticleHash[0],
               m_worldOrigin_global,
               m_gridSize_global,
               m_cellSize,
               m_gridSizel_global,
               m_cellSizel,
               m_numParticles);
  // printf("calcHash rank=%8d it=%8d\n",rank,it);

// sort particles based on hash----------------------------------------------
    RadixSort((KeyValuePair *) m_dParticleHash[0],
                (KeyValuePair *) m_dParticleHash[1],m_numParticles, 32);

    findBoundryNum(m_hprehead,m_dParticleHash[0],m_numParticles,
                  m_Ngridrbs,m_Ngridrb,m_Ngridlb,m_Ngridrh,m_Ngridlh,
                 &m_hPidrbs,&m_hPidlb,&m_hPidrb,&m_hPidlh,&m_hPidrh,
                 &m_numParIn); //newly added

   // printf("findBoundryNum rank=%8d it=%8d,m_numParticles=%8d,m_numParIn=%8d\n",
	//rank,it,m_numParticles,m_numParIn);

	TreatBoundHalo();		 //newly added
   //printf("TreatBoundHalo rank=%8d it=%8d\n",rank,it);
		  
	CheckIfEnoughAllocateNum();
  // printf("CheckIfEnoughAllocateNum rank=%8d it=%8d\n",rank,it);

	ReOrderParticleArrays();
   //printf("ReOrderParticleArrays rank=%8d it=%8d\n",rank,it);

	if(allocatebuf<allocatenum)
	{
	FreeMPIParticleBuffer();
	AllocateMPIParticleBuffer(allocatenum);
	//printf("ReAllocateMPIBuffer rank=%8d it=%8d\n",rank,it);
	}

	/*if(m_maxIlistbuf<m_maxIlist)
	{
	FreeMPINeigListBuffer();
	AllocateMPINeigListBuffer(m_maxIlist);
	//printf("ReAllocateMPIBuffer rank=%8d it=%8d\n",rank,it);
	}*/

	mpiDataSendRecv();

	 //  printf("after mpiDataSendRecv rank=%8d it=%8d\n",rank,it);
    calculateCN(m_hprehead,m_numParticles,nsize,
            m_drad,
            m_dsize,
            m_dmaxCn,
            m_dCN 
            );

   prefixsumCN(m_numParticles, m_dCN, m_dsumCN, &totalCN);

    if(m_numParticles>0 && totalCN>allocatenei && totalCN<m_numParticles*1000 && totalCN>=0) //not enough space allocate for Anei
   {
   allocatenei=(uint)(1.1*totalCN);
   cpuErrchk(cudaFree(m_dAnei));
  // printf("it=%8d,rank=%4d,m_numParticles=%8d,totalCN=%8d,allocatenei=%8d\n",
	//   it,rank,m_numParticles,totalCN,allocatenei);
   cpuErrchk(cudaMalloc((void**)&m_dAnei,allocatenei*sizeof(uint)));
   }
	 //  printf("before findBCellStart rank=%8d it=%8d\n",rank,it);AllocateCPUNeigListArray

   // find start of each cell, parallel for all particles
   findBCellStart(NParticlesrank,m_gridStartl,m_Ngridlb,m_Ngridrbs,
	           Nlh1,Nlb1,Nrh1,Nrb1,					
	           m_hPidlh,m_hPidrh,m_hPidlb,m_hPidrb,
                  m_nGridCellsl,
                  m_dParticleHash[0],//return
                  m_dCellStart,
		    m_dCellStartB);
	 //  printf("after findBCellStart rank=%8d it=%8d\n",rank,it);

   // generate neighbour list Anei for particle i
    neighborarray(m_hprehead,
		m_numParticles,
             // idelete,
              m_gridStartl,
              m_Ngridlb,m_Ngridrbs,
              m_Ngridlh,m_Ngridrb,
              m_hPidlh,m_hPidrbs,m_hPidrb,m_hPidrh,
              Nlh1,Nlb1,Nrb,Nrh,

              m_dsize,
              m_dCN,
              m_dsumCN,
		m_dParticleHash[0],
		m_dCellStart,
		m_dCellStartB,
             m_worldOrigin_global,
			   m_gridSize_global,
               m_gridSizel_global,
               m_cellSize,
               m_cellSizel,

              maxCnPerParticle,
              allocatenei,
              m_dpos,
              m_drad,
             // return value
              m_dAnei,
              m_dnjgi,
              m_dnjli);
	 //  printf("after neighborarray rank=%8d it=%8d\n",rank,it);

   cpuErrchk(cudaMemset(m_doldsjgi,0xffffffff, m_numParticles*sizeof(uint)));
   cpuErrchk(cudaMemcpy(m_doldsjgi,m_dsjgi,m_oldnumParticles*sizeof(uint),cudaMemcpyDeviceToDevice));

      prefixsum(m_numParticles,
              m_dnjgi,
              // return value
              m_dsjgi,
              &m_totalcontacts);
 	 //  printf("after prefixsum rank=%8d it=%8d\n",rank,it);

    if(m_numParticles>0 && m_totalcontacts==0xffffffff) 
   { 
    // printf("m_numParticles>0 && m_totalcontacts==0xffffffff,it=%8d,rank=%8d,m_numParticles=%8d\n",
      //      it,rank,m_numParticles);
     m_totalcontacts=m_maxIlist;
    }
//----------------------------------------------------------------------------
// if not enough space for arraies 
      if(m_maxIlist<m_totalcontacts)
      {
        ReAllocateAneiArrays();
      }

	 pairlist(m_hprehead,
                  m_hPidlh,
                  m_hPidrb,
                  Nlb1,
                  Nrb1,

             m_oldnumParticles+m_holdprehead,
             //idelete,
             m_dCN,
             m_dsumCN,

             m_doldIP,
             m_dParticleHash[0],
             m_dsjgi,
	      m_doldsjgi,
             m_doldLn,

             m_dAnei,
             m_dnjgi,
             m_dnjli,
	      m_dolddispt,
             m_doldfricp,
             maxCnPerParticle,
             m_numParticles,
             m_maxIlist,
			 m_doldEngdisppair,
			 m_dEngdisppair,
			  m_doldcontactEndpair,
			   m_dcontactEndpair,
             // return value
             m_dLp,
             m_dLn,
	      m_ddispt,
             m_dfricp);

    // updated neighbour list
       m_hiall=0;
       if(idelete==1)idelete=0;
    }

  // printf("pairlist rank=%8d it=%8d\n",rank,it);

 // calculate contact force and torque for each candidate pair
	calculateforcepair(m_hprehead,m_dpos,m_dpdis,m_dangv,m_drad,m_drmass,m_dmatId,	
					   m_dLp,
					   m_dLn,
                                      m_doldLn,
		               m_ddispt,
                       m_dfricp,
                       m_dolddispt,
                       m_doldfricp,
					   m_doldEngdisppair,
					   m_doldcontactEndpair,
					   m_dEngdisppair,
				       m_dEngdispVarpair,//input && output
					   diam,
                      //hamaker,cutoffdis,fac,diam,icalculateFvdw,
                      //bondforceratio,bondmaxgap,icalculateFbond,
					   xmin, hmin, hmax,
                       dt,
                       // return value
                       m_dforcepair,
                       m_dtorqpairi,
					   m_dtorqpairj,
					   m_dcontactpair,
                       m_dpospairi,
                       m_dpospairj,
					   m_dcontactEndpair,
				   m_totalcontacts);
 	// printf("after calculateforcepair rank=%8d it=%8d\n",rank,it);

	if(hminflag==1) //particles drop to ship top, only record belt outlet to ship top section
	{
    calculateContactFrequency(m_hprehead, m_totalcontacts,
				   nsize, fac, diam,
				   Emin,Emax,gap,Nsect,
				   hmin,
				   hmax,
                   m_dpos,
	               m_drad,				   
				   m_dsize,
				   m_dcontactEndpair,
				   m_dEngdisppair,
                   m_dpospairi,
                   m_dpospairj,
				   m_dcontactSizepair);
	}

 // sum all the forces on particle i
    sumallforcei( m_hprehead,maxCnPerParticle,
                  m_numParticles,m_totalcontacts,
                   m_dCN,
                   m_dsumCN,
                  m_dnjgi,
                  m_dnjli,

                  m_dforcepair,
                  m_dtorqpairi,
		  m_dtorqpairj,
		  m_dcontactpair,
                  m_dAnei,
				  fmat, Wmin, rmin,diam,
				  m_dcontactEndpair,
                  m_dEngdisppair,
				  m_drad,
		  // return
                  m_dforcei,
                  m_dtorquei,
		          m_dcontacti,
		          m_dcontactEnd,
                  m_dEngdisp,
                  m_dBriProb);
 	//  printf("after sumallforcei rank=%8d it=%8d\n",rank,it);

// calculate particle-wall forces ---------------------
  calculateforcepw(m_hprehead,m_numParticles,
	    Nmesh,
         m_dmeshSize,
         m_dmeshOrigin,
         m_dmeshEnd,

        m_cellSize, 
        m_dpos,
	    m_dpdis,
	    m_drad,
	    m_drmass,
	    m_dangv,

        // wall input ------------
        m_dFaceHashStart,
        m_dFaceHash,
        m_dFnodeid,
        m_dPosnode,
        m_dsharePE,
        m_dsharePV,
        m_dEdgeHead,
        m_dVertHead, 

        m_dNfacestart,
        m_dNnodestart,
        m_dFaceHstart,
        m_dMeshstart,

        NumOfRotation,NumOfTranslation,NumOfVibration,
         it,dt,real_dt, diam,
		fmat, Wmin, rmin,
		m_dEngdispw,
		m_dEngdispVarw,
		xmin, hmin, hmax,
        // return value-------------
        m_ddisptw,
        m_dfricpw,

        m_dforcei,
        m_dtorquei,
		m_dcontactEnd,
        m_dEngdisp,
        m_dBriProb);

  // printf("calculateforcepw rank=%8d it=%8d\n",rank,it);

// update particle position, displacement,orientation,----------
// angv, and angvb, q0, etc.------------------------------------
   m_hminhz=hmax;
	updateposition(m_hprehead,m_numParticles,
		       m_dpos,
                     m_dpdis,
                       m_dangv,
                    m_drad,
                       m_drmass,
		       m_dinert,

                       m_dforcei,
                       m_dtorquei,
                       m_dqi,
                       m_dkeng,
                    radmin,
		       dt,
			   diam,
                hmin, hmax, xmin,
                   // return value
		       &m_hiall,
		       &m_hminhz,
			   		&m_hmaxhz,
                     &averkeng);

// if any of the process need to update neighbour list, update all the process
 MPI_CHECK(MPI_Allreduce(&m_hiall,&m_hiall_total,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD)); 
 MPI_CHECK(MPI_Allreduce(&m_hminhz,&m_hminhz_total,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD)); 

if(m_hminhz_total<m_meshEnd[1][2] && hminflag==0)
   {
		  itstarthmin=it;
		  hminflag=1;
		  if(rank==0) 
		  printf("particles drop to the ship bin, itstarthmin=%8d\n",itstarthmin);
	}
 //printf("updateposition rank=%8d it=%8d m_hiall_total=%8d %8d\n",rank,it,m_hiall_total,m_numParticles);
 MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
//**********************************************************************
//         store data of particle positons at different time           *
//**********************************************************************
	if(it==1 || (it%freqppor)==0 && m_numParticles>=0)
	{
          WriteDataToPpor();
		 if(m_hminhz_total<hmax-2.0 && hminflag==1)
		  WriteOutEngdsp();
	}
//-----------------------------------------------------------------
   if(it==1 || it%freqparticle==0) //freqoutput
        {
       WriteParticleToTecplot();
	   //exit(1);
       }
  if(it%freqcontfreq==0 && hminflag==1)
  {
	  WriteOutContactFrequency();
	 // exit(1);
  }
/*
   if(it==1 || it%freqftv==0)
     {
     WriteParticleForceTree();
    }*/

 //printf("WriteParticleToTecplot rank=%8d it=%8d\n",rank,it);

//-----------------------------------------------------------------
//    *      write preflow data for re-calculation       *
//-----------------------------------------------------------------
      finish = clock();
      tret=(double)(finish-start)/CLOCKS_PER_SEC;

      if((it%10000)==0)
	  {
         WriteTimesteps();
      }

      if(it%freqpreflow==0)
	  {
         WriteParticleForRerun();

        realtime=realtime+freqpreflow;

       //  WriteDataToRestart();
      }

//-----------------------------------------------------------------
//    *      stop programme temporarily if running time is 	      *
//    *       beyond 7 hours                          *             
//-----------------------------------------------------------------
       if(it>=tstop)
	   { 
        FreeGPUData();

        FreeCPUData();

        exit(1);
       } 

} // while
   JobFinish();
   exit(1);
   return 0;
}
