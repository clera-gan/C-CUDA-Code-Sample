#include "dempacking.h"
//#include "WriteData.h"
//#include "Movement.h"

void WriteDataToPpor()
{
	 MPI_CHECK(MPI_Allreduce(&m_hmaxhz,&m_hmaxhz_total,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD)); 
	
	if(rank==0)
	{
    if((fp1=fopen("ppor.dat","ab+"))==NULL)printf("Cannot open file ppor.dat!");
    tsd=(it)*real_dt;
    disrate=m_htotalout*pi*diam*diam*diam*denp/6/tsd;
    fprintf(fp1,"%8d %15.6g %8d %15.6g %15.6g %15.6g %15.6g\n",
		it,tsd,m_numParticles_total,disrate,m_hminhz_total,m_hmaxhz_total,averkeng*gg*diam);
    printf(     "%8d %15.6g %8d %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
		it,tsd,m_numParticles_total,disrate,m_hminhz_total,m_hmaxhz_total,averkeng*gg*diam,hmin,hmax,xmin);
    fclose(fp1);
	}
}

void WriteDataToRestart()
{
	if(rank==0)
	{
        if((fp1=fopen("input/restart.dat","rb+"))==NULL){
        printf("Cannot open file restart.dat!");
        exit(0);
		}
        fprintf(fp1,"%8d %8d",0,nparticledat); 
        fclose(fp1);
	}
}

void WriteParticleForceTree()
{
    double *outputforcepair,*hforcepair_total;
	double *outputpospairi,*hpospairi_total,*outputpospairj,*hpospairj_total;

     int *displs,*count;
	 uint m_totalcontactsSum;
	 MPI_Allreduce(&m_totalcontacts,&m_totalcontactsSum,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
     displs=(int *)malloc(numprocs*sizeof(int));
     count=(int *)malloc(numprocs*sizeof(int));

// forcepair---------------------------------
     outputforcepair=(double *)malloc((m_totalcontacts)*3*sizeof(double));
     if(rank==0)hforcepair_total=(double *)malloc((m_totalcontactsSum)*3*sizeof(double));
     int sendcount=(int)m_numParticles*3;
     int displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
  	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 } 
     cpuErrchk(cudaMemcpy(outputforcepair,m_dforcepair+m_hprehead*3,m_totalcontacts*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputforcepair,m_totalcontacts*3,MPI_DOUBLE,hforcepair_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputforcepair);

// pospairi---------------------------------
     outputpospairi=(double *)malloc((m_totalcontacts)*3*sizeof(double));
     if(rank==0)hpospairi_total=(double *)malloc((m_totalcontactsSum)*3*sizeof(double));
   
     cpuErrchk(cudaMemcpy(outputpospairi,m_dpospairi+m_hprehead*3,m_totalcontacts*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputpospairi,m_totalcontacts*3,MPI_DOUBLE,hpospairi_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputpospairi);

// pospairj---------------------------------
     outputpospairj=(double *)malloc((m_totalcontacts)*3*sizeof(double));
     if(rank==0)hpospairj_total=(double *)malloc((m_totalcontactsSum)*3*sizeof(double));
   
     cpuErrchk(cudaMemcpy(outputpospairj,m_dpospairj+m_hprehead*3,m_totalcontacts*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputpospairj,m_totalcontacts*3,MPI_DOUBLE,hpospairj_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputpospairj);

     free(displs);
     free(count);

      if(rank==0)
      {
      if(m_totalcontactsSum>0)
      {
//--------------------------------------------------------------------------
// force tree ----------------------
//--------------------------------------------------------------------------
      tsd=(it)*real_dt;
      sprintf(str1,"%06.3f",tsd);
      sprintf(str,"output/ftv%s.fdat",str1);

      if((fp1=fopen(str,"wb+"))==NULL){
      printf("Cannot open file forcetree.dat! %d\n",it);
      exit(0);
	  }

       for(i=0;i<m_totalcontactsSum;i++)
       {
       forcepairmag=sqrt(hforcepair_total[i*3]*hforcepair_total[i*3]+ \
		                 hforcepair_total[i*3+1]*hforcepair_total[i*3+1]+ \
                         hforcepair_total[i*3+2]*hforcepair_total[i*3+2]);

       if(fabs(forcepairmag)>1e-20){
       fprintf(fp1,"%15.6g %15.6g %15.6g %6.3g %15.6g %15.6g %15.6g  %6.3g %15.6g\n",
                hpospairi_total[i*3],hpospairi_total[i*3+1],hpospairi_total[i*3+2],0.5,
                hpospairj_total[i*3],hpospairj_total[i*3+1],hpospairj_total[i*3+2],0.5,forcepairmag);}
       }

	   free(hforcepair_total);
	   free(hpospairi_total);
	   free(hpospairj_total);

       fclose(fp1);
      }
	  }
}

void WriteParticleToTecplot()
{
	 double *outputpos,*hpos_total,*outputrad,*hrad_total;
	 double *outputVel,*hVel_total;
	 double *outputBriProb,*hBriProb_total,*outputEngdisp,*hEngdisp_total;
	 uint *outputcontactEnd,*hcontactEnd_total;

     int displsi,sendcount;
     int *displs,*count;
      outputpos=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hpos_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));


     displs=(int *)malloc(numprocs*sizeof(int));
     count=(int *)malloc(numprocs*sizeof(int));

     sendcount =(int)m_numParticles*3;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }
	
     cpuErrchk(cudaMemcpy(outputpos,m_dpos+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputpos,m_numParticles*3,MPI_DOUBLE,hpos_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputpos);

	 outputVel=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hVel_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));

	 cpuErrchk(cudaMemcpy(outputVel,m_dpdis+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputVel,m_numParticles*3,MPI_DOUBLE,hVel_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputVel);

     outputrad=(double *)malloc((m_numParticles)*sizeof(double));
     if(rank==0)hrad_total=(double *)malloc((m_numParticles_total)*sizeof(double));

     sendcount =(int)m_numParticles;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }

     cpuErrchk(cudaMemcpy(outputrad,m_drad+m_hprehead,m_numParticles*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputrad,m_numParticles,MPI_DOUBLE,hrad_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputrad);
//-----------------
	  /*outputBriProb=(double *)malloc((m_numParticles)*sizeof(double));
     if(rank==0)hBriProb_total=(double *)malloc((m_numParticles_total)*sizeof(double));

     cpuErrchk(cudaMemcpy(outputBriProb,m_dBriProb+m_hprehead,m_numParticles*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputBriProb,m_numParticles,MPI_DOUBLE,hBriProb_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputBriProb);*/

	 outputcontactEnd=(uint *)malloc((m_numParticles)*sizeof(uint));
     if(rank==0)
	 {
	 hcontactEnd_total=(uint *)malloc((m_numParticles_total)*sizeof(uint));
	 }
	 cpuErrchk(cudaMemcpy(outputcontactEnd,m_dcontactEnd+m_hprehead,m_numParticles*sizeof(uint),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputcontactEnd,m_numParticles,MPI_UNSIGNED,hcontactEnd_total,count,displs,MPI_UNSIGNED,0,MPI_COMM_WORLD)); 
     free(outputcontactEnd);
	 //---------------------------------------------------------------
	 sendcount =(int)m_numParticles*3;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }
     outputEngdisp=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hEngdisp_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));

     cpuErrchk(cudaMemcpy(outputEngdisp,m_dEngdisp+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputEngdisp,m_numParticles*3,MPI_DOUBLE,hEngdisp_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputEngdisp);

//-----------------------------------
     free(displs);
     free(count);

      if(rank==0)
      {
     totalNode=0;
     totalFace=0;

    for(i=0;i<Nmesh;i++)
    {
    totalNode+=Nnode[i];
    totalFace+=Nface[i];
    }
      if((fp1=fopen("output/particle.dat","wb+"))==NULL){
      printf("Cannot open file particle01.dat! %d\n",it);
      exit(0);
	  }

	  tsd=(it)*real_dt;
	  sprintf(str1,"%06.3f",tsd);

       if(it==1)
       {
       fprintf(fp1,"TITLE     = \"GAMBIT to Fluent File\"\n");
       fprintf(fp1,"VARIABLES = \"X\"\n");
       fprintf(fp1,"\"Y\"\n");
       fprintf(fp1,"\"Z\"\n");
	   fprintf(fp1,"\"Velocity\"\n");
       fprintf(fp1,"\"Engdisp\"\n");
	   fprintf(fp1,"\"Engdispiiw\"\n");
	   fprintf(fp1,"\"Engdispiw\"\n");
       fprintf(fp1,"\"Radius\"\n");
       fprintf(fp1,"\"CollisionNum\"\n"); 

       fprintf(fp1,"DATASETAUXDATA Common.VectorVarsAreVelocity=\"TRUE\"\n");
       }

     for(i=0;i<Nmesh;i++)
     {
      fprintf(fp1,"ZONE T=\"ZONEW%s\"\n",m_hmeshFile[i].name);
       fprintf(fp1," N=%8d, E=%8d, ZONETYPE=FEQuadrilateral\n",Nnode[i],Nface[i]);
       fprintf(fp1," DATAPACKING=POINT\n");
       fprintf(fp1," AUXDATA Common.BoundaryCondition=\"Wall\"\n");
       fprintf(fp1," AUXDATA Common.IsBoundaryZone=\"TRUE\"\n");
       fprintf(fp1," DT=(DOUBLE DOUBLE DOUBLE )\n");

       for(j=0;j<Nnode[i];j++)
       { 
   // treat particle position when they are near or inside the rotary part
      px=m_hPosnode[i][j*3];
      py=m_hPosnode[i][j*3+1];
      pz=m_hPosnode[i][j*3+2];

        int RotateId=IsTimeToRotateHost(it, real_dt,NumOfRotation);
        if(RotateId<NumOfRotation &&  m_hRotation[RotateId].meshId==i)
        {
        double px1,py1,pz1;
        double rotaterate =-m_hRotation[RotateId].rotatespeed*2*pi/60*real_dt/dt; //reduced rotate rate
        double theta=rotaterate*it*dt;
        RotationPositionConvertHost(px,py,pz,m_hRotation[RotateId].axis.x,m_hRotation[RotateId].axis.y,m_hRotation[RotateId].axis.z,theta,&px1,&py1,&pz1);
        fprintf(fp1,"%10.4f%10.4f%10.4f%2d%2d%2d%2d%2d%2d\n",px1*diam,py1*diam,pz1*diam,0,0,0,0,0,0);
        }
        else
        fprintf(fp1,"%10.4f%10.4f%10.4f%2d%2d%2d%2d%2d%2d\n",px*diam,py*diam,pz*diam,0,0,0,0,0,0);
       }

       for(j=0;j<Nface[i];j++)
       {
       fprintf(fp1,"%8d%8d%8d%8d\n",m_hFnodeid[i][j*4]+1,m_hFnodeid[i][j*4+1]+1,m_hFnodeid[i][j*4+2]+1,m_hFnodeid[i][j*4+3]+1);
       }

       } //for(i=0;

	   if(m_numParticles_total>0)
       {
       fprintf(fp1,"ZONE T=\"ZONEP%s\"\n",str1);

       for(i=0;i<m_numParticles_total;i++)
        {
        px=hpos_total[3*i];
        py=hpos_total[3*i+1];
        pz=hpos_total[3*i+2];
		double velmag=sqrt(hVel_total[3*i]*hVel_total[3*i]+hVel_total[3*i+1]*hVel_total[3*i+1]+hVel_total[3*i+2]*hVel_total[3*i+2]);
        velmag=velmag*sqrt(gg*diam)/dt;

        fprintf(fp1,"%10.4f %10.4f %10.4f %13.6g %15.6g %15.6g %15.6g %9.6f %8d\n",
		px*diam,py*diam,pz*diam,velmag,hEngdisp_total[i*3]*denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0,\
		hEngdisp_total[i*3+1]*denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0,hEngdisp_total[i*3+2]*denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0,hrad_total[i]*2.0*diam,hcontactEnd_total[i]);
	    }

	   free(hpos_total);
	   free(hrad_total);

	   free(hVel_total);
	   //free(hBriProb_total);
	   free(hEngdisp_total);
	   free(hcontactEnd_total);
	  }

       fclose(fp1);
	  }
} 

void WriteOutEngdsp()
{
	 double *outputpos,*hpos_total,*outputrad,*hrad_total;
	 double *outputVel,*hVel_total;
	 double *outputBriProb,*hBriProb_total,*outputEngdisp,*hEngdisp_total;
	 uint *outputcontactEnd,*hcontactEnd_total;

     int displsi,sendcount;
     int *displs,*count;
      outputpos=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hpos_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));

     displs=(int *)malloc(numprocs*sizeof(int));
     count=(int *)malloc(numprocs*sizeof(int));

     sendcount =(int)m_numParticles*3;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }
     cpuErrchk(cudaMemcpy(outputpos,m_dpos+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputpos,m_numParticles*3,MPI_DOUBLE,hpos_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputpos);

	 outputVel=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hVel_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));

	 cpuErrchk(cudaMemcpy(outputVel,m_dpdis+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputVel,m_numParticles*3,MPI_DOUBLE,hVel_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputVel);

     outputrad=(double *)malloc((m_numParticles)*sizeof(double));
     if(rank==0)hrad_total=(double *)malloc((m_numParticles_total)*sizeof(double));

     sendcount =(int)m_numParticles;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }
     cpuErrchk(cudaMemcpy(outputrad,m_drad+m_hprehead,m_numParticles*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputrad,m_numParticles,MPI_DOUBLE,hrad_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputrad);
//-----------------
	/*  outputBriProb=(double *)malloc((m_numParticles)*sizeof(double));
     if(rank==0)hBriProb_total=(double *)malloc((m_numParticles_total)*sizeof(double));

     cpuErrchk(cudaMemcpy(outputBriProb,m_dBriProb+m_hprehead,m_numParticles*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputBriProb,m_numParticles,MPI_DOUBLE,hBriProb_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputBriProb);*/



	 //-------------------------------------------------------------------------------
	 sendcount =(int)m_numParticles*3;
	 displsi=0;

     MPI_CHECK(MPI_Gather(&sendcount,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD));
	 for(i=0;i<numprocs;i++)
	 {
	 displs[i]=displsi;
	 displsi +=count[i];
	 }
     outputEngdisp=(double *)malloc((m_numParticles)*3*sizeof(double));
     if(rank==0)hEngdisp_total=(double *)malloc((m_numParticles_total)*3*sizeof(double));

     cpuErrchk(cudaMemcpy(outputEngdisp,m_dEngdisp+m_hprehead*3,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputEngdisp,m_numParticles*3,MPI_DOUBLE,hEngdisp_total,count,displs,MPI_DOUBLE,0,MPI_COMM_WORLD)); 
     free(outputEngdisp);

	 outputcontactEnd=(uint *)malloc((m_numParticles)*3*sizeof(uint));
     if(rank==0)
	 {
	 hcontactEnd_total=(uint *)malloc((m_numParticles_total)*3*sizeof(uint));
	 }
	 cpuErrchk(cudaMemcpy(outputcontactEnd,m_dcontactEnd+m_hprehead*3,m_numParticles*3*sizeof(uint),cudaMemcpyDeviceToHost));
     MPI_CHECK(MPI_Gatherv(outputcontactEnd,m_numParticles*3,MPI_UNSIGNED,hcontactEnd_total,count,displs,MPI_UNSIGNED,0,MPI_COMM_WORLD)); 
     free(outputcontactEnd);
//-----------------------------------
     free(displs);
     free(count);

	  if(rank==0)
      {
	   uint icount=0,icountw=0,icountwi=0;
	   uint totalcontactEnd=0,totalcontactEndii=0,totalcontactEndiw=0;
	   double totalEdisp=0.0,totalEdispw=0.0,totalEdispwi=0.0;
	   tsd=(it)*real_dt;

      uint *NumInsize =(uint *)malloc(nsize*sizeof(uint));
      double *EngdispInsize =(double *)malloc(nsize*sizeof(double));
	  uint *contactEndInsize =(uint *)malloc(nsize*sizeof(uint));
	  memset(NumInsize,0,nsize*sizeof(uint));
	  memset(EngdispInsize,0.0,nsize*sizeof(double));
	  memset(contactEndInsize,0,nsize*sizeof(uint));
//--------------------------------------------------------------
	  int nsect=(int) ((m_meshOrigin[2][2]-m_meshOrigin[1][2])*diam/0.5)+1;
	  uint *NumInsect =(uint *)malloc(nsect*sizeof(uint));
      double *EngdispInsect =(double *)malloc(nsect*sizeof(double));
	  double *MeanxInsect =(double *)malloc(nsect*sizeof(double));
	  uint *contactEndInsect =(uint *)malloc(nsect*sizeof(uint));

	  memset(NumInsect,0,nsect*sizeof(uint));
	  memset(EngdispInsect,0.0,nsect*sizeof(double));
	  memset(MeanxInsect,0.0,nsect*sizeof(double));
	  memset(contactEndInsect,0,nsect*sizeof(uint));

	  sprintf(str1,"%09.6f",tsd);

	   if(m_numParticles_total>0)
       {
       for(i=0;i<m_numParticles_total;i++)
        {
        px=hpos_total[3*i];
        py=hpos_total[3*i+1];
        pz=hpos_total[3*i+2];
		double velmag=sqrt(hVel_total[3*i]*hVel_total[3*i]+hVel_total[3*i+1]*hVel_total[3*i+1]+hVel_total[3*i+2]*hVel_total[3*i+2]);
        velmag=velmag*sqrt(gg*diam)/dt; 
		hEngdisp_total[i*3] *= denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0; //unit:J
		hEngdisp_total[i*3+1] *= denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0;
		hEngdisp_total[i*3+2] *= denp*pi*pow(2.0*hrad_total[i]*diam,3)/6.0;

		if(velmag<=1.0 && pz<m_meshEnd[1][2] && px>m_meshOrigin[1][0] && hEngdisp_total[i*3]>1e-20) //settle down
		{
		icount++;
		totalEdisp += hEngdisp_total[i*3];
		totalcontactEnd +=hcontactEnd_total[3*i];
		totalcontactEndii +=hcontactEnd_total[3*i+1];
		totalcontactEndiw +=hcontactEnd_total[3*i+2];

		for(k=0;k<nsize;k++)
		{
			if(fabs(hrad_total[i]-Size[k])<1e-4)
			{
             NumInsize[k] +=1;
			 EngdispInsize[k] +=hEngdisp_total[i*3];
			 contactEndInsize[k] += hcontactEnd_total[3*i];
			 break;
			}
		}

		if(hEngdisp_total[i*3+1]>1e-20)
		{
		icountw++;
		totalEdispw += hEngdisp_total[i*3+1];
		}

		if(hEngdisp_total[i*3+2]>1e-20)
		{
		icountwi++;
		totalEdispwi += hEngdisp_total[i*3+2];
		}

		} //if(velmag<1e-2 
	    //---------------------------------------	
		if(velmag>1e-2 && px>xmin )//collision
		{
		int isect= (int) ((pz-m_meshOrigin[1][2])*diam/0.5);
		if(isect<nsect && isect>=0)
		{
		NumInsect[isect] +=1;
		EngdispInsect[isect] +=hEngdisp_total[i*3];
		MeanxInsect[isect] +=px;
	    contactEndInsect[isect] += hcontactEnd_total[3*i];
		}
		}
		//-------------------------
	    }

	   free(hpos_total);
	   free(hrad_total);

	   free(hVel_total);
	   free(hEngdisp_total);
	  // free(hBriProb_total);	
	   free(hcontactEnd_total);
	  } //if(m_numParticles_total>0)

//----------------------------------------------------------
	   double averEdisp=0.0;
	   double averEdispw=0.0;
	   double averEdispwi=0.0;
	   double avercontactEnd=0.0;
	   double avercontactEndii=0.0;
	   double avercontactEndiw=0.0;

       if(icount<=0)icount=1;
	   averEdisp=totalEdisp/icount;
	   avercontactEnd=1.0*totalcontactEnd/icount;
	   

	   if(icountw<=0)icountw=1;
	   averEdispw=totalEdispw/icountw;
       avercontactEndii=1.0*totalcontactEndii/icountw;

	   if(icountwi<=0)icountwi=1;
	   averEdispwi=totalEdispwi/icountwi;
       avercontactEndiw=1.0*totalcontactEndiw/icountwi;

      if((fp1=fopen("output/AverEngdisp.dat","ab+"))==NULL){
      printf("Cannot open file AverEngdisp.dat! %d\n",it);
      exit(0);
	  }
       fprintf(fp1,"%10.6g %8d %8d %8d %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
		tsd,icount,icountw,icountwi,avercontactEnd,avercontactEndii,avercontactEndiw,averEdisp/fmax(1,icount*avercontactEnd),\
		averEdispw/fmax(1,icountw*avercontactEndii),averEdispwi/fmax(1,icountwi*avercontactEndiw),averEdisp,averEdispw,averEdispwi);

       fclose(fp1);
//----------------------------------------------------------
      if((fp1=fopen("output/DisEngdisp.dat","ab+"))==NULL){
      printf("Cannot open file DisEngdisp.dat! %d\n",it);
      exit(0);
	  }

	  fprintf(fp1,"ZONE T=\"ZONEP%s\"\n",str1);

	   	for(k=0;k<nsize;k++)
		{
		if(NumInsize[k]<1)NumInsize[k]=1;
		EngdispInsize[k] /=NumInsize[k];
		contactEndInsize[k] /=NumInsize[k];
		fprintf(fp1,"%10.6f %8d %8d %15.6g\n",
		2.0*Size[k]*diam,NumInsize[k],contactEndInsize[k],EngdispInsize[k]);
		}
	   free(NumInsize);	
	   free(EngdispInsize);	

	   fclose(fp1);
//----------------------------------------------------------
      if((fp1=fopen("output/AxialEngdisp.dat","ab+"))==NULL){
      printf("Cannot open file AxialEngdisp.dat! %d\n",it);
      exit(0);
	  }

	  fprintf(fp1,"ZONE T=\"ZONEP%s\"\n",str1);

	   	for(k=0;k<nsect;k++)
		{
		if(NumInsect[k]<1)NumInsect[k]=1;
		EngdispInsect[k] /=NumInsect[k];
		MeanxInsect[k] /=NumInsect[k];
		contactEndInsect[k] /=NumInsect[k];
		double dh=(k+1)*0.5-0.25-(m_meshOrigin[2][2]-m_meshOrigin[1][2])*diam;
		if(dh>0.0)dh=0.0;
		fprintf(fp1,"%10.6f %8d %8d %15.6g %15.6g\n",
		dh,NumInsect[k],contactEndInsect[k],EngdispInsect[k],MeanxInsect[k]*diam);
		}
	   free(NumInsect);	
	   free(EngdispInsect);	

	   fclose(fp1);
	  }
	
}



void WriteOutContactFrequency()
{
    uint *m_hcontactSizepair_total;
     if(rank==0)m_hcontactSizepair_total=(uint *)malloc((6*Nsect)*sizeof(uint));

     cpuErrchk(cudaMemcpy(m_hcontactSizepair,m_dcontactSizepair,6*Nsect*sizeof(uint),cudaMemcpyDeviceToHost));
     MPI_Reduce(m_hcontactSizepair,m_hcontactSizepair_total,6*Nsect,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);

      if(rank==0)
      {
      if((fp1=fopen("output/ContactFrequency.dat","ab+"))==NULL){
      printf("Cannot open file ContactFrequency.dat! %d\n",it);
      exit(0);
	  }

      double tdrop=(it-itstarthmin+1)*real_dt;
	  fprintf(fp1,"ZONE T=\"ZONEW%06.3f\"\n",tsd);

	  //printf("tsd=%15.6f tdrop=%15.6f\n",tsd,tdrop);

	  for(i=0;i<Nsect;i++)
       {
		double Binstart=pow(10,i*gap+log10(Emin));
        fprintf(fp1,"%15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
		   Binstart,m_hcontactSizepair_total[0*Nsect+i]/tdrop,m_hcontactSizepair_total[1*Nsect+i]/tdrop,\
		   m_hcontactSizepair_total[2*Nsect+i]/tdrop,m_hcontactSizepair_total[3*Nsect+i]/tdrop,\
		   m_hcontactSizepair_total[4*Nsect+i]/tdrop,m_hcontactSizepair_total[5*Nsect+i]/tdrop);
	   }

       fclose(fp1);
	   free(m_hcontactSizepair_total);
	  }
}

void WriteTimesteps()
{
	 if(rank==0)
      {
      if((fp1=fopen("output/timesteps.dat","ab+"))==NULL)
      printf("Cannot open file timesteps.dat!");
      tsd=(it)*real_dt;
      fprintf(fp1,"%15.6g %15.6g %15.6g\n", tret,it*1.0,tsd);
      fclose(fp1);
	 }
}

void WriteParticleForRerun()
{
       sprintf(str1,"%02d",rank);
       sprintf(str,"preflow%s.dat",str1);

        if((fp1=fopen(str,"wb+"))==NULL)printf("Cannot open file %s!",str);
        fprintf(fp1,"%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n", 
			it,m_numParticles,m_totalcontacts,itime, \
			FeedId,m_htotalout,NallocateH,NallocateM,NallocateL,Feedstart,hminflag,itstarthmin);

		if(m_numParticles>0)
		{
  CopyParticleArrayDeviceToHost(m_hprehead, m_numParticles);
  CopyNeigListArrayDeviceToHost(m_totalcontacts);

        for(i=0;i<m_numParticles;i++){
        fprintf(fp1,"%15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %4d\n",
			m_hpos[3*i],  m_hpos[3*i+1],  m_hpos[3*i+2],
			m_hpdis[3*i], m_hpdis[3*i+1], m_hpdis[3*i+2],
            m_hangv[3*i], m_hangv[3*i+1], m_hangv[3*i+2],
            m_hrad[i], m_hrmass[i], m_hinert[i], m_hmatId[i]);
		}

        for(i=0;i<m_numParticles;i++){
         fprintf(fp1,"%8u %15.6g %15.6g %15.6g %15.6g\n",
			 m_holdsjgi[i],m_holddisptw[3*i],m_holddisptw[3*i+1],m_holddisptw[3*i+2],m_holdfricpw[i]);
		}

       for(i=0;i<m_numParticles;i++){
         fprintf(fp1,"%8u %8u %8u %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
			 m_hcontactEnd[3*i],m_hcontactEnd[3*i+1],m_hcontactEnd[3*i+2],m_hEngdisp[i*3],m_hEngdisp[i*3+1],m_hEngdisp[i*3+2],m_hEngdispw[i],\
			 m_hEngdispVarw[4*i],m_hEngdispVarw[4*i+1],m_hEngdispVarw[4*i+2],m_hEngdispVarw[4*i+3]);
		}

        for(i=0;i<m_totalcontacts;i++){
        fprintf(fp1,"%8u %15.6g %15.6g %15.6g %15.6g",
				m_holdLn[i],m_holddispt[3*i],m_holddispt[3*i+1],m_holddispt[3*i+2],m_holdfricp[i]);
        fprintf(fp1,"\n");
        }

		for(i=0;i<m_totalcontacts;i++){
        fprintf(fp1,"%8u %15.6g %15.6g %15.6g %15.6g %15.6g",\
         m_holdcontactEndpair[i],m_hEngdisppair[i],m_hEngdispVarpair[4*i],\
		 m_hEngdispVarpair[4*i+1],m_hEngdispVarpair[4*i+2],m_hEngdispVarpair[4*i+3]);
        fprintf(fp1,"\n");
        }

       cpuErrchk(cudaMemcpy(m_hcontactSizepair,m_dcontactSizepair,6*Nsect*sizeof(uint),
		 cudaMemcpyDeviceToHost));

		for(i=0;i<6*Nsect;i++)
			fprintf(fp1,"%8u ", m_hcontactSizepair[i]);

		fprintf(fp1,"\n");

		}
        fclose(fp1);
}
