#include "dempacking.h"
//#include "Feed.h"

int IsTimeToFeed()
{
 for(i=0;i<NumOfFeed;i++)
 {
  if(it> (int)(m_hfeedpro[i].starttime/real_dt) && it<= (int)(m_hfeedpro[i].endtime/real_dt) )
 return i;
 }
 return NumOfFeed;
}


void TreatSizeDistrbutionCN()
{
    k=0;

    for(i=0;i<NumOfFeed;i++)
    {
     j=m_hfeedpro[i].matId;
   
    for(uint m=0;m<m_hmatSize[j].sizeNum;m++)
    {

     double hrad=0.5*m_hmatSize[j].sizedis[m].dia/diam;

    if(k>0)
    {
    for(ir=0;ir<k;ir++)
    {
    if(Size[ir]==hrad)break;
    }

    if(ir==k) 
     {
     Size[k++]=hrad;
     radmin=fmin(radmin,hrad);
     }

    }
    else if(k==0)
     {
    Size[k++]=hrad;
     radmin=fmin(radmin,hrad);
     }
     } //for(m=0;m<m_hmatpro[j].sizeNum;m++)

     } // for(i=0;i<NumOfFeed;i++)

     nsize=k;

//---- size distribution -----------------
     m_hmaxCn =(uint *)malloc(nsize*sizeof(uint));

     treatSizedistribution(nsize,Size,radmin,m_cellSizeh,m_cellSizel,
                           m_hmaxCn,&radhmin,&radhmax,&radlmin,&radlmax);//return
      if(radlmin>radlmax)
      {
      radlmin=radlmax;
      m_hUSE_HIERARCHY=0;
      }
      else
      m_hUSE_HIERARCHY=1;

    /* m_hRlim.radhmin=radhmin;
     m_hRlim.radlmin=radlmin;
     m_hRlim.radhmax=radhmax;
     m_hRlim.radlmax=radlmax;*/
     m_hRlim.radmin=radmin;

}

 void calFeedVolFraction(double Dxmin, double Dymin, double Dzmin,
                          double Dxmax, double Dymax, double Dzmax,
                          double Fxmin, double Fymin, double Fzmin,
                          double Fxmax, double Fymax, double Fzmax,
                          double *deltav, double *deltax, double *deltay, double *deltaz,
                          double *FVfeed,
                          double *FeedOrigin,
                          double *FeedEnd)
{

    double FeedVolume;
    *deltav=0;
    *deltax=0;
    *deltay=0;
    *deltaz=0;
    *FVfeed=0;

     FeedVolume=(Fxmax -Fxmin)*(Fymax -Fymin)*(Fzmax -Fzmin);

     printf("Dxmin=%6.3f %6.3f %6.3f Dxmax=%6.3f %6.3f %6.3f Fxmin=%6.3f %6.3f %6.3f Fxmax=%6.3f %6.3f %6.3f\n",
             Dxmin,Dymin,Dzmin,Dxmax,Dymax,Dzmax,Fxmin,Fymin,Fzmin,Fxmax,Fymax,Fzmax);

    if(Dxmin<=Fxmin  && Fxmin<=Dxmax && Fxmax >Dxmax )
    { 
     *deltax=Dxmax -Fxmin;
     Fxmax=Dxmax;  
    }
    if( Fxmin< Dxmin && Dxmin <= Fxmax && Fxmax<=Dxmax  ) 
    { 
    *deltax=Fxmax -Dxmin;
    Fxmin= Dxmin;
    }

    if( Dxmin<=Fxmin  &&  Fxmax <=Dxmax) *deltax=Fxmax -Fxmin;

    if( Dymin<=Fymin  && Fymin<=Dymax && Fymax >Dymax) 
    { 
     *deltay=Dymax -Fymin;
     Fymax=Dymax;
    }
    if( Fymin< Dymin && Dymin <= Fymax && Fymax<=Dymax) 
    { 
    *deltay=Fymax -Dymin;
    Fymin=Dymin;
    }
    if( Dymin<=Fymin  &&  Fymax <=Dymax) *deltay=Fymax -Fymin; 

    if( Dzmin<=Fzmin  && Fzmin<=Dzmax && Fzmax >Dzmax )
    {
     *deltaz=Dzmax -Fzmin;
     Fzmax=Dzmax;
    }
    if( Fzmin< Dzmin && Dzmin <= Fzmax && Fzmax<=Dzmax  )
    { 
    *deltaz=Fzmax -Dzmin;
    Fzmin=Dzmin;
    }
    if( Dzmin<=Fzmin  &&  Fzmax <=Dzmax) *deltaz=Fzmax -Fzmin; 
     *deltav=(*deltax)*(*deltay)*(*deltaz);

    *FVfeed=*deltav/FeedVolume;

    FeedOrigin[0]=Fxmin;
    FeedOrigin[1]=Fymin;
    FeedOrigin[2]=Fzmin;

    FeedEnd[0]=Fxmin+*deltax;
    FeedEnd[1]=Fymin+*deltay;
    FeedEnd[2]=Fzmin+*deltaz;
}

void WhetherFeedinDomain()
{
    Feedindomain=0;
    FVfeed=0;
	deltav=0,deltax=0,deltay=0,deltaz=0;   

	for(i=0;i<NumOfFeed;i++)
    {
    calFeedVolFraction(Dxmin,Dymin,Dzmin,Dxmax,Dymax,Dzmax,m_hfeedpro[i].xmin,m_hfeedpro[i].ymin,m_hfeedpro[i].zmin,
                       m_hfeedpro[i].xmax,m_hfeedpro[i].ymax,m_hfeedpro[i].zmax,&deltav,&deltax,&deltay,&deltaz,&FVfeed,FeedOrigin,FeedEnd);

    if(FVfeed!=0)
    {
    Feedindomain=1;
	 m_hfeedlim[i].xmin = FeedOrigin[0];
     m_hfeedlim[i].xmax = FeedEnd[0];
     m_hfeedlim[i].ymin = FeedOrigin[1];
     m_hfeedlim[i].ymax = FeedEnd[1];
     m_hfeedlim[i].zmin = FeedOrigin[2];
     m_hfeedlim[i].zmax = FeedEnd[2];

	//nrate_coke= (int)(disrate_coke*FVfeed)+1; //reduced number of particles generated in a time step,num/s
    printf("rank=%4d,FVfeed=%6.3f,FeedOrigin=%6.3f %6.3f %6.3f,FeedEnd=%6.3f %6.3f %6.3f\n",
            rank,FVfeed,FeedOrigin[0],FeedOrigin[1],FeedOrigin[2],FeedEnd[0],FeedEnd[1],FeedEnd[2]);
	break;
	}

	/*
	MPI_Exscan(&nrate_coke,&nrate_coke_total,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);

    if(nrate_coke!=0) printf("rank=%4d,FVfeed=%6.3f,nrate_coke=%8d,nrate_ore=%8d,nrate_coke_total=%8d\n",
    rank,FVfeed,nrate_coke,nrate_ore,nrate_coke_total);

   if(nrate_coke_total+nrate_coke>disrate_coke) 
   { 
   if(nrate_coke>0) nrate_coke =disrate_coke-nrate_coke_total;
   }*/

	}
}

void SetFeedGridNumber()
{

   FeedgridSizeH=(uint **)malloc(NumOfFeed*sizeof(uint *));
   FeedgridSizeM=(uint **)malloc(NumOfFeed*sizeof(uint *));
   FeedgridSizeL=(uint **)malloc(NumOfFeed*sizeof(uint *));

   m_numgridsH =(uint *)malloc(NumOfFeed*sizeof(uint ));
   m_numgridsM =(uint *)malloc(NumOfFeed*sizeof(uint ));
   m_numgridsL =(uint *)malloc(NumOfFeed*sizeof(uint ));


    printf("rank=%4d,NumOfFeed=%4d,m_radHmax[0]=%6.4f,m_radLmax[0]=%6.4f\n",rank,NumOfFeed,m_radHmax[0],m_radLmax[0]);

   for(i=0;i<NumOfFeed;i++)
     {
        FeedgridSizeH[i]=(uint *)malloc(3*sizeof(uint));

         FeedgridSizeH[i][0]= (uint)((m_hfeedlim[i].xmax-m_hfeedlim[i].xmin)/(2.05*m_radHmax[i]));
          FeedgridSizeH[i][1]= (uint)((m_hfeedlim[i].ymax-m_hfeedlim[i].ymin)/(2.05*m_radHmax[i]));
          FeedgridSizeH[i][2]= (uint)((m_hfeedlim[i].zmax-m_hfeedlim[i].zmin)/(2.05*m_radHmax[i]));
         /*  FeedgridSizeH[i][0]= (uint)((deltax)/(2.05*m_radHmax[i]));
          FeedgridSizeH[i][1]= (uint)((deltay)/(2.05*m_radHmax[i]));
          FeedgridSizeH[i][2]= (uint)((deltaz)/(2.05*m_radHmax[i]));*/

          m_numgridsH[i]=FeedgridSizeH[i][0]*FeedgridSizeH[i][1]*FeedgridSizeH[i][2];

          printf("rank=%4d,Feed %4d FeedgridSizeH=%4d %4d %4d,m_numgridsH=%6d,m_radHmax=%6.3f\n",
          rank,i,FeedgridSizeH[i][0],FeedgridSizeH[i][1],FeedgridSizeH[i][2],m_numgridsH[i],m_radHmax[i]);

           if(m_hUSE_HIERARCHY==1)
           {
		  FeedgridSizeM[i]=(uint *)malloc(3*sizeof(uint));
         FeedgridSizeM[i][0]= (uint)((m_hfeedlim[i].xmax-m_hfeedlim[i].xmin)/(2.05*m_radMmax[i]));
          FeedgridSizeM[i][1]= (uint)((m_hfeedlim[i].ymax-m_hfeedlim[i].ymin)/(2.05*m_radMmax[i]));
          FeedgridSizeM[i][2]= (uint)((m_hfeedlim[i].zmax-m_hfeedlim[i].zmin)/(2.05*m_radMmax[i]));
           /*FeedgridSizeM[i][0]= (uint)((deltax)/(2.05*m_radMmax[i]));
          FeedgridSizeM[i][1]= (uint)((deltay)/(2.05*m_radMmax[i]));
          FeedgridSizeM[i][2]= (uint)((deltaz)/(2.05*m_radMmax[i]));*/

          m_numgridsM[i]=FeedgridSizeM[i][0]*FeedgridSizeM[i][1]*FeedgridSizeM[i][2];
          printf("Feed %4d FeedgridSizeM=%4d %4d %4d,m_numgridsL=%6d,m_radMmax=%6.3f\n",
          i,FeedgridSizeM[i][0],FeedgridSizeM[i][1],FeedgridSizeM[i][2],m_numgridsM[i],m_radMmax[i]);

          FeedgridSizeL[i]=(uint *)malloc(3*sizeof(uint));
         FeedgridSizeL[i][0]= (uint)((m_hfeedlim[i].xmax-m_hfeedlim[i].xmin)/(2.05*m_radLmax[i]));
          FeedgridSizeL[i][1]= (uint)((m_hfeedlim[i].ymax-m_hfeedlim[i].ymin)/(2.05*m_radLmax[i]));
          FeedgridSizeL[i][2]= (uint)((m_hfeedlim[i].zmax-m_hfeedlim[i].zmin)/(2.05*m_radLmax[i]));
         /*  FeedgridSizeL[i][0]= (uint)((deltax)/(2.05*m_radLmax[i]));
          FeedgridSizeL[i][1]= (uint)((deltay)/(2.05*m_radLmax[i]));
          FeedgridSizeL[i][2]= (uint)((deltaz)/(2.05*m_radLmax[i]));*/

          m_numgridsL[i]=FeedgridSizeL[i][0]*FeedgridSizeL[i][1]*FeedgridSizeL[i][2];

          printf("Feed %4d FeedgridSizeL=%4d %4d %4d,m_numgridsL=%6d,m_radLmax=%6.3f\n",
          i,FeedgridSizeL[i][0],FeedgridSizeL[i][1],FeedgridSizeL[i][2],m_numgridsL[i],m_radLmax[i]);
           }
    }

}

void SetFeedMaxRad()
{
   FeedMaxRad=(double *)malloc(NumOfFeed*sizeof(double));

   for(i=0;i<NumOfFeed;i++)
     {
      k=m_hfeedpro[i].matId;
    FeedMaxRad[i]=MatMaxRad[k];

    printf("rank=%4d,Feed %4d's matId=%4d,maxrad=%6.3f\n",rank,i,k,FeedMaxRad[i]);
    }
}

void AllocateCPUParticleBatchData()
{
   m_hpdisbatch =(double **)malloc(NumOfFeed*sizeof(double *));
   m_hangvbatch =(double **)malloc(NumOfFeed*sizeof(double *));

   m_hradbatch =(double **)malloc(NumOfFeed*sizeof(double *));
   m_hrmassbatch =(double **)malloc(NumOfFeed*sizeof(double *));
   m_hinertbatch =(double **)malloc(NumOfFeed*sizeof(double *));
   m_hmatIdbatch =(uint **)malloc(NumOfFeed*sizeof(uint *));
}

void InitialFeedBatchArray()
{
   printf("rank=%4d,InitialFeedBatchArray,FVfeed=%6.3f\n",rank,FVfeed);
   AllocateCPUParticleBatchData();

   m_NumAddH =(uint **)malloc(NumOfFeed*sizeof(uint *));
   m_NumAddM =(uint **)malloc(NumOfFeed*sizeof(uint *));
   m_NumAddL =(uint **)malloc(NumOfFeed*sizeof(uint *));

   m_radHmax=(double *)malloc(NumOfFeed*sizeof(double));
   m_radMmax=(double *)malloc(NumOfFeed*sizeof(double));
   m_radLmax=(double *)malloc(NumOfFeed*sizeof(double));

   m_radHmin=(double *)malloc(NumOfFeed*sizeof(double));
   m_radLmin=(double *)malloc(NumOfFeed*sizeof(double));

    Nbatch=(uint *)malloc(NumOfFeed*sizeof(uint));
   AllocateNum=(uint *)malloc(NumOfFeed*sizeof(uint));
   Nallocatebatch=(uint *)malloc(NumOfFeed*sizeof(uint));

   uint GlobalAllocateNum,AllocateNum_total;

   for(i=0;i<NumOfFeed;i++)
   {
    uint m=m_hfeedpro[i].matId;

     Nbatch[i]=(int) ceil( (double) m_hfeedpro[i].freq/0.002);

     AllocateNum[i]=m_hfeedpro[i].feedParNum;

	if(FVfeed!=0)
    {
	uint AllocateNum_Fv=min(AllocateNum[i],(uint)(AllocateNum[i]*FVfeed)+1);	 
    AllocateNum_total=0;
  /*  if(FVfeed!=0) printf("before rank=%4d,FVfeed=%6.3f,AllocateNum[i]=%8d,AllocateNum_Fv=%8d,AllocateNum_total=%8d\n",
    rank,FVfeed,AllocateNum[i],AllocateNum_Fv,AllocateNum_total);*/

    MPI_Exscan(&AllocateNum_Fv,&AllocateNum_total,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);

    if(FVfeed!=0) printf("rank=%4d,FVfeed=%6.3f,AllocateNum[i]=%8d,AllocateNum_Fv=%8d,AllocateNum_total=%8d\n",
    rank,FVfeed,AllocateNum[i],AllocateNum_Fv,AllocateNum_total);

   if(AllocateNum_total+AllocateNum_Fv>AllocateNum[i]) 
   { 
   AllocateNum_Fv=AllocateNum[i]-AllocateNum_total;
   }

   GlobalAllocateNum=AllocateNum[i];
   AllocateNum[i]=AllocateNum_Fv;
	}
	else
	{
	AllocateNum_total=0;
    GlobalAllocateNum=AllocateNum[i];
	}
     Nallocatebatch[i]= (int) ceil((float)AllocateNum[i]/Nbatch[i]);

	printf("rank=%4d,FVfeed=%6.3f,AllocateNum[i]=%8d,GlobalAllocateNum=%8d,Nallocatebatch[i]=%8d\n",
    rank,FVfeed,AllocateNum[i],GlobalAllocateNum,Nallocatebatch[i]);

	double *m_hradglobal =(double *)malloc(GlobalAllocateNum*sizeof(double));

	ip=0;
	 while(ip<GlobalAllocateNum)
     {
    for(j=0;j<m_hmatSize[m].sizeNum;j++)
    {
     for(k=0;k<m_hmatSize[m].sizedis[j].Num;k++)
     { 
     m_hradglobal[ip]=0.5*m_hmatSize[m].sizedis[j].dia/diam;
     ip++;
     } //for(k=0;k<m_hmatpro[m].sizedis[j].Num;k++)
	}
	 }

    for(j=0;j<5;j++)
    {
    std::random_shuffle(m_hradglobal,m_hradglobal+GlobalAllocateNum);
    }

    m_hpdisbatch[i] =(double *)malloc(AllocateNum[i]*3*sizeof(double));
    m_hangvbatch[i] =(double *)malloc(AllocateNum[i]*3*sizeof(double));

    m_hradbatch[i] =(double *)malloc(AllocateNum[i]*sizeof(double));
    m_hrmassbatch[i] =(double *)malloc(AllocateNum[i]*sizeof(double));
    m_hinertbatch[i] =(double *)malloc(AllocateNum[i]*sizeof(double));
    m_hmatIdbatch[i] =(uint *)malloc(AllocateNum[i]*sizeof(uint));

    m_NumAddH[i]=(uint *)malloc(Nbatch[i]*sizeof(uint));
	m_NumAddM[i]=(uint *)malloc(Nbatch[i]*sizeof(uint));
    m_NumAddL[i]=(uint *)malloc(Nbatch[i]*sizeof(uint));

    memset(m_NumAddH[i],0,Nbatch[i]*sizeof(uint));
    memset(m_NumAddM[i],0,Nbatch[i]*sizeof(uint));
	memset(m_NumAddL[i],0,Nbatch[i]*sizeof(uint));

    ip=0;
     while(ip<GlobalAllocateNum)
     {
    for(j=0;j<m_hmatSize[m].sizeNum;j++)
    {
      if(ip==GlobalAllocateNum)break;

     for(k=0;k<m_hmatSize[m].sizedis[j].Num;k++)
     { 
     if(ip==GlobalAllocateNum)break;

	 if(ip>=AllocateNum_total)
	 {
     m_hradbatch[i][ip-AllocateNum_total]=m_hradglobal[ip];
	 }
       ip++;
     } //for(k=0;k<m_hmatpro[m].sizedis[j].Num;k++)

     }// for(j=0;j<m_hmatpro[m].sizeNum;j++)
     }// do while(ip<AllocateNum[i])

	free(m_hradglobal);


    for(j=0;j<5;j++)
    {
    std::random_shuffle(m_hradbatch[i],m_hradbatch[i]+AllocateNum[i]);
    }

	if(AddOnceForAll!=1)
    {
    for(uint ib=0;ib<Nbatch[i];ib++)
    {
    for(j=ib*Nallocatebatch[i];j<(ib+1)*Nallocatebatch[i];j++)
    {
    for(k=j+1;k<(ib+1)*Nallocatebatch[i];k++)
    {
    if(m_hradbatch[i][j]< m_hradbatch[i][k])
     {
     double temp=m_hradbatch[i][j];
     m_hradbatch[i][j]=m_hradbatch[i][k];
     m_hradbatch[i][k]=temp;
     }
    }
    } //for(j=ib*Nbatch
     } //for(uint ib=0
	}
     else
     {
    for(j=0;j<1;j++)
    {
    std::random_shuffle(m_hradbatch[i],m_hradbatch[i]+AllocateNum[i]);
    }
     }

    for(j=0;j<AllocateNum[i];j++)
    {
     m_hrmassbatch[i][j]=8*m_hradbatch[i][j]*m_hradbatch[i][j]*m_hradbatch[i][j]*m_hmatpro[m].denp/denp;
     m_hinertbatch[i][j]=0.4*m_hrmassbatch[i][j]*(m_hradbatch[i][j]*m_hradbatch[i][j])*1.0;

       rno= rand()/(double)(RAND_MAX);
       ang=pi*(rno*.33333+.33333);
	  // vel0 /=sqrt(9.81*diam);

       m_hpdisbatch[i][j*3+2]=-0.75*vel0*sin(ang)*dt;
	   // m_hpdisbatch[i][j*3]= 0.1*vel0*sin(ang)*dt;

       dxy=-vel0*cos(ang)*dt;
       rno= rand()/(double)(RAND_MAX);
       ang=2.0*pi*rno;

      // m_hpdisbatch[i][j*3]=fabs(0.7*dxy*cos(ang));
	   m_hpdisbatch[i][j*3]=0.1*dxy*cos(ang);
       m_hpdisbatch[i][j*3+1]=0.1*dxy*sin(ang);

	m_hangvbatch[i][3*j]=0.0;
	m_hangvbatch[i][3*j+1]=0.0;
	m_hangvbatch[i][3*j+2]=0.0;

	m_hmatIdbatch[i][j]=m;
     }

	short Hindex=0, Mindex=0,Lindex=0;

	 for(j=0;j<AllocateNum[i];j++)
     {
      ibatch=(uint) floor(1.0*j/Nallocatebatch[i]);

     if(m_hradbatch[i][j]<FeedMaxRad[i]*0.5 && Hindex==0)
     { m_NumAddH[i][ibatch] =j-ibatch*Nallocatebatch[i];
      Hindex=1;
     }
     else if(m_hradbatch[i][j]<FeedMaxRad[i]*0.25 && Mindex==0)
     { m_NumAddM[i][ibatch] =j-ibatch*Nallocatebatch[i]-m_NumAddH[i][ibatch];
	   m_NumAddL[i][ibatch] = min((ibatch+1)*Nallocatebatch[i],AllocateNum[i])-\
		                     ibatch*Nallocatebatch[i]-m_NumAddM[i][ibatch]-m_NumAddH[i][ibatch];
	  Mindex=1;
	  Lindex=1;
     }
     if(Hindex==1 && Mindex==1)
	 {
	 Hindex=0;
	 Mindex=0;
	 Lindex=0;
	 j=(ibatch+1)*Nallocatebatch[i]-1;
	 //if(rank==0 || rank==1)
	 //printf("rank=%4d,i= %8d,ibatch=%8d,m_NumAddH=%8d,m_NumAddM=%8d,m_NumAddL=%8d\n",
	//	 rank,i,ibatch,m_NumAddH[i][ibatch],m_NumAddM[i][ibatch],m_NumAddL[i][ibatch]);
	 continue;
	 }

	 if(j==AllocateNum[i]-1)
     { 
		 if(m_hradbatch[i][j]>FeedMaxRad[i]*0.5)
		 {
	   m_NumAddH[i][ibatch]   = AllocateNum[i]-ibatch*Nallocatebatch[i];
	   m_NumAddM[i][ibatch] = 0;
       m_NumAddL[i][ibatch] = 0;
		 }
		 else if(m_hradbatch[i][j]<FeedMaxRad[i]*0.5 && m_hradbatch[i][j]>FeedMaxRad[i]*0.25)
		 {
       m_NumAddM[i][ibatch] = AllocateNum[i]-ibatch*Nallocatebatch[i]-m_NumAddH[i][ibatch];
       m_NumAddL[i][ibatch] = 0;
		 }
		 else //if(m_hradbatch[i][j]<=FeedMaxRad[i]*0.25)
           m_NumAddL[i][ibatch] = AllocateNum[i]-ibatch*Nallocatebatch[i]- \
		                         m_NumAddH[i][ibatch] -m_NumAddM[i][ibatch] ;
     break;
     }
    } // for(j=0;j<AllocateNum[i];j++)

	m_radHmax[i]=FeedMaxRad[i]*0.5;
    m_radMmax[i]=FeedMaxRad[i]*0.25;
    m_radLmax[i]=radlmin;

    m_radHmin[i]=radhmax;
    m_radLmin[i]=radlmin;

    for(j=0;j<AllocateNum[i];j++)
    {
    if(m_hradbatch[i][j]>m_radHmax[i])
     m_radHmax[i]=m_hradbatch[i][j];
     else if(m_hradbatch[i][j]>m_radMmax[i] && m_hradbatch[i][j]<FeedMaxRad[i]*0.5 )
     m_radMmax[i]=m_hradbatch[i][j];
     else if(m_hradbatch[i][j]>m_radLmax[i] && m_hradbatch[i][j]<FeedMaxRad[i]*0.25)
     m_radLmax[i]=m_hradbatch[i][j];

	if(m_hradbatch[i][j]<m_radHmin[i] && m_hradbatch[i][j]>radlmax)
		m_radHmin[i]=m_hradbatch[i][j];

    }//for(j=0;j<ip;j++)

	if(rank==0)
     printf("feed %4d m_radHmax=%6.3f m_radMmax=%6.3f m_radLmax=%6.3f,m_radHmin=%6.3f,m_radLmin=%6.3f,Nbatch[i]=%4d,Nallocatebatch[i]=%4d,AllocateNum[i]=%8d\n",
              i,m_radHmax[i],m_radMmax[i],m_radLmax[i],m_radHmin[i],m_radLmin[i],Nbatch[i],Nallocatebatch[i],AllocateNum[i]);

   } //for(i=0;i<NumOfFeed;i++)
  // exit(1);
}

void SetRandomFeedArray()
{

   NrandomH =(uint **)malloc(NumOfFeed*sizeof(uint *));

   if(m_hUSE_HIERARCHY==1)
   {
   NrandomM =(uint **)malloc(NumOfFeed*sizeof(uint *));
   NrandomL =(uint **)malloc(NumOfFeed*sizeof(uint *));
   }

   for(i=0;i<NumOfFeed;i++)
   {
    NrandomH[i] =(uint *)malloc(m_numgridsH[i]*sizeof(uint));

    for(j=0;j<m_numgridsH[i];j++)
    {
    NrandomH[i][j]=j;
    }

     if(m_hUSE_HIERARCHY==1)
    {
     NrandomM[i] =(uint *)malloc(m_numgridsM[i]*sizeof(uint));
    for(j=0;j<m_numgridsM[i];j++)
    {
    NrandomM[i][j]=j;
    }

     NrandomL[i] =(uint *)malloc(m_numgridsL[i]*sizeof(uint));
    for(j=0;j<m_numgridsL[i];j++)
    { 
    NrandomL[i][j]=j;
    }
    }

   }
}

void SetFeedArrays()
{
SetFeedMaxRad();
TreatSizeDistrbutionCN();

WhetherFeedinDomain();//newly added

InitialFeedBatchArray();

SetFeedGridNumber();
SetRandomFeedArray();

}

void AddParticlesBatch()
{
	if(Feedindomain==1) 
    {
          NallocateH=m_NumAddH[FeedId][ibatch];
		  NallocateM=m_NumAddM[FeedId][ibatch];
          NallocateL=m_NumAddL[FeedId][ibatch];
          Nallocate=NallocateH+NallocateM+NallocateL;

         /* radfeedLmax=m_radLmax[FeedId];
          radfeedHmax=m_radHmax[FeedId];*/
          //radfeedLmax=m_radLmin[FeedId];
          //radfeedHmax=m_radHmin[FeedId];
  // printf("add pars,rank=%8d,it=%8d,FeedId=%4d,ibatch=%4d,NallocateH=%8d,NallocateM=%8d,NallocateL=%8d,Nallocatebatch=%8d,ntot_total=%8d\n",
   //       rank,it,FeedId,ibatch,NallocateH,NallocateM,NallocateL,Nallocatebatch[FeedId],ntot_total);

          if(m_numParticles+Nallocate>ntot_total) //exceed number
          {
          Nallocate =ntot_total-m_numParticles;

          if(Nallocate<NallocateH)
          { NallocateH=Nallocate;
		    NallocateM=0;
            NallocateL=0;
          }
          else if(Nallocate>=NallocateH && Nallocate<NallocateH+NallocateM)
          { 
            NallocateM=Nallocate-NallocateH;
			NallocateL=0;
          }
		  else if(Nallocate>=NallocateH+NallocateM)
		  {
	        NallocateL=Nallocate-NallocateH-NallocateM;	  
		  }

          }
          m_numParticles=m_numParticles+Nallocate;

   m_maxcontacts=m_numParticles*maxCnPerParticle;  
   m_maxIlist=(m_numParIn+Nallocate)*AverIlist;
   if(m_maxIlist<m_totalcontacts)m_maxIlist=m_totalcontacts;
   allocatenei=m_numParticles*maxCnPerParticle;

   //printf("after add pars,rank=%8d,it=%8d,FeedId=%4d,ibatch=%4d,NallocateH=%8d,NallocateM=%8d,NallocateL=%8d,Nallocatebatch=%8d\n",
    //      rank,it,FeedId,ibatch,NallocateH,NallocateM,NallocateL,Nallocatebatch[FeedId]);
	}

 
}

void SetDataForParticles()
{
          if(Feedindomain==1 && m_numParticles>0)
           { 
           if(DivX==1)
           {
           allocatenum=(uint) (m_numParticles*1.2*(m_gridSize[0]+2.0)/(m_gridSize[0]-2.0)); //modify here
           if(allocatenum<m_oldnumParticles+Nallocate)
           allocatenum=(uint) ((m_oldnumParticles+Nallocate)*1.2*(m_gridSize[0]+2.0)/(m_gridSize[0]-2.0));
           }
           else if(DivZ==1)
           {
           allocatenum=(uint) (m_numParticles*1.2*(m_gridSize[2]+2.0)/(m_gridSize[2]-2.0)); //modify here
           if(allocatenum<m_oldnumParticles)
           allocatenum=(uint) (m_oldnumParticles*1.2*(m_gridSize[2]+2.0)/(m_gridSize[2]-2.0));
          }
		   if(allocatenum==0)allocatenum=1;
		   if(numprocs==1)allocatenum=m_numParticles;

          if(allocatenum>ntot_total)allocatenum=ntot_total;
           m_maxIlist=(m_numParIn+Nallocate)*6;
           }
           else //if(Feedindomain==0)
           {
           if(m_oldnumParticles==0) allocatenum= (uint)(m_numParticles_total*1.0/numprocs)+1;
           else
           {
           if(DivX)
           allocatenum=(uint) (m_oldnumParticles*1.2*(m_gridSize[0]+2.0)/(m_gridSize[0]-2.0));
           else if(DivZ)
           allocatenum=(uint) (m_oldnumParticles*1.2*(m_gridSize[2]+2.0)/(m_gridSize[2]-2.0));
           }
		   if(allocatenum==0)allocatenum=1;
		  if(numprocs==1)allocatenum=m_numParticles;
          if(allocatenum>ntot_total)allocatenum=ntot_total;
           m_maxIlist=(allocatenum)*6;
           }

		  //printf("rank=%8d,it=%8d,allocatenum=%8d,m_numParticles=%8d\n",rank,it,allocatenum,m_numParticles);

          if(m_maxIlist> 15*m_numParticles && m_numParticles>0)m_maxIlist=15*m_numParticles; //we do not want to allocate too much memory at the beginning
		  else if(m_numParticles==0)m_maxIlist=10;

          if(m_maxIlist<m_totalcontacts && m_totalcontacts!= 0xffffffff) 
           {
          /* printf("m_maxIlist<m_totalcontacts,rank=%8d,it=%8d,m_maxIlist=%8d,m_totalcontacts=%8d\n",
           rank,it,m_maxIlist,m_totalcontacts);*/

           m_maxIlist=m_totalcontacts;
           }

	if(allocatenum<m_oldnumParticles+Nallocate) allocatenum=m_oldnumParticles+Nallocate;	   
   allocatenei=allocatenum*maxCnPerParticle;

   AllocateParticleDeviceArrays(allocatenum,allocatenei); //note oldIp and Oldpos array also include
   AllocateNeigListDeviceArrays(m_maxIlist);

  //printf("old number! rank=%4d,it=%8d,m_oldnumParticles=%8d,m_numParticles=%8d,NallocateH=%8d,NallocateM=%8d,NallocateL=%8d,ibatch=%8d,Nallocatebatch[FeedId]=%8d\n",
	//  rank,it,m_oldnumParticles,m_numParticles,NallocateH,NallocateM,NallocateL,ibatch,Nallocatebatch[FeedId]);

  CopyParticlePropertyHostToDevice(m_oldnumParticles);
  CopyNeigListArrayHostToDevice(m_maxIlist);

  if(Feedindomain==1)
  {
	  
   if(Feedstart+Nallocate<AllocateNum[FeedId])
   {
	//  printf("ibatch*Nallocatebatch[FeedId]+Nallocate<AllocateNum[FeedId],rank=%4d,it=%8d,FeedId=%4d,m_oldnumParticles=%8d,m_numParticles=%8d,ibatch=%8d,Nallocatebatch[FeedId]=%8d %8d\n",
	 // rank,it,FeedId,m_oldnumParticles,m_numParticles,ibatch,Nallocatebatch[FeedId],AllocateNum[FeedId]);
	  
   CopyParticleBatchHostToDevice(Nallocate,m_oldnumParticles,Feedstart);
   Feedstart +=Nallocate;
   }
   else
   {
	  //printf("later ibatch*Nallocatebatch[FeedId]+Nallocate<AllocateNum[FeedId],rank=%4d,it=%8d,m_oldnumParticles=%8d,m_numParticles=%8d,ibatch=%8d,Nallocatebatch[FeedId]=%8d %8d\n",
	 // rank,it,m_oldnumParticles,m_numParticles,ibatch,Nallocatebatch[FeedId],AllocateNum[FeedId]);
	  
   CopyParticleBatchHostToDevice(AllocateNum[FeedId]-Feedstart,m_oldnumParticles,Feedstart);
   CopyParticleBatchHostToDevice(Nallocate-(AllocateNum[FeedId]-Feedstart),\
	   m_oldnumParticles+(AllocateNum[FeedId]-Feedstart),0);
   Feedstart=Nallocate-(AllocateNum[FeedId]-Feedstart);
   }

   if(NallocateH>0)
   {
   cudaMalloc((void**)&NumParticleInHCell, m_numgridsH[FeedId]*sizeof(uint));
   if(m_oldnumParticles==0)cudaMemset(NumParticleInHCell, 0, m_numgridsH[FeedId]*sizeof(uint));

   calcParticleInFeedAreaH(it,m_oldnumParticles,
          m_numgridsH[FeedId],
          m_dpos,
          m_drad,
		  m_radHmax[FeedId], 
          FeedgridSizeH[FeedId],
          FeedId,
          NumParticleInHCell //return value
          );

     std::random_shuffle(NrandomH[FeedId],NrandomH[FeedId]+m_numgridsH[FeedId]);
     cpuErrchk(cudaMalloc( (void**)&m_dNrandomH, m_numgridsH[FeedId]*sizeof(uint)));
     cpuErrchk(cudaMemcpy(m_dNrandomH,NrandomH[FeedId],m_numgridsH[FeedId]*sizeof(uint),cudaMemcpyHostToDevice));
  
     RandomAddParticlesH( 
         m_numgridsH[FeedId],
          m_oldnumParticles,
          m_dNrandomH,
          NallocateH,
          FeedgridSizeH[FeedId],
          m_drad,
          m_radHmax[FeedId],
          FeedId,
          NumParticleInHCell, 

          m_dpos,
          m_dmatId,
          Addcount);//return value

  cudaFree(m_dNrandomH);

  cudaMemcpy(&m_hAddcount,Addcount, sizeof(uint),cudaMemcpyDeviceToHost);

  if(m_hAddcount<NallocateH)
  {
   printf("not enough H particles generated,rank=%8d,NallocateH=%8d,Addcount=%8d,it=%8d,ibatch=%8d\n",rank,NallocateH,m_hAddcount,it,ibatch);
   exit(1);
   }
}

if(NallocateM>0)
{
   cudaMalloc((void**)&NumParticleInMCell, m_numgridsM[FeedId]*sizeof(uint));
   if(m_oldnumParticles==0)cudaMemset(NumParticleInMCell, 0, m_numgridsM[FeedId]*sizeof(uint));
   
   calcParticleInFeedAreaL(it,m_oldnumParticles+NallocateH,
          m_numgridsM[FeedId],
          m_dpos,
          m_drad,
          m_radMmax[FeedId],
          FeedgridSizeM[FeedId],
          FeedId,
          NumParticleInMCell//return value
          );

        std::random_shuffle(NrandomM[FeedId],NrandomM[FeedId]+m_numgridsM[FeedId]);
        cpuErrchk(cudaMalloc( (void**)&m_dNrandomM, m_numgridsM[FeedId]*sizeof(uint)));
        cpuErrchk(cudaMemcpy(m_dNrandomM,NrandomM[FeedId],m_numgridsM[FeedId]*sizeof(uint),cudaMemcpyHostToDevice));

		RandomAddParticlesL( 
          m_numgridsM[FeedId],
          m_oldnumParticles+NallocateH,
          m_dNrandomM,
          NallocateM,

          FeedgridSizeM[FeedId],
          m_radMmax[FeedId],
          FeedId,
          NumParticleInMCell,
          m_dpos,
          m_dmatId,
          Addcount); //return value 
  cudaFree(m_dNrandomM);
  cudaMemcpy(&m_hAddcount,Addcount, sizeof(uint),cudaMemcpyDeviceToHost);

  if(m_hAddcount<NallocateM)
  {
   printf("not enough M particles generated,NallocateH=%8d,Addcount=%8d,it=%8d\n",NallocateM,m_hAddcount,it);
   exit(1);
   }
  }

if(NallocateL>0)
{
   cudaMalloc((void**)&NumParticleInLCell, m_numgridsL[FeedId]*sizeof(uint));
   if(m_oldnumParticles==0)cudaMemset(NumParticleInLCell, 0, m_numgridsL[FeedId]*sizeof(uint));
   
   calcParticleInFeedAreaL(it,m_oldnumParticles+NallocateH+NallocateM,
          m_numgridsL[FeedId],
          m_dpos,
          m_drad,
          m_radLmax[FeedId],
          FeedgridSizeL[FeedId],
          FeedId,
          NumParticleInLCell//return value
          );

        std::random_shuffle(NrandomL[FeedId],NrandomL[FeedId]+m_numgridsL[FeedId]);
        cpuErrchk(cudaMalloc( (void**)&m_dNrandomL, m_numgridsL[FeedId]*sizeof(uint)));
        cpuErrchk(cudaMemcpy(m_dNrandomL,NrandomL[FeedId],m_numgridsL[FeedId]*sizeof(uint),cudaMemcpyHostToDevice));

		RandomAddParticlesL( 
          m_numgridsL[FeedId],
          m_oldnumParticles+NallocateH+NallocateM,
          m_dNrandomL,
          NallocateL,

          FeedgridSizeL[FeedId],
		  m_radLmax[FeedId],
          FeedId,
          NumParticleInLCell,
          m_dpos,
          m_dmatId,
          Addcount); //return value 
  cudaFree(m_dNrandomL);
  cudaMemcpy(&m_hAddcount,Addcount, sizeof(uint),cudaMemcpyDeviceToHost);

  if(m_hAddcount<NallocateL)
  {
   printf("not enough L particles generated,NallocateH=%8d,Addcount=%8d,it=%8d\n",NallocateL,m_hAddcount,it);
   exit(1);
   }
  }
 } //  if(Feedindomain==1)

  
 /*  if(Feedindomain==1 && it==1 && rank==0)
  {
   cpuErrchk(cudaMemcpy(m_hpos,m_dpos,m_numParticles*3*sizeof(double),cudaMemcpyDeviceToHost));
   cpuErrchk(cudaMemcpy(m_hrad,m_drad, m_numParticles*sizeof(double),cudaMemcpyDeviceToHost));

   for(i=0;i<m_numParticles;i++)
   {
	if(m_hrad[i]<radmin)
   printf("rank=%4d,it=%8d,i=%8d,pos=%6.3f %6.3f %6.3f,rad=%6.3f,m_oldnumParticles=%8d\n",
           rank,it,i,m_hpos[i*3],m_hpos[i*3+1],m_hpos[i*3+2],m_hrad[i],m_oldnumParticles);
   }
   }*/

  ibatch++;
  if(ibatch==Nbatch[FeedId])ibatch=0;
}
