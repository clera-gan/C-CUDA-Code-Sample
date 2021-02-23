#include "dempacking.h"
//#include "ReadData.h"

//    ******************************************************************
//    *	         read initial values to parameters	              *
//    ******************************************************************

void readinputparameters()
{
		printf("readinputparameters. rank=%4d\n",rank);
//-------------------------
      if((fp1=fopen("input/hopp3d.inp","rb"))==NULL){
      printf("Cannot open file hopp3d.inp !");
      exit(1);
	   }
//-----define the youngs modulus, major semi-axis, and density
       fscanf(fp1, "%d", &iop_fn);	
//-------------------------------------
       fscanf(fp1, "%d %lf %lf", &ntot_total,&denp,&p_jpolfac);
       fscanf(fp1," %lf %lf", &dt, &tstop);
       fscanf(fp1," %lf %lf", &vel0,&tnewpar);
       fscanf(fp1," %le %le %d",&hamaker,&cutoffdis,&icalculateFvdw);
       fscanf(fp1," %lf %lf %d\n",&bondforceratio,&bondmaxgap,&icalculateFbond);

	   fscanf(fp1," %lf %lf %lf %lf %lf %lf",&m_houtmin[0],&m_houtmin[1],&m_houtmin[2],&m_houtmax[0],&m_houtmax[1],&m_houtmax[2]);
       fscanf(fp1," %d %d %d %d",&devnum,&maxCnPerParticle,&maxFacePerGrid,&AddOnceForAll);
       fscanf(fp1," %d %d %d",&DivX,&DivY,&DivZ);
       fscanf(fp1," %d %d %d %d %d",&freqppor,&freqparticle,&freqcontfreq,&freqftv,&freqpreflow);
        fclose(fp1);
//---------output data for check----------------------------
		if(rank==0)
		{
       printf("%d\n", iop_fn);	
//-------------------------------------
       printf("%d %lf %lf\n", ntot_total,denp,p_jpolfac);
       printf("%lf %lf\n", dt, tstop);
       printf("%lf %lf\n", vel0,tnewpar);
       printf("%le %le %d\n",hamaker,cutoffdis,icalculateFvdw);
       printf("%lf %lf %d\n",bondforceratio,bondmaxgap,icalculateFvdw);
       printf("%lf %lf %lf %lf %lf %lf",m_houtmin[0],m_houtmin[1],m_houtmin[2],m_houtmax[0],m_houtmax[1],m_houtmax[2]);
       printf("%d  %d %d %d\n",devnum,maxCnPerParticle,maxFacePerGrid,AddOnceForAll);
	   printf("%d %d %d\n",DivX,DivY,DivZ);
       printf("%d %d %d %d %d\n",freqppor,freqparticle,freqcontfreq,freqftv,freqpreflow);
		}
		
		printf("readinputparameters completed. rank=%4d\n",rank);
//-------------------------------------------------------
}

 void readmaterialproperties()
{
	printf("readmaterialproperties. rank=%4d\n",rank);
   if((fp1=fopen("input/materials.dat","rb"))==NULL){
       printf("Cannot open file materials.dat!\n");
       exit(1);}
    if(rank==0)printf("materials.dat:\n");

    fscanf(fp1,"%d", &NumOfMat);
    if(rank==0)printf("NumOfMat=%8d\n",NumOfMat);
	m_hmatpro =(Matproperties *)malloc(NumOfMat*sizeof(Matproperties));
	m_hmatSize =(MatSizeDist *)malloc(NumOfMat*sizeof(MatSizeDist));
     radmin=1.0;
    radmax=0.0;

    for(i=0;i<NumOfMat;i++)
    {
    fscanf(fp1,"%s %lf %lf %lf %lf %lf %lf %lf %d", 
		&m_hmatSize[i].name,&m_hmatpro[i].denp,&m_hmatpro[i].estarp,&m_hmatpro[i].ufric,&m_hmatpro[i].dampn,&m_hmatpro[i].dampt,&m_hmatpro[i].vp,&m_hmatpro[i].upp,&m_hmatSize[i].sizeNum);
    if(rank==0)printf("%s %lf %lf %lf %lf %lf %lf %lf %8d\n", 
		m_hmatSize[i].name,m_hmatpro[i].denp,m_hmatpro[i].estarp,m_hmatpro[i].ufric,m_hmatpro[i].dampn,m_hmatpro[i].dampt,m_hmatpro[i].vp,m_hmatpro[i].upp,m_hmatSize[i].sizeNum);

    for(j=0;j<m_hmatSize[i].sizeNum;j++)
    {
    fscanf(fp1,"%lf ",&m_hmatSize[i].sizedis[j].dia);
    if(rank==0)printf("mat %4d size %4d's diam=%f \n",i,j,m_hmatSize[i].sizedis[j].dia);

     radmin= fmin(radmin,m_hmatSize[i].sizedis[j].dia);
     radmax= fmax(radmax,m_hmatSize[i].sizedis[j].dia);
    }
     uint totalNum=0;
    for(j=0;j<m_hmatSize[i].sizeNum;j++)
    {
    fscanf(fp1,"%d ",&m_hmatSize[i].sizedis[j].Num);
    if(rank==0)printf("mat %4d size %4d's No.=%8d \n",i,j,m_hmatSize[i].sizedis[j].Num);
    totalNum+=m_hmatSize[i].sizedis[j].Num;
    }
     m_hmatSize[i].TotalAddNum=totalNum;

    } //for(i=0;i<NumOfMat;i++)

    diam=radmax;
    radmax /=diam;
    radmin /=diam;
    if(rank==0)printf("materials maxrad=%6.3f,minrad=%6.3f\n",radmax,radmin);

    for(i=0;i<NumOfMat;i++)
    {
    if(iop_fn ==1)
    {   
    m_hmatpro[i].estarp=6.0*m_hmatpro[i].estarp/(3.1415926*m_hmatpro[i].denp*9.81*diam);
    m_hmatpro[i].estarp=1.3333*m_hmatpro[i].estarp/(2.0*(1.0-m_hmatpro[i].vp*m_hmatpro[i].vp));
    }
     }

   fclose(fp1); 

    m_hParW.estarw=10.0*m_hmatpro[0].estarp;
    m_hParW.ufricw=m_hmatpro[0].ufric;
    m_hParW.dampnw=m_hmatpro[0].dampn;
    m_hParW.damptw=m_hmatpro[0].dampt;
    m_hParW.vw=m_hmatpro[0].vp;
    m_hParW.upw=m_hmatpro[0].upp;
    if(rank==0)printf("wall property: estarw  ufricw  dampnw  damptw  vw  upw:\n");
    if(rank==0)printf("%15.6f%15.6f%15.6f%15.6f%15.6f%15.6f:\n",
              m_hParW.estarw,m_hParW.ufricw,m_hParW.dampnw,m_hParW.damptw,m_hParW.vw,m_hParW.upw);

	printf("readmaterialproperties completed. rank=%4d\n",rank);
}

void readfeedinfo()
{
    char matname[20];
    if((fp1=fopen("input/feed.dat","rb"))==NULL){
       printf("Cannot open file feed.dat!\n");
       exit(1);}
    if(rank==0)printf("feed.dat:\n");

    fscanf(fp1,"%d", &NumOfFeed);
    if(rank==0)printf(" total feed batch number=%8d\n",NumOfFeed);
    
	m_hfeedpro=(Feedproperties *)malloc(NumOfFeed*sizeof(Feedproperties));
	m_hfeedlim=(Feedlimit *)malloc(NumOfFeed*sizeof(Feedlimit));

    for(i=0;i<NumOfFeed;i++)
    {
    fscanf(fp1,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", 
      &matname,&m_hfeedpro[i].xmin,&m_hfeedpro[i].xmax,&m_hfeedpro[i].ymin,&m_hfeedpro[i].ymax,&m_hfeedpro[i].zmin,&m_hfeedpro[i].zmax,&m_hfeedpro[i].starttime,&m_hfeedpro[i].endtime,&m_hfeedpro[i].freq,&m_hfeedpro[i].feedParNum);

	if(rank==0)
	printf("feed property %4d: %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
    i,matname,m_hfeedpro[i].xmin,m_hfeedpro[i].xmax,m_hfeedpro[i].ymin,m_hfeedpro[i].ymax,m_hfeedpro[i].zmin,m_hfeedpro[i].zmax,m_hfeedpro[i].starttime,m_hfeedpro[i].endtime,m_hfeedpro[i].freq,m_hfeedpro[i].feedParNum);

     m_hfeedpro[i].xmin /=diam;
     m_hfeedpro[i].xmax /=diam;
     m_hfeedpro[i].ymin /=diam;
     m_hfeedpro[i].ymax /=diam;
     m_hfeedpro[i].zmin /=diam;
     m_hfeedpro[i].zmax /=diam; 

    for(j=0;j<NumOfMat;j++) 
     {
     if(strcmp(matname, m_hmatSize[j].name)==0)
     m_hfeedpro[i].matId=j;
    if(rank==0) printf("feed %4d's matId is %4d\n",i,m_hfeedpro[i].matId);
     break;
     }
     if(i==NumOfFeed)m_hfeedpro[i].matId=0;

	 m_hfeedlim[i].matId=m_hfeedpro[i].matId;
	 m_hfeedlim[i].xmin =m_hfeedpro[i].xmin;
     m_hfeedlim[i].xmax =m_hfeedpro[i].xmax;
     m_hfeedlim[i].ymin =m_hfeedpro[i].ymin;
     m_hfeedlim[i].ymax =m_hfeedpro[i].ymax;
     m_hfeedlim[i].zmin =m_hfeedpro[i].zmin;
     m_hfeedlim[i].zmax =m_hfeedpro[i].zmax;

	if(rank==0)
	printf("feed limit reduced %4d: %lf %lf %lf %lf %lf %lf \n",
    i,m_hfeedlim[i].xmin,m_hfeedlim[i].xmax,m_hfeedlim[i].ymin,m_hfeedlim[i].ymax,m_hfeedlim[i].zmin,m_hfeedlim[i].zmax);
    }
    fclose(fp1);
}

void readmeshinfo()
{
	pxmax=0.0;
	pymax=0.0;
	pzmax=0.0;

	pxmin=1000.0;
	pymin=1000.0;
	pzmin=1000.0;

       if((fp=fopen("input/meshinfo.dat","rb"))==NULL){
      printf("Cannot open file meshinfo.dat!\n");
	  exit(1);}
      fscanf(fp,"%d", &Nmesh);
      if(rank==0)printf("Nmesh=%d\n", Nmesh);

      m_hmeshFile = (MeshFile *)malloc(Nmesh*sizeof(MeshFile));
      Nnode=(uint *)malloc(Nmesh*sizeof(uint));
      Nface=(uint *)malloc(Nmesh*sizeof(uint));

      m_meshOrigin=(double **)malloc(Nmesh*sizeof(double *));
      m_meshEnd=(double **)malloc(Nmesh*sizeof(double *));

      m_meshSize=(uint **)malloc(Nmesh*sizeof(uint *));
      m_nmeshCells=(uint *)malloc(Nmesh*sizeof(uint));

	  m_hPosnode=(double **)malloc(Nmesh*3*sizeof(double *));
       if(m_hPosnode==NULL)exit(1);
	  m_hFnodeid=(uint **)malloc(Nmesh*4*sizeof(uint *));
       if(m_hFnodeid==NULL)exit(1);

      for(i=0;i<Nmesh;i++)
      {
        m_meshOrigin[i]=(double *)malloc(3*sizeof(double));
       m_meshEnd[i]=(double *)malloc(3*sizeof(double));
       m_meshSize[i]=(uint *)malloc(3*sizeof(uint));

       m_meshOrigin[i][0]=10000.0;
       m_meshOrigin[i][1]=10000.0;
       m_meshOrigin[i][2]=10000.0;
       m_meshEnd[i][0]=-10000.0;
       m_meshEnd[i][1]=-10000.0;
       m_meshEnd[i][2]=-10000.0;

       //char meshfilename[25];
	   //double mx=0,my=0,mz=0;

       fscanf(fp,"%s %lf  %lf  %lf  %lf",&m_hmeshFile[i].name,&m_hmeshFile[i].disp.x,&m_hmeshFile[i].disp.y,&m_hmeshFile[i].disp.z,&m_hmeshFile[i].scalefactor);
       if(rank==0)
		   printf("mesh %d name: %s.dat, disp= %lf  %lf  %lf, scalefactor= %lf\n",
		   i,m_hmeshFile[i].name,m_hmeshFile[i].disp.x,m_hmeshFile[i].disp.y,m_hmeshFile[i].disp.z,m_hmeshFile[i].scalefactor);
 
	   sprintf(str,"input/meshFiles/%s",m_hmeshFile[i].name);

       if((fp1=fopen(str,"rb"))==NULL){
       printf("Cannot open file %s.dat!\n",m_hmeshFile[i].name);
	   exit(1);}


       fscanf(fp1,"%d %d",&Nnode[i], &Nface[i]);
       if(rank==0)printf("mesh %4d,node=%8d,face=%8d\n",i,Nnode[i],Nface[i]);

	m_hPosnode[i]=(double *)malloc(Nnode[i]*3*sizeof(double));
       if(m_hPosnode[i]==NULL)exit(1);
	m_hFnodeid[i]=(uint *)malloc(Nface[i]*4*sizeof(uint));
       if(m_hFnodeid[i]==NULL)exit(1);

        //node
	 for(k=0;k<Nnode[i];k++)
        {
          fscanf(fp1,"%lf %lf %lf",
           &m_hPosnode[i][k*3],&m_hPosnode[i][k*3+1],&m_hPosnode[i][k*3+2]);
		   m_hPosnode[i][k*3]   +=m_hmeshFile[i].disp.x;
           m_hPosnode[i][k*3+1] +=m_hmeshFile[i].disp.y;
           m_hPosnode[i][k*3+2] +=m_hmeshFile[i].disp.z;

           m_hPosnode[i][k*3]   *=m_hmeshFile[i].scalefactor/diam;
           m_hPosnode[i][k*3+1] *=m_hmeshFile[i].scalefactor/diam;
           m_hPosnode[i][k*3+2] *=m_hmeshFile[i].scalefactor/diam;

	     m_meshOrigin[i][0]=fmin(m_meshOrigin[i][0],m_hPosnode[i][k*3]);
	     m_meshOrigin[i][1]=fmin(m_meshOrigin[i][1],m_hPosnode[i][k*3+1]);
	     m_meshOrigin[i][2]=fmin(m_meshOrigin[i][2],m_hPosnode[i][k*3+2]);

	     m_meshEnd[i][0]=fmax(m_meshEnd[i][0],m_hPosnode[i][k*3]);
	     m_meshEnd[i][1]=fmax(m_meshEnd[i][1],m_hPosnode[i][k*3+1]);
	     m_meshEnd[i][2]=fmax(m_meshEnd[i][2],m_hPosnode[i][k*3+2]);

		pxmax=fmax(pxmax,m_hPosnode[i][k*3]);
		pymax=fmax(pymax,m_hPosnode[i][k*3+1]);
		pzmax=fmax(pzmax,m_hPosnode[i][k*3+2]);

		pxmin=fmin(pxmin,m_hPosnode[i][k*3]);
		pymin=fmin(pymin,m_hPosnode[i][k*3+1]);
		pzmin=fmin(pzmin,m_hPosnode[i][k*3+2]);

        }

       m_meshOrigin[i][0] -=0.5;
       m_meshOrigin[i][1] -=0.5;
       m_meshOrigin[i][2] -=0.5;
       m_meshEnd[i][0] +=0.5;
       m_meshEnd[i][1] +=0.5;
       m_meshEnd[i][2] +=0.5;

       //face
       for(k=0;k<Nface[i];k++)
        { 
          if(i==0)
          fscanf(fp1,"%d %d %d %d",
          &m_hFnodeid[i][k*4],&m_hFnodeid[i][k*4+1],&m_hFnodeid[i][k*4+2],&m_hFnodeid[i][k*4+3]);
          else
          fscanf(fp1,"%d %d %d %d",
          &m_hFnodeid[i][k*4+1],&m_hFnodeid[i][k*4],&m_hFnodeid[i][k*4+3],&m_hFnodeid[i][k*4+2]);

          // add to number
          m_hFnodeid[i][k*4]   -=1;
          m_hFnodeid[i][k*4+1] -=1;
          m_hFnodeid[i][k*4+2] -=1;
          m_hFnodeid[i][k*4+3] -=1; 

        } // for(i=0;i<

       fclose(fp1);
	   
//write mesh to file for checking ---------------
      sprintf(str,"%s.vtk",m_hmeshFile[i].name);
	  if(rank==0)
	  {
       if((fp1=fopen(str,"wb"))==NULL){
      printf("Cannot open file %s!\n",str);
      exit(0);
	   }
       fprintf(fp1,"# vtk DataFile Version 2.0\n");
       fprintf(fp1,"Mesh geometry\n");
       fprintf(fp1,"ASCII\n");
       fprintf(fp1,"DATASET UNSTRUCTURED_GRID\n");
       fprintf(fp1,"POINTS %d double\n",Nnode[i]);
       for(j=0;j<Nnode[i];j++)
       { 
       fprintf(fp1,"%15.6f %15.6f %15.6f\n",
       m_hPosnode[i][j*3]*diam,m_hPosnode[i][j*3+1]*diam,m_hPosnode[i][j*3+2]*diam);
       }
       fprintf(fp1,"CELLS %d %d\n",Nface[i],Nface[i]*5);
       for(j=0;j<Nface[i];j++)
       {
       fprintf(fp1,"%d %d %d %d %d\n",
       4,m_hFnodeid[i][j*4],m_hFnodeid[i][j*4+1],m_hFnodeid[i][j*4+2],m_hFnodeid[i][j*4+3]);
       }
      fprintf(fp1,"CELL_TYPES %d\n",Nface[i]);
       for(j=0;j<Nface[i];j++)
       {
       fprintf(fp1,"%d\n",9);
       }
      fclose(fp1);
	   }
       } // for(i=0;i<Nmesh;i++)
       fclose(fp);

	 // used for energy dissipation calculation
	   hmin=m_meshOrigin[1][2]-0.5;
	   hmax=m_meshOrigin[2][2]+1.5;
	   xmin=m_meshEnd[2][0];
	   if(rank==0)printf("%10.6f %10.6f %10.6f\n",hmin,hmax,xmin);
}

void readmovementinfo()
{
    if((fp1=fopen("input/movements.dat","rb"))==NULL){
       printf("Cannot open file movements.dat!\n");
       exit(1);}
    if(rank==0)printf("movements.dat:\n");

    fscanf(fp1,"%d %d %d ", &NumOfTranslation,&NumOfRotation,&NumOfVibration);
    if(rank==0)printf("NumOfTranslation=%8d,NumOfRotation=%8d,NumOfVibration=%8d\n",NumOfTranslation,NumOfRotation,NumOfVibration);

    if(NumOfTranslation>0)
    {
    m_hTranslation=(Translation *)malloc(NumOfTranslation*sizeof(Translation));

    for(i=0;i<NumOfTranslation;i++)
    {
    fscanf(fp1,"%d %lf %lf %lf %lf %lf",&m_hTranslation[i].meshId, &m_hTranslation[i].starttime,&m_hTranslation[i].endtime,&m_hTranslation[i].speed.x,&m_hTranslation[i].speed.y,&m_hTranslation[i].speed.z);
    if(rank==0)printf("%d %lf %lf %lf %lf %lf\n",
          m_hTranslation[i].meshId, m_hTranslation[i].starttime,m_hTranslation[i].endtime,m_hTranslation[i].speed.x,m_hTranslation[i].speed.y,m_hTranslation[i].speed.z);

	 m_hTranslation[i].speed.x *= dt/sqrt(gg*diam); //here speed is reduced to pdis
     m_hTranslation[i].speed.y *= dt/sqrt(gg*diam);
     m_hTranslation[i].speed.z *= dt/sqrt(gg*diam);
     printf("Reduced Translation[%4d].speed.x=%15.6f %15.6f %15.6f\n",i,m_hTranslation[i].speed.x ,m_hTranslation[i].speed.y,m_hTranslation[i].speed.z);
 
    }
    }

    if(NumOfRotation>0)
    {
    m_hRotation=(Rotation *)malloc(NumOfRotation*sizeof(Rotation));

    for(i=0;i<NumOfRotation;i++)
    {
    fscanf(fp1,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&m_hRotation[i].meshId,&m_hRotation[i].starttime,&m_hRotation[i].endtime,&m_hRotation[i].axis.x,&m_hRotation[i].axis.y,&m_hRotation[i].axis.z,&m_hRotation[i].pointonaxis.x,&m_hRotation[i].pointonaxis.y,&m_hRotation[i].pointonaxis.z,&m_hRotation[i].rotatespeed);
     if(rank==0)printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
     m_hRotation[i].meshId,m_hRotation[i].starttime,m_hRotation[i].endtime,m_hRotation[i].axis.x,m_hRotation[i].axis.y,m_hRotation[i].axis.z,m_hRotation[i].pointonaxis.x,m_hRotation[i].pointonaxis.y,m_hRotation[i].pointonaxis.z,m_hRotation[i].rotatespeed);

    }
    }

    if(NumOfVibration>0)
    {
    m_hVibration=(Vibration *)malloc(NumOfVibration*sizeof(Vibration));

    for(i=0;i<NumOfVibration;i++)
    {
    fscanf(fp1,"%d %lf %lf %lf %lf %lf",&m_hVibration[i].meshId,&m_hVibration[i].starttime,&m_hVibration[i].endtime,&m_hVibration[i].amplitude,&m_hVibration[i].frequency,&m_hVibration[i].startangle);
    if(rank==0)printf("%d %lf %lf %lf %lf %lf\n",
     m_hVibration[i].meshId,m_hVibration[i].starttime,m_hVibration[i].endtime,m_hVibration[i].amplitude,m_hVibration[i].frequency,m_hVibration[i].startangle);

    }
    }

    fclose(fp1);
}


void ReadRestartData()
{
   if((fp1=fopen("input/restart.dat","rb+"))==NULL)
    {
    printf("Cannot open file restart.dat!");
    exit(0);
	}
    fscanf(fp1,"%d %d",&nrestart,&nparticledat);
	fclose(fp1);
}

void ReadParticleData()
{
//**********************************************************************
//    *    if the value of nrestart is zero, { read data from the   *
//    *     file of preflow.dat                                        *
//**********************************************************************
	 sprintf(str1,"%02d",rank);
     sprintf(str,"preflow%s.dat",str1);

     if((fp1=fopen(str,"rb+"))==NULL){
     printf("Cannot open file %s!\n",str);
     exit(0);
	 }
      fscanf(fp1,"%ld %ld %ld %d %d %d %d %d %d %d %d %d",
		  &it,&m_numParticles,&m_totalcontacts,&itime,&FeedId,&m_htotalout,&NallocateH,&NallocateM,&NallocateL,&Feedstart,&hminflag,&itstarthmin);

	  	if(m_numParticles>0)
		{
	  if((fp2=fopen("view.dat","ab+"))==NULL){
      printf("Cannot open file view.dat 1#! %d\n");
      exit(0);
	  }
      if (m_numParticles!=ntot)fprintf(fp2, "n is unequal to ntot !\n");
        fclose(fp2);

        for(i=0;i<m_numParticles;i++){
        fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
			&m_hpos[3*i],&m_hpos[3*i+1],&m_hpos[3*i+2],&m_hpdis[3*i],&m_hpdis[3*i+1],&m_hpdis[3*i+2],
			&m_hangv[3*i],&m_hangv[3*i+1],&m_hangv[3*i+2],&m_hrad[i],&m_hrmass[i],&m_hinert[i],&m_hmatId[i]);
		}

        for(i=0;i<m_numParticles;i++){
         fscanf(fp1,"%u %lf %lf %lf %lf",
			 &m_holdsjgi[i],&m_holddisptw[3*i],&m_holddisptw[3*i+1],&m_holddisptw[3*i+2],&m_holdfricpw[i]);
		}

       for(i=0;i<m_numParticles;i++){
         fscanf(fp1,"%u %u %u %lf %lf %lf %lf %lf %lf %lf %lf",
			 &m_hcontactEnd[3*i],&m_hcontactEnd[3*i+1],&m_hcontactEnd[3*i+2],&m_hEngdisp[i*3],&m_hEngdisp[i*3+1],&m_hEngdisp[i*3+2],&m_hEngdispw[i],\
			 &m_hEngdispVarw[4*i],&m_hEngdispVarw[4*i+1],&m_hEngdispVarw[4*i+2],&m_hEngdispVarw[4*i+3]);
		}

        for(i=0;i<m_totalcontacts;i++){
        fscanf(fp1,"%u %lf %lf %lf %lf",
				&m_holdLn[i],&m_holddispt[3*i],&m_holddispt[3*i+1],&m_holddispt[3*i+2],&m_holdfricp[i]);
        }

	for(i=0;i<m_totalcontacts;i++){
        fscanf(fp1,"%u %lf %lf %lf %lf %lf",\
         &m_holdcontactEndpair[i],&m_hEngdisppair[i],&m_hEngdispVarpair[4*i],\
		&m_hEngdispVarpair[4*i+1],&m_hEngdispVarpair[4*i+2],&m_hEngdispVarpair[4*i+3]);
        }

		for(i=0;i<6*Nsect;i++)
		fscanf(fp1,"%u", &m_hcontactSizepair[i]);
        cpuErrchk(cudaMemcpy(m_dcontactSizepair,m_hcontactSizepair,6*Nsect*sizeof(uint),cudaMemcpyHostToDevice));
	   }
       fclose(fp1);
		
}
