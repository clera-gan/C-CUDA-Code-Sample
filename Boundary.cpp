#include "dempacking.h"
//#include "Boundary.h"
#include <mpi.h>
void SetGeometryBoundary()
{
    m_cellSize[0]=1.55;
    m_cellSize[1]=1.55;
    m_cellSize[2]=1.55;

    m_cellSizeh[0]=1.55;
    m_cellSizeh[1]=1.55;
    m_cellSizeh[2]=1.55;

    m_cellSizel[0]=m_cellSize[0]/3.0;
    m_cellSizel[1]=m_cellSize[1]/3.0;
    m_cellSizel[2]=m_cellSize[2]/3.0;

    m_hGlim.posxmin= pxmin-0.5;
    m_hGlim.posxmax= pxmax+0.5;
    m_hGlim.posymin= pymin-0.5;
    m_hGlim.posymax= pymax+0.5;
    m_hGlim.poszmin= pzmin-0.5;
    m_hGlim.poszmax= pzmax+0.5;

	m_worldSize_global[0] = m_hGlim.posxmax-m_hGlim.posxmin;
    m_worldSize_global[1] = m_hGlim.posymax-m_hGlim.posymin;
    m_worldSize_global[2] = m_hGlim.poszmax-m_hGlim.poszmin;

    m_worldOrigin_global[0] = m_hGlim.posxmin;
    m_worldOrigin_global[1] = m_hGlim.posymin;
    m_worldOrigin_global[2] = m_hGlim.poszmin;

    if(rank==0)
    {
    printf("global domain limit: x direction:%6.3f,%6.3f,y direction:%6.3f,%6.3f,z direction:%6.3f,%6.3f\n",
    m_hGlim.posxmin,m_hGlim.posxmax,m_hGlim.posymin,m_hGlim.posymax,m_hGlim.poszmin,m_hGlim.poszmax);
	}
	if(DivX==1)
	{
    m_gridSize[0] = (int)((int)((m_hGlim.posxmax-m_hGlim.posxmin)/m_cellSize[0]+1)/numprocs)+1; //local value
   // m_gridSize[0] = (int)((m_hGlim.posxmax-m_hGlim.posxmin)/m_cellSize[0])+1;
    m_gridSize[1] = (int)((m_hGlim.posymax-m_hGlim.posymin)/m_cellSize[1])+1;
    m_gridSize[2] = (int)((m_hGlim.poszmax-m_hGlim.poszmin)/m_cellSize[2])+1;

	if((m_gridSize[0]-1)*numprocs+rank+1>(int)((m_hGlim.posxmax-m_hGlim.posxmin)/m_cellSize[0])+1)
    m_gridSize[0]-=1;
    if(m_gridSize[0]<4)printf("not enough grid x,rank=%4d,m_gridSize[0]=%4d\n",rank,m_gridSize[0]);

    m_worldSize[0] = m_gridSize[0]*m_cellSize[0];  // devide region in x direction into numprocs parts
    m_worldSize[1] = m_hGlim.posymax-m_hGlim.posymin;
    m_worldSize[2] = m_hGlim.poszmax-m_hGlim.poszmin;

	MPI_Exscan(&m_worldSize[0],&m_worldOrigin[0],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(rank==numprocs-1)m_worldSize[0]=(m_hGlim.posxmax-m_hGlim.posxmin)-m_worldOrigin[0];

	// m_worldOrigin[0] = m_hGlim.posxmin; 
    m_worldOrigin[0] += m_hGlim.posxmin; 
    m_worldOrigin[1] = m_hGlim.posymin;
    m_worldOrigin[2] = m_hGlim.poszmin;
 
	m_gridSize_global[0] = (int)(m_worldSize_global[0]/m_cellSize[0])+1+2;
    m_gridSize_global[1] = (int)(m_worldSize_global[1]/m_cellSize[1])+1;
    m_gridSize_global[2] = (int)(m_worldSize_global[2]/m_cellSize[2])+1;

	printf("rank=%4d,m_gridSize=%4d %4d %4d,m_worldSize=%7.4f %7.4f %7.4f,m_worldOrigin=%7.4f %7.4f %7.4f,m_gridSize_global[0]\n",
		    rank,m_gridSize[0],m_gridSize[1],m_gridSize[2],m_worldSize[0],m_worldSize[1],m_worldSize[2],\
			m_worldOrigin[0],m_worldOrigin[1],m_worldOrigin[2],m_gridSize_global[0],m_gridSize_global[0],m_gridSize_global[0]);
	}
	else if(DivZ==1)
	{
     m_gridSize[0] = (int)((m_hGlim.posxmax-m_hGlim.posxmin)/m_cellSize[0])+1;
    m_gridSize[1] = (int)((m_hGlim.posymax-m_hGlim.posymin)/m_cellSize[1])+1;
    m_gridSize[2] = (int)((int)((m_hGlim.poszmax-m_hGlim.poszmin)/m_cellSize[2])+1/numprocs)+1;

	if((m_gridSize[2]-1)*numprocs+rank+1>(int)((m_hGlim.poszmax-m_hGlim.poszmin)/m_cellSize[2])+1)
    m_gridSize[2]-=1;
    if(m_gridSize[2]<4)printf("not enough grid z,rank=%4d,m_gridSize[2]=%4d\n",rank,m_gridSize[2]);

	m_worldSize[0] = m_hGlim.posxmax-m_hGlim.posxmin;
    m_worldSize[1] = m_hGlim.posymax-m_hGlim.posymin;
    m_worldSize[2] = m_gridSize[2]*m_cellSize[2]; // devide region in z direction into numprocs parts

	 MPI_Exscan(&m_worldSize[2],&m_worldOrigin[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(rank==numprocs-1)m_worldSize[2]=(m_hGlim.poszmax-m_hGlim.poszmin)-m_worldOrigin[2];

	m_worldOrigin[0] = m_hGlim.posxmin; 
    //m_worldOrigin[0] += m_hGlim.posxmin; 
    m_worldOrigin[1] = m_hGlim.posymin;
    m_worldOrigin[2] += m_hGlim.poszmin;

	m_gridSize_global[0] = (int)(m_worldSize_global[0]/m_cellSize[0])+1;
    m_gridSize_global[1] = (int)(m_worldSize_global[1]/m_cellSize[1])+1;
    m_gridSize_global[2] = (int)(m_worldSize_global[2]/m_cellSize[2])+1+2;
	}

	Dxmin=m_worldOrigin[0];
    Dymin=m_worldOrigin[1];
    Dzmin=m_worldOrigin[2];

    Dxmax=m_worldOrigin[0]+m_worldSize[0];
    Dymax=m_worldOrigin[1]+m_worldSize[1];
    Dzmax=m_worldOrigin[2]+m_worldSize[2];
	printf("domain limit: rank=%4d,x direction:%6.3f,%6.3f,y direction:%6.3f,%6.3f,z direction:%6.3f,%6.3f\n",
    rank,Dxmin,Dxmax,Dymin,Dymax,Dzmin,Dzmax);

    m_gridStart=0;

    if(DivX)
    {
    MPI_Exscan(&m_gridSize[0],&m_gridStart,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
    m_gridStart *=m_gridSize[1]*m_gridSize[2];
    m_gridSize[0]+=2;
    }
    else if(DivZ)
    {
    MPI_Exscan(&m_gridSize[2],&m_gridStart,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
    m_gridStart *=m_gridSize[0]*m_gridSize[1];
    m_gridSize[2]+=2; 
    }

	m_gridSizeh[0] = m_gridSize[0];
    m_gridSizeh[1] = m_gridSize[1];
    m_gridSizeh[2] = m_gridSize[2];

    m_gridSizel[0] = m_gridSizeh[0]*3;
    m_gridSizel[1] = m_gridSizeh[1]*3;
    m_gridSizel[2] = m_gridSizeh[2]*3;

    m_gridSizel_global[0] = m_gridSize_global[0]*3;
    m_gridSizel_global[1] = m_gridSize_global[1]*3;
    m_gridSizel_global[2] = m_gridSize_global[2]*3;

    m_nGridCellsh = m_gridSizeh[0]*m_gridSizeh[1]*m_gridSizeh[2];
    m_nGridCellsl = m_gridSizel[0]*m_gridSizel[1]*m_gridSizel[2];
    m_nGridCells  = m_nGridCellsh+m_nGridCellsl;

//    m_nGridCells = m_gridSize[0]*m_gridSize[1]*m_gridSize[2];
 
    printf("pxmin=%6.3f,pymin=%6.3f,pzmin=%6.3f\n",pxmin,pymin,pzmin);
    printf("pxmax=%6.3f,pymax=%6.3f,pzmax=%6.3f\n",pxmax,pymax,pzmax);

    m_houtmin[0] /=diam;
    m_houtmax[0] /=diam;
    m_houtmin[1] /=diam;
    m_houtmax[1] /=diam;
    m_houtmin[2] /=diam;
    m_houtmax[2] /=diam;

//----------------------------------------------------------------------
// two hierarchy grids
//----------------------------------------------------------------------
	m_cellSizeh[0]=1.55;
    m_cellSizeh[1]=1.55;
    m_cellSizeh[2]=1.55;

    m_cellSizel[0]=1.55/3.0;
    m_cellSizel[1]=1.55/3.0;
    m_cellSizel[2]=1.55/3.0;

    m_gridSizeh[0] = m_gridSize[0];
    m_gridSizeh[1] = m_gridSize[1];
    m_gridSizeh[2] = m_gridSize[2];

    m_gridSizel[0] = m_gridSizeh[0]*3;
    m_gridSizel[1] = m_gridSizeh[1]*3;
    m_gridSizel[2] = m_gridSizeh[2]*3;

    m_gridSizel_global[0] = m_gridSize_global[0]*3;
    m_gridSizel_global[1] = m_gridSize_global[1]*3;
    m_gridSizel_global[2] = m_gridSize_global[2]*3;

	 if(DivX==1)
    {
    m_gridStartl=m_gridStart*3*3*3;
    m_Ngridlh= m_gridStartl+3*m_gridSizel[1]*m_gridSizel[2]-1; //global value,point to the start of Lh (before Lh)
    m_Ngridlb= m_gridStartl+6*m_gridSizel[1]*m_gridSizel[2]-1;
    m_Ngridrbs=m_gridStartl+(m_gridSizel[0]-6)*m_gridSizel[1]*m_gridSizel[2]; //global value,point to the start of rb (the beginning of rb)
    m_Ngridrb= m_gridStartl+(m_gridSizel[0]-3)*m_gridSizel[1]*m_gridSizel[2]-1;
    m_Ngridrh= m_gridStartl+m_gridSizel[0]*m_gridSizel[1]*m_gridSizel[2]-1;
    }
	 else if(DivZ==1)
    {
    m_gridStartl=m_gridStart*3*3*3;
    m_Ngridlh= m_gridStartl+3*m_gridSizel[0]*m_gridSizel[1]-1; //global value,point to the start of Lh (before Lh)
    m_Ngridlb= m_gridStartl+6*m_gridSizel[0]*m_gridSizel[1]-1;
    m_Ngridrbs=m_gridStartl+(m_gridSizel[2]-6)*m_gridSizel[0]*m_gridSizel[1]; //global value,point to the start of rb (the beginning of rb)
    m_Ngridrb= m_gridStartl+(m_gridSizel[2]-3)*m_gridSizel[0]*m_gridSizel[1]-1;
    m_Ngridrh= m_gridStartl+m_gridSizel[2]*m_gridSizel[0]*m_gridSizel[1]-1;
    }

    printf("low level rank=%8d,m_gridStartl=%8d,m_Ngridlh=%8d,m_Ngridlb=%8d,Ngridrbs=%8d,Ngridrb=%8d,Ngridrh=%8d,m_nGridCellsl=%8d\n",
           rank,m_gridStartl,m_Ngridlh,m_Ngridlb,m_Ngridrbs,m_Ngridrb,m_Ngridrh,m_nGridCellsl);

    printf("rank=%8d,worldOrigin=%6.3f,%6.3f,%6.3f\n",rank,m_worldOrigin[0],m_worldOrigin[1],m_worldOrigin[2]);
    printf("rank=%8d,global worldOrigin=%6.3f,%6.3f,%6.3f\n",rank,m_worldOrigin_global[0],m_worldOrigin_global[1],m_worldOrigin_global[2]);

    printf("grid at high grid level=%8d, low level=%8d,total grids=%8d\n",m_nGridCellsh,m_nGridCellsl,m_nGridCells);
    printf("high level gridSizeh[0]=%8d,gridSizeh[1]=%8d,gridSizeh[2]=%8d\n",m_gridSizeh[0],m_gridSizeh[1],m_gridSizeh[2]);
    printf("low level gridSizel[0]=%8d,gridSizel[1]=%8d,gridSizel[2]=%8d\n",m_gridSizel[0],m_gridSizel[1],m_gridSizel[2]);

	if(DivX==1)
       ntot= (int) ((m_gridSize[0]+4.0)/(m_gridSize[0]-2.0)*1.2*ntot_total/numprocs); // need to modify when using loading balance
    else if(DivY==1)
       ntot= (int) ((m_gridSize[1]+4.0)/(m_gridSize[1]-2.0)*1.2*ntot_total/numprocs);
    else if(DivZ==1)
       ntot= (int) ((m_gridSize[2]+4.0)/(m_gridSize[2]-2.0)*1.2*ntot_total/numprocs);
}


//--- treat mesh data, once for all-----------------------
// calculate wether share plane for each node or edge
// calculate smallest face index for each sharing node or edge
// cellstartWall[grid],WallHash[grid][faceid]

void TreatMesh()
{
//--------------------------------------------------------
// treat meshes on GPU devices
// mesh data -------------------------------------
 printf("treat mesh:");

 Nnodestart=(uint *)malloc(Nmesh*sizeof(uint));
 Nfacestart=(uint *)malloc(Nmesh*sizeof(uint));
 FaceHstart=(uint *)malloc(Nmesh*sizeof(uint));
 Meshstart=(uint *)malloc(Nmesh*sizeof(uint));


  for(imesh=0;imesh<Nmesh;imesh++)
  {
  Nnodestart[imesh] =0;
  Nfacestart[imesh] =0;
  FaceHstart[imesh] =0;
  Meshstart[imesh]  =0;
  }

  for(imesh=1;imesh<Nmesh;imesh++)
  {
  Nnodestart[imesh] = Nnodestart[imesh-1]+ Nnode[imesh-1];
  Nfacestart[imesh] = Nfacestart[imesh-1]+ Nface[imesh-1];

   printf("imesh=%4d,Nnodestart[imesh]=%6d,Nfacestart[imesh]=%6d,Nnode[imesh-1]=%6d,Nface[imesh-1]=%6d\n",
           imesh,Nnodestart[imesh],Nfacestart[imesh],Nnode[imesh-1],Nface[imesh-1]);
  }


  cudaMalloc((void**)&m_dNnodestart,Nmesh*sizeof(uint));
  cudaMalloc((void**)&m_dNfacestart,Nmesh*sizeof(uint));
  cudaMalloc((void**)&m_dFaceHstart,Nmesh*sizeof(uint));
  cudaMalloc((void**)&m_dMeshstart,Nmesh*sizeof(uint));

  //cudaMalloc((void**)&m_Gridinv,Nmesh*sizeof(uint));
   m_Gridinv=(uint *)malloc(Nmesh*sizeof(uint));

 uint Nnodetotal=Nnodestart[Nmesh-1]+Nnode[Nmesh-1];
 uint Nfacetotal=Nfacestart[Nmesh-1]+Nface[Nmesh-1];

//----------------------------------------------------------
 /*
   for(uint imesh=0;imesh<Nmesh;imesh++)
   {
   m_meshSize[imesh][0]=(m_meshEnd[imesh][0]-m_meshOrigin[imesh][0])/m_cellSize[0]+1;
   m_meshSize[imesh][1]=(m_meshEnd[imesh][1]-m_meshOrigin[imesh][1])/m_cellSize[1]+1;
   m_meshSize[imesh][2]=(m_meshEnd[imesh][2]-m_meshOrigin[imesh][2])/m_cellSize[2]+1;
 
   m_nmeshCells[imesh]=m_meshSize[imesh][0]*m_meshSize[imesh][1]*m_meshSize[imesh][2];
   if(imesh>0) Meshstart[imesh]  = Meshstart[imesh-1] + m_nmeshCells[imesh-1];

   printf("mesh %4d,m_meshOrigin=%6.3f,%6.3f,%6.3f,m_meshEnd=%6.3f,%6.3f,%6.3f\n",
   imesh,m_meshOrigin[imesh][0],m_meshOrigin[imesh][1],m_meshOrigin[imesh][2],m_meshEnd[imesh][0],m_meshEnd[imesh][1],m_meshEnd[imesh][2]);

   printf("mesh %4d, meshSize.x=%4d,%4d,%4d,m_nmeshGridCells=%8d,Meshstart[imesh]=%8d\n",
   imesh,m_meshSize[imesh][0],m_meshSize[imesh][1],m_meshSize[imesh][2],m_nmeshCells[imesh],Meshstart[imesh]); 
   }*/
 //----------------------------------------------------------
 // FOR Multiple processes
 //----------------------------------------------------------
   for(imesh=0;imesh<Nmesh;imesh++)
   {
     if(DivX==1)
     {
      if(m_meshOrigin[imesh][0]>Dxmax || m_meshEnd[imesh][0]<Dxmin) //out of domain
      {
      m_meshOrigin[imesh][0]=Dxmin;
      m_meshEnd[imesh][0]=Dxmin;
      }
      if(m_meshOrigin[imesh][0]<Dxmin)m_meshOrigin[imesh][0]=Dxmin;
      if(m_meshEnd[imesh][0]>Dxmax)m_meshEnd[imesh][0]=Dxmax;
      if(i==0)m_meshOrigin[imesh][2]=m_hGlim.poszmin;

     /* m_meshOrigin[imesh][0]=m_worldOrigin[0];*/
      m_meshSize[imesh][0]=(m_meshEnd[imesh][0]-m_meshOrigin[imesh][0])/m_cellSize[0]+1;
      m_meshSize[imesh][1]=(m_meshEnd[imesh][1]-m_meshOrigin[imesh][1])/m_cellSize[1]+1;
      m_meshSize[imesh][2]=(m_meshEnd[imesh][2]-m_meshOrigin[imesh][2])/m_cellSize[2]+1;
     }
    else if(DivZ==1)
     {
      if(m_meshOrigin[imesh][2]>Dzmax || m_meshEnd[imesh][2]<Dzmin) //out of domain
      {
      m_meshOrigin[imesh][2]=Dzmin;
      m_meshEnd[imesh][2]=Dzmin;
      }
      if(m_meshOrigin[imesh][2]<Dzmin)m_meshOrigin[imesh][2]=Dzmin;
      if(m_meshEnd[imesh][2]>Dzmax)m_meshEnd[imesh][2]=Dzmax;
      if(i==0)m_meshOrigin[imesh][2]=m_hGlim.poszmin;

    // m_meshOrigin[imesh][2]=m_worldOrigin[2];
     m_meshSize[imesh][0]=(m_meshEnd[imesh][0]-m_meshOrigin[imesh][0])/m_cellSize[0]+1;
     m_meshSize[imesh][1]=(m_meshEnd[imesh][1]-m_meshOrigin[imesh][1])/m_cellSize[1]+1;
     m_meshSize[imesh][2]=(m_meshEnd[imesh][2]-m_meshOrigin[imesh][2])/m_cellSize[2]+1;
     }

   printf("rank=%8d,mesh %4d,m_meshOrigin=%6.3f,%6.3f,%6.3f,m_meshEnd=%6.3f,%6.3f,%6.3f\n",
   rank,imesh,m_meshOrigin[imesh][0],m_meshOrigin[imesh][1],m_meshOrigin[imesh][2],m_meshEnd[imesh][0],m_meshEnd[imesh][1],m_meshEnd[imesh][2]);

   printf("rank=%8d,mesh %4d, meshSize.x=%4d,%4d,%4d,m_nmeshGridCells=%8d,Meshstart[imesh]=%8d\n",
   rank,imesh,m_meshSize[imesh][0],m_meshSize[imesh][1],m_meshSize[imesh][2],m_nmeshCells[imesh],Meshstart[imesh]);

   m_nmeshCells[imesh]=m_meshSize[imesh][0]*m_meshSize[imesh][1]*m_meshSize[imesh][2];
   if(imesh>0) Meshstart[imesh]  = Meshstart[imesh-1] + m_nmeshCells[imesh-1];

   printf("rank=%8d,mesh %4d,m_meshOrigin=%6.3f,%6.3f,%6.3f,m_meshEnd=%6.3f,%6.3f,%6.3f\n",
   rank,imesh,m_meshOrigin[imesh][0],m_meshOrigin[imesh][1],m_meshOrigin[imesh][2],m_meshEnd[imesh][0],m_meshEnd[imesh][1],m_meshEnd[imesh][2]);
   printf("rank=%8d,mesh %4d, meshSize.x=%4d,%4d,%4d,m_nmeshGridCells=%8d,Meshstart[imesh]=%8d\n",
   rank,imesh,m_meshSize[imesh][0],m_meshSize[imesh][1],m_meshSize[imesh][2],m_nmeshCells[imesh],Meshstart[imesh]); 
   }
 
   cudaMalloc((void**)&m_dmeshSize,Nmesh*3*sizeof(uint));
   cudaMalloc((void**)&m_dmeshOrigin,Nmesh*3*sizeof(double));
   cudaMalloc((void**)&m_dmeshEnd,Nmesh*3*sizeof(double));

   cudaMalloc((void**)&m_dPosnode, Nnodetotal*3*sizeof(double)); 
   cudaMalloc((void**)&m_dFnodeid, Nfacetotal*4*sizeof(uint)); 

   cudaMalloc((void**)&m_dFaceHash,(Meshstart[Nmesh-1]+m_nmeshCells[Nmesh-1])*2*sizeof(uint)); 
   cudaMalloc((void**)&m_dFaceHashStart,(Meshstart[Nmesh-1]+m_nmeshCells[Nmesh-1])*sizeof(uint)); 
   cudaMemset(m_dFaceHashStart, 0xffffffff, (Meshstart[Nmesh-1]+m_nmeshCells[Nmesh-1])*sizeof(uint));

   cudaMalloc((void**)&m_dsharePE, Nfacetotal*4*sizeof(uint));
   cudaMalloc((void**)&m_dEdgeHead, Nfacetotal*4*sizeof(uint));
   cudaMalloc((void**)&m_dsharePV, Nnodetotal*sizeof(uint));
   cudaMalloc((void**)&m_dVertHead, Nnodetotal*sizeof(uint));

   cudaMemset(m_dsharePE, 0, Nfacetotal*4*sizeof(uint));
   cudaMemset(m_dEdgeHead, 0xffffffff, Nfacetotal*4*sizeof(uint));
   cudaMemset(m_dsharePV, 0, Nnodetotal*sizeof(uint));
   cudaMemset(m_dVertHead, 0xffffffff, Nnodetotal*sizeof(uint));

  for(imesh=0;imesh<Nmesh;imesh++)
  {
  cpuErrchk(cudaMemcpy(m_dmeshSize+imesh*3,m_meshSize[imesh],3*sizeof(uint),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dmeshOrigin+imesh*3,m_meshOrigin[imesh],3*sizeof(double),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dmeshEnd+imesh*3,m_meshEnd[imesh],3*sizeof(double),cudaMemcpyHostToDevice));

 cpuErrchk(cudaMemcpy(m_dPosnode+Nnodestart[imesh]*3, m_hPosnode[imesh],Nnode[imesh]*3*sizeof(double),cudaMemcpyHostToDevice));
 cpuErrchk(cudaMemcpy(m_dFnodeid+Nfacestart[imesh]*4, m_hFnodeid[imesh],Nface[imesh]*4*sizeof(uint),cudaMemcpyHostToDevice));

 // if(Nnode[imesh]==0) continue;
  double *dPosnode;
  uint *dFnodeid;
  uint WallHashSize=Nface[imesh]*30;

  cpuErrchk(cudaMalloc((void**)&dPosnode, Nnode[imesh]*3*sizeof(double)));
  cpuErrchk(cudaMalloc((void**)&dFnodeid, Nface[imesh]*4*sizeof(uint)));

  cpuErrchk(cudaMemcpy(dPosnode, m_hPosnode[imesh],Nnode[imesh]*3*sizeof(double),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(dFnodeid, m_hFnodeid[imesh],Nface[imesh]*4*sizeof(uint),cudaMemcpyHostToDevice));

 cpuErrchk(cudaMalloc((void**)&sharePE, Nface[imesh]*4*sizeof(uint)));
 cpuErrchk(cudaMalloc((void**)&EdgeHead, Nface[imesh]*4*sizeof(uint)));
 cpuErrchk(cudaMalloc((void**)&sharePV, Nnode[imesh]*sizeof(uint)));
 cpuErrchk(cudaMalloc((void**)&VertHead, Nnode[imesh]*sizeof(uint)));

 cpuErrchk(cudaMalloc((void**)&m_dcellStartWall, m_nmeshCells[imesh]*sizeof(uint)));
 printf("imesh=%8d,Nnode[imesh]=%8d,m_nmeshCells[imesh]=%8d\n",imesh,Nnode[imesh],m_nmeshCells[imesh]);

 cpuErrchk(cudaMemset(sharePE, 0, Nface[imesh]*4*sizeof(uint)));
 cpuErrchk(cudaMemset(EdgeHead, 0xffffffff, Nface[imesh]*4*sizeof(uint)));

 cpuErrchk(cudaMemset(sharePV, 0, Nnode[imesh]*sizeof(uint)));
 cpuErrchk(cudaMemset(VertHead, 0xffffffff, Nnode[imesh]*sizeof(uint)));

 calEdgeInfo(dFnodeid,Nface[imesh],dPosnode,sharePE,EdgeHead);
 calNodeInfo(dFnodeid,Nface[imesh],Nnode[imesh],dPosnode,sharePV,VertHead);

 cpuErrchk(cudaMalloc((void**)&m_dWallHash, WallHashSize*sizeof(uint)));
 cpuErrchk(cudaMemset(m_dWallHash,0xffffffff,WallHashSize*sizeof(uint)));

 uint m_count=0;
 uint *dcount;
 cpuErrchk(cudaMalloc((void**)&dcount, sizeof(uint)));
 cpuErrchk(cudaMemset(dcount, 0, sizeof(uint)));

 //calculate WallHash [grid, FaceID], for a specific face, which faces it will affect  
 calcWallHash(dFnodeid,Nface[imesh],dPosnode,m_meshSize[imesh],
              m_meshOrigin[imesh],m_cellSize,m_dWallHash,WallHashSize,dcount);
 cpuErrchk(cudaMemcpy(&m_count,dcount,sizeof(uint),cudaMemcpyDeviceToHost));
 cpuErrchk(cudaFree(dcount));

 printf("rank=%8d,imesh=%4d,m_count=%8d,WallHashSize=%8d,Nface=%8d\n",rank,imesh,m_count,WallHashSize,Nface[imesh]);

 if(m_count*2>Nface[imesh]*30)
 {
  WallHashSize=m_count*2;

 printf("m_count*2>WallHashSize,rank=%8d,m_count=%8d,Nface[imesh]=%8d\n",
         rank,m_count,Nface[imesh]);

 cpuErrchk(cudaFree(m_dWallHash));
 cpuErrchk(cudaMalloc((void**)&m_dWallHash, WallHashSize*sizeof(uint)));
 cpuErrchk(cudaMemset(m_dWallHash,0xffffffff,WallHashSize*sizeof(uint)));

 cpuErrchk(cudaMalloc((void**)&dcount, sizeof(uint)));
 cpuErrchk(cudaMemset(dcount, 0, sizeof(uint)));
 calcWallHash(dFnodeid,Nface[imesh],dPosnode,m_meshSize[imesh],m_meshOrigin[imesh],m_cellSize,m_dWallHash,WallHashSize,dcount);
 cpuErrchk(cudaMemcpy(&m_count,dcount,sizeof(uint),cudaMemcpyDeviceToHost));
 cpuErrchk(cudaFree(dcount));

 printf("reallocated m_count=%8d,Nface=%8d\n",m_count,Nface[imesh]);
 }

 cpuErrchk(cudaMalloc((void**)&m_dWallHash0, m_count*2*sizeof(uint)));
 cpuErrchk(cudaMalloc((void**)&m_dWallHash1, m_count*2*sizeof(uint)));

 cpuErrchk(cudaMemcpy(m_dWallHash0,m_dWallHash,m_count*2*sizeof(uint),cudaMemcpyDeviceToDevice));
 cpuErrchk(cudaFree(m_dWallHash));

 //sort particles based on hash, after sorting, obtained a sorted hash for [grid, FaceID]
 RadixSort((KeyValuePair *) m_dWallHash0,(KeyValuePair *) m_dWallHash1, m_count, 32);

 cpuErrchk(cudaMemset(m_dcellStartWall,0xffffffff,m_nmeshCells[imesh]*sizeof(uint)));
 findCellStartWall(m_dWallHash0,m_dcellStartWall,m_count);
 cpuErrchk(cudaFree(m_dWallHash1));

 m_Gridinv[imesh]=0;
 uint FashHashSize,Gridv=0;

 uint *dGridv;
 cpuErrchk(cudaMalloc((void**)&dGridv, sizeof(uint)));
 cpuErrchk(cudaMemset(dGridv, 0, sizeof(uint)));

 FashHashSize=m_nmeshCells[imesh]*2;

 cpuErrchk(cudaMalloc((void**)&m_dFaceHash0,FashHashSize*sizeof(uint)));
 cpuErrchk(cudaMemset(m_dFaceHash0, 0xffffffff,FashHashSize*sizeof(uint)));

 //calculate when a particle is locaed in Grid i, which neighbouring faces will be affected (similar to neighbour list creating)
 AffectGrid(m_meshSize[imesh],m_nmeshCells[imesh],maxFacePerGrid,m_dWallHash0,m_dcellStartWall,m_dFaceHash0,FashHashSize,dGridv);


 cpuErrchk(cudaMemcpy(&Gridv,dGridv,sizeof(uint),cudaMemcpyDeviceToHost));
 cpuErrchk(cudaFree(dGridv));
 printf("rank=%4d,imesh=%8d,Gridv=%8d\n",rank,imesh,Gridv);
 
 m_Gridinv[imesh]=Gridv;

 if(FashHashSize<m_Gridinv[imesh]*2)
 {
 
 printf("rank=%4d,m_nmeshCells[imesh]<m_Gridinv*2, imesh=%4d,FashHashSize=%8d,m_nmeshCells[imesh]=%4d,m_Gridinv=%4d\n",
 rank,imesh,FashHashSize,m_nmeshCells[imesh],m_Gridinv[imesh]);

 Gridv=0;
 cpuErrchk(cudaMalloc((void**)&dGridv, sizeof(uint)));
 cpuErrchk(cudaMemset(dGridv, 0, sizeof(uint)));

 FashHashSize=m_Gridinv[imesh]*2;

 cpuErrchk(cudaFree(m_dFaceHash0));
 cpuErrchk(cudaMalloc((void**)&m_dFaceHash0,FashHashSize*sizeof(uint)));
 cpuErrchk(cudaMemset(m_dFaceHash0, 0xffffffff,FashHashSize*sizeof(uint)));

 AffectGrid(m_meshSize[imesh],m_nmeshCells[imesh],maxFacePerGrid,m_dWallHash0,m_dcellStartWall,m_dFaceHash0,FashHashSize,dGridv);
 cpuErrchk(cudaMemcpy(&Gridv,dGridv,sizeof(uint),cudaMemcpyDeviceToHost));
 cpuErrchk(cudaFree(dGridv));
 m_Gridinv[imesh]=Gridv;

 printf("rank=%4d,reallocate FaceHash0, imesh=%4d,FashHashSize=%8d,m_nmeshCells[imesh]=%4d,m_Gridinv=%4d\n",
 rank,imesh,FashHashSize,m_nmeshCells[imesh],m_Gridinv[imesh]);
 }

//exit(1);
  cpuErrchk(cudaFree(m_dWallHash0));
  cpuErrchk(cudaFree(m_dcellStartWall));

 if(imesh>0)
 {
 FaceHstart[imesh] = FaceHstart[imesh-1] +m_Gridinv[imesh-1];
 }

 cpuErrchk(cudaMalloc((void**)&FaceHash,m_Gridinv[imesh]*2*sizeof(uint)));
 cpuErrchk(cudaMemcpy(FaceHash,m_dFaceHash0,m_Gridinv[imesh]*2*sizeof(uint),cudaMemcpyDeviceToDevice));
 cpuErrchk(cudaFree(m_dFaceHash0));

 cpuErrchk(cudaMalloc((void**)&m_dFaceHash01, m_Gridinv[imesh]*2*sizeof(uint)));

 //sort particles based on hash
 RadixSort((KeyValuePair *) FaceHash,(KeyValuePair *) m_dFaceHash01, m_Gridinv[imesh], 32);
 cpuErrchk(cudaFree(m_dFaceHash01));

 cpuErrchk(cudaMalloc((void**)&FaceHashStart, m_nmeshCells[imesh]*sizeof(uint)));
 cpuErrchk(cudaMemset(FaceHashStart,0xffffffff,m_nmeshCells[imesh]*sizeof(uint)));

 findCellStartWall(FaceHash,FaceHashStart,m_Gridinv[imesh]);

 if((FaceHstart[imesh] +m_Gridinv[imesh])*2 > (Meshstart[Nmesh-1]+m_nmeshCells[Nmesh-1])*2)
 {

 printf("not enough memory allocated for FaceHash,Gridinv=%8d,imesh=%8d,FaceHstart[imesh]=%8d %8d %8d\n",
 m_Gridinv[imesh],imesh,FaceHstart[imesh],Meshstart[Nmesh-1],m_nmeshCells[Nmesh-1]);

 uint *dTemp;
 cpuErrchk(cudaMalloc((void**)&dTemp,FaceHstart[imesh]*2*sizeof(uint)));
 cpuErrchk(cudaMemcpy(dTemp,m_dFaceHash,FaceHstart[imesh]*2*sizeof(uint),cudaMemcpyDeviceToDevice));

 cpuErrchk(cudaFree(m_dFaceHash));
 cpuErrchk(cudaMalloc((void**)&m_dFaceHash,(FaceHstart[imesh] +m_Gridinv[imesh])*2*sizeof(uint)));
 cpuErrchk(cudaMemcpy(m_dFaceHash,dTemp,FaceHstart[imesh]*2*sizeof(uint),cudaMemcpyDeviceToDevice));
 cpuErrchk(cudaFree(dTemp));
 }

 cpuErrchk(cudaMemcpy(m_dFaceHash+FaceHstart[imesh]*2,FaceHash,m_Gridinv[imesh]*2*sizeof(uint),cudaMemcpyDeviceToDevice));

 cpuErrchk(cudaMemcpy(m_dFaceHashStart+Meshstart[imesh],FaceHashStart,m_nmeshCells[imesh]*sizeof(uint),cudaMemcpyDeviceToDevice));

 cpuErrchk(cudaMemcpy(m_dsharePE+Nfacestart[imesh]*4,sharePE,Nface[imesh]*4*sizeof(uint),cudaMemcpyDeviceToDevice));
 cpuErrchk(cudaMemcpy(m_dEdgeHead+Nfacestart[imesh]*4,EdgeHead,Nface[imesh]*4*sizeof(uint),cudaMemcpyDeviceToDevice));

 cpuErrchk(cudaMemcpy(m_dsharePV+Nnodestart[imesh],sharePV,Nnode[imesh]*sizeof(uint),cudaMemcpyDeviceToDevice));
 cpuErrchk(cudaMemcpy(m_dVertHead+Nnodestart[imesh],VertHead,Nnode[imesh]*sizeof(uint),cudaMemcpyDeviceToDevice));

 cpuErrchk(cudaFree(FaceHash));
 cpuErrchk(cudaFree(FaceHashStart));
 cpuErrchk(cudaFree(dPosnode));
 cpuErrchk(cudaFree(dFnodeid));

 cpuErrchk(cudaFree(sharePE));
 cpuErrchk(cudaFree(EdgeHead));
 cpuErrchk(cudaFree(sharePV));
 cpuErrchk(cudaFree(VertHead));
 } // for(imesh=0;imesh<Nmesh;imesh++)

  cpuErrchk(cudaMemcpy(m_dNnodestart,Nnodestart,Nmesh*sizeof(uint),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dNfacestart,Nfacestart,Nmesh*sizeof(uint),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dFaceHstart,FaceHstart,Nmesh*sizeof(uint),cudaMemcpyHostToDevice));
  cpuErrchk(cudaMemcpy(m_dMeshstart,Meshstart,Nmesh*sizeof(uint),cudaMemcpyHostToDevice));

  if(rank==0)printf("treat mesh successful!\n");
   MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
}

void TreatBoundHalo()
{
	if(m_numParticles>0)
	{
   if(m_hPidlh != 0xffffffff)
      {
       Nlh=max(0,m_hPidlh+1-m_hprehead);
       Nlb=max(0,m_hPidlb-m_hPidlh);
       }     
   else
    {
      Nlh=0;
      if(m_hPidlb != 0xffffffff)
        Nlb=max(0,m_hPidlb+1-m_hprehead);
       else
        Nlb=0;
     }

   if(m_hPidrh != 0xffffffff)
   { 
    if(m_hPidrb != 0xffffffff) //m_hPidrb include prehead
    {
    if(m_hPidrbs != 0xffffffff)
    Nrb=max(0,m_hPidrb-m_hPidrbs);
    else
    Nrb=max(0,m_hPidrb-m_hprehead+1);

    Nrh=max(0,m_hPidrh-m_hPidrb);
    }
    else
    {
    Nrb=0;
    Nrh=max(0,m_hPidrh-m_hprehead+1);
    }
   }
   else
   {
    if(m_hPidrb != 0xffffffff) //m_hPidrb include prehead
    {
    if(m_hPidrbs != 0xffffffff)
    Nrb=max(0,m_hPidrb-m_hPidrbs);
    else
    Nrb=max(0,m_hPidrb-m_hprehead+1);
    Nrh=max(0,m_numParIn+m_hprehead-1-m_hPidrb);
    }
    else
    {
    if(m_hPidrbs != 0xffffffff)
    Nrb=max(0,m_numParIn+m_hprehead-1-m_hPidrbs);
    else
    Nrb=0;
    Nrh=0;
    }
   }
  
   if(Nrb>m_numParIn){
   printf("rank=%8d,Nlh=%8d,Nlb=%8d,Nrb=%8d,Nrh=%8d,m_numParIn=%8d\n",rank,Nlh,Nlb,Nrb,Nrh,m_numParIn);
   printf("rank=%8d,Nlh1=%8d,Nlb1=%8d,Nrb1=%8d,Nrh1=%8d\n",rank,Nlh1,Nlb1,Nrb1,Nrh1);
   printf("rank=%4d,m_nGridCellsh=%8d,m_nGridCellsl=%8d,m_hPidlh=%8d,m_hPidlb=%8d,m_hPidrbs=%8d,m_hPidrb=%8d,m_hPidrh=%8d\n",
        rank,m_nGridCellsh,m_nGridCellsl,m_hPidlh,m_hPidlb,m_hPidrbs,m_hPidrb,m_hPidrh);
   }
   } //m_numParticles>0
	else
	{
	Nlh=0;
    Nlb=0;
	Nrb=0;
    Nrh=0;	
	}


   MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
// define the particle number to be receive from neighbour process and save in left and right bounary and halo 
   MPI_Sendrecv(&Nrb,1,MPI_INT,right,1,&Nlh1,1,MPI_INT,left,1,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(&Nrh,1,MPI_INT,right,2,&Nlb1,1,MPI_INT,left,2,MPI_COMM_WORLD,&status);

// date move left
   MPI_Sendrecv(&Nlb,1,MPI_INT,left,3,&Nrh1,1,MPI_INT,right,3,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(&Nlh,1,MPI_INT,left,4,&Nrb1,1,MPI_INT,right,4,MPI_COMM_WORLD,&status);

// renumbering particle id,reserve array size Nlh1+Nlb1 for particles from left side
// note: m_dParticleHash[0] index is saved from prehead to prehead+m_numParticles-1

   if(m_numParIn !=m_numParticles && m_numParticles>0) //flow out fo geometry limit
    {
     idelete=1; 
    // printf("m_numParIn !=m_numParticles,it=%8d,rank=%8d,m_numParticles=%8d,m_numParIn=%8d\n",it,rank,m_numParticles,m_numParIn);
    } 

   m_hPidlh=  Nlh+Nlh1;
   m_hPidlb=  m_hPidlh+Nlb1+Nlb;

   m_hPidrbs=max(0,Nlh1+Nlb1+m_numParIn-Nrb-Nrh);
   if(m_hPidrbs <0)m_hPidrbs=m_hprehead; 

   m_hPidrb=  m_hPidrbs+Nrb+Nrb1;
   m_hPidrh=  m_hPidrb+Nrh+Nrh1; //note that after treatment, Pidrbs, rb, rh point to the start of rb, rh, end of array

   if((m_hPidrbs ==0xffffffff || m_hPidlh==0xffffffff || m_hPidlb==0xffffffff || m_hPidrb==0xffffffff || m_hPidrh==0xffffffff)&& m_numParticles>0)
   printf("after rank=%4d,m_nGridCells=%8d,m_hPidlh=%8d,m_hPidlb=%8d,m_hPidrb=%8d,m_hPidrh=%8d,m_numParticles=%8d\n",
          rank,m_nGridCells,m_hPidlh,m_hPidlb,m_hPidrb,m_hPidrh, m_numParticles);

   if(m_hPidrh< m_numParIn+Nlh1+Nlb1+Nrb1 && m_numParIn>0) 
   printf("m_hPidrh<ParIn,rank=%4d,Pidlh=%8d,Pidlb=%8d,Pidrb=%8d,Pidrh=%8d,Nlh1=%8d,Nlb1=%8d,Nrb1=%8d,Nrh1=%8d,m_numParIn=%8d,allocatenum=%8d\n",
          rank,m_hPidlh,m_hPidlb,m_hPidrb,m_hPidrh,Nlh1,Nlb1,Nrb1,Nrh1,m_numParIn,allocatenum);

   //printf("after rank=%4d,m_hPidlh=%8d,m_hPidlb=%8d,m_hPidrb=%8d,m_hPidrh=%8d,m_numParticles=%8d,Nlh1=%8d,Nlb1=%8d,Nrb1=%8d,Nrh1=%8d\n",
   //       rank,m_hPidlh,m_hPidlb,m_hPidrb,m_hPidrh,m_numParticles,Nlh1,Nlb1,Nrb1,Nrh1);
}