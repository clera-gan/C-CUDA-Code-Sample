#include "dempacking.h"
#include "dempacking.cuh"
//#include "GPUSet.h"

void InitialCUDA()
{
    // Choose which GPU to run on, change this on a multi-GPU system.
    //int   ndev=0;
    //cudaGetDeviceCount(&ndev);
    //cudaSetDevice(devnum);
    printf("Now Initialize GPU successful, rank=%8d\n",
		 rank);
	SetCUDADevice(lr%num_devices);
    printf("Initialize GPU successful, rank=%8d, set device index to %4d,lr=%4d,num_devices=%4d\n",
		 rank,lr%num_devices,lr,num_devices);

    //cudaRuntimeGetVersion(&ndev);
    //cudaDriverGetVersion(&ndev);

}

void AllocateGPUArray()
{
  if(m_hUSE_HIERARCHY==1)
   {cpuErrchk(cudaMalloc( (void**)&m_dCellStart, m_nGridCells*sizeof(uint)));

   if(DivX==1)
    allocateArray( (void **)&m_dCellStartB,(m_gridSizel[1]*m_gridSizel[2])*12*sizeof(uint)); //devide in x direction
    else if(DivZ==1)
    allocateArray( (void **)&m_dCellStartB,(m_gridSizel[0]*m_gridSizel[1])*12*sizeof(uint)); //devide in z direction  
   
  }
  else
   {
	cpuErrchk(cudaMalloc( (void**)&m_dCellStart, m_nGridCellsh*sizeof(uint)));

    if(DivX==1)
    allocateArray( (void **)&m_dCellStartB,(m_gridSizeh[1]*m_gridSizeh[2])*12*sizeof(uint)); //devide in x direction
    else if(DivZ==1)
    allocateArray( (void **)&m_dCellStartB,(m_gridSizeh[0]*m_gridSizeh[1])*12*sizeof(uint)); //devide in z direction  
  }

    cpuErrchk(cudaMalloc( (void**)&Addcount, sizeof(uint)));

    cudaMalloc( (void**)&m_dmaxCn, nsize*sizeof(uint));
    cudaMemcpy(m_dmaxCn, m_hmaxCn,nsize*sizeof(uint),cudaMemcpyHostToDevice);

    cudaMalloc( (void**)&m_dsize, nsize*sizeof(double));
    cudaMemcpy(m_dsize,Size,nsize*sizeof(double),cudaMemcpyHostToDevice);


}



void CopyConstantsToCUDA()
{
	/*
size_t size;
size=sizeof(GeoDivide)+sizeof(Wallproperties)+sizeof(Geolimit)+sizeof(Radlimit)+\
	NumOfMat*sizeof(Matproperties)+NumOfFeed*sizeof(Feedlimit)+\
	NumOfTranslation*sizeof(Translation)+NumOfRotation*sizeof(Rotation)+\
	NumOfVibration*sizeof(Vibration)+sizeof(uint);

printf("size of constant memory is %8d\n",size);*/
	//printf("feed limit GPU %4d: %lf %lf %lf %lf %lf %lf \n",
   // i,m_hfeedlim[0].xmin,m_hfeedlim[0].xmax,m_hfeedlim[0].ymin,m_hfeedlim[0].ymax,m_hfeedlim[0].zmin,m_hfeedlim[0].zmax);

	 copyHostToConstmem(&m_hGDiv,sizeof(GeoDivide),
		                &m_hParW,sizeof(Wallproperties),
                        &m_hGlim,sizeof(Geolimit),
                        &m_hRlim,sizeof(Radlimit),

                        m_hmatpro,NumOfMat*sizeof(Matproperties),
                        m_hfeedlim,NumOfFeed*sizeof(Feedlimit),

						m_hTranslation,NumOfTranslation*sizeof(Translation),
                        //m_hRotation,NumOfRotation*sizeof(Rotation),
                        //m_hVibration,NumOfVibration*sizeof(Vibration),

                        &m_hUSE_HIERARCHY, sizeof(uint));
}

