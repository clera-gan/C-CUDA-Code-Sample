#include "dempacking.h"
//#include "Materials.h"

void SetMatMaxRad()
{
	printf("SetMatMaxRad. rank=%4d\n",rank);
    MatMaxRad=(double *)malloc(NumOfMat*sizeof(double));

    for(i=0;i<NumOfMat;i++)
    {
    double radmax=0.0;
    for(j=0;j<m_hmatSize[i].sizeNum;j++)
    {
     radmax= fmax(radmax,m_hmatSize[i].sizedis[j].dia);
    }
     MatMaxRad[i]=0.5*radmax/diam;
     printf("The maximum rad of mat %4d is %15.6f\n",i,MatMaxRad[i]);
     }
}

void AllocateMaterialArray()
{
    m_hMatrad=(double **)malloc(NumOfMat*sizeof(double *));
    m_hMatrmass=(double **)malloc(NumOfMat*sizeof(double *));
    m_hMatinert=(double **)malloc(NumOfMat*sizeof(double *));
}

void InitialMaterialArray()
{
    m_hMatrad=(double **)malloc(NumOfMat*sizeof(double *));
    m_hMatrmass=(double **)malloc(NumOfMat*sizeof(double *));
    m_hMatinert=(double **)malloc(NumOfMat*sizeof(double *));


    for(i=0;i<NumOfMat;i++)
    {
    m_hMatrad[i]=(double *)malloc(m_hmatSize[i].TotalAddNum*sizeof(double));
    m_hMatrmass[i]=(double *)malloc(m_hmatSize[i].TotalAddNum*sizeof(double));
    m_hMatinert[i]=(double *)malloc(m_hmatSize[i].TotalAddNum*sizeof(double));

    ip=0;
    for(j=0;j<m_hmatSize[i].sizeNum;j++)
    {
    for(k=0;k<m_hmatSize[i].sizedis[j].Num;k++)
      m_hMatrad[i][ip++]=m_hmatSize[i].sizedis[j].dia*0.5/diam;
     
    }
    std::random_shuffle(m_hMatrad[i],m_hMatrad[i]+m_hmatSize[i].TotalAddNum);

    for(j=0;j<m_hmatSize[i].TotalAddNum;j++)
    {
     m_hMatrmass[i][j]=8*m_hMatrad[i][j]*m_hMatrad[i][j]*m_hMatrad[i][j]*m_hmatpro[i].denp/denp; //check here denp
     m_hMatinert[i][j]=0.4*m_hMatrmass[i][j]*(m_hMatrad[i][j]*m_hMatrad[i][j])*1.0;
     /* printf("i=%4d,j=%4d,rad=%6.3f,rmass=%6.3f,inert=%6.3f\n",
              i,j,m_hMatrad[i][j],m_hMatrmass[i][j],m_hMatinert[i][j]);*/
    }
    }
}
