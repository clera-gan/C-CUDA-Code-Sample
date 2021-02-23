#include "dempacking.h"
//#include "mpiFunctions.h"

void StartupMPI()  
{
   char processor_name[MPI_MAX_PROCESS_NAME];

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   rank=0;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank); // global rank, eg. 0,1,2,3,4,5,6,7,8,....
   MPI_Get_processor_name(processor_name,&namelen);
  // fprintf(stdout,"Process %d of %d is on %s \n",rank,numprocs,processor_name);
//   fflush(stdout); 
 //  exit(1);

// GPU start 
   num_devices=CountCUDADevice();  // num_devices= ranks per node=4 
   lr = rank %num_devices;   // local rank, eg. 0,1,2,3
   printf("rank=%4d,num_devices=%4d,lr=%4d,devnum=%4d\n",rank,num_devices,lr,lr%num_devices);
}

void mpiReduceParNum()
 {
	 MPI_Allreduce(&m_numParticles,&m_numParticles_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
     if(m_numParticles_total>ntot_total && rank==0)m_numParticles_total=ntot_total; //modify here

     //if any of the process need to update neighbour list, update all the process
     MPI_Allreduce(&m_hiall,&m_hiall_total,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);  
 }

void mpiDataSendRecv()
{
MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
	//send and receive date
if(Nrb !=0)
{
cpuErrchk(cudaMemcpy(m_hParticleHashsbuf+(Nlh+Nlb)*2,m_dParticleHash[0]+(m_hPidrbs)*2,(Nrb)*2*sizeof(uint),cudaMemcpyDeviceToHost)); //right data	
MPI_CHECK(MPI_Send(m_hParticleHashsbuf+(Nlh+Nlb)*2,Nrb*2,MPI_UNSIGNED,right,10,MPI_COMM_WORLD)); 
}
if(Nrh !=0 )
{
cpuErrchk(cudaMemcpy(m_hParticleHashsbuf+(Nlh+Nlb+Nrb)*2,m_dParticleHash[0]+(m_hPidrb)*2,(Nrh)*2*sizeof(uint),cudaMemcpyDeviceToHost));  //right data
MPI_Send(m_hParticleHashsbuf+(Nlh+Nlb+Nrb)*2,Nrh*2,MPI_UNSIGNED,right,11,MPI_COMM_WORLD); 
}
if(Nlh1 !=0){
MPI_Recv(m_hParticleHashrbuf,Nlh1*2,MPI_UNSIGNED,left,10,MPI_COMM_WORLD,&status);
cpuErrchk(cudaMemcpy(m_dParticleHash[0],m_hParticleHashrbuf,Nlh1*2*sizeof(uint),cudaMemcpyHostToDevice));
}
if(Nlb1 !=0 ){
MPI_Recv(m_hParticleHashrbuf+Nlh1*2,Nlb1*2,MPI_UNSIGNED,left,11,MPI_COMM_WORLD,&status);
cpuErrchk(cudaMemcpy(m_dParticleHash[0]+(Nlh1+Nlh)*2,m_hParticleHashrbuf+Nlh1*2,Nlb1*2*sizeof(uint),cudaMemcpyHostToDevice));
}

// m_dParticleHash[0] data move to left: 
if(Nlh !=0)
{
cpuErrchk(cudaMemcpy(m_hParticleHashsbuf,m_dParticleHash[0]+Nlh1*2,(Nlh)*2*sizeof(uint),cudaMemcpyDeviceToHost)); //left data
MPI_Send(m_hParticleHashsbuf,Nlh*2,MPI_UNSIGNED,left,12,MPI_COMM_WORLD);
}

if(Nlb !=0)
{
cpuErrchk(cudaMemcpy(m_hParticleHashsbuf+(Nlh)*2,m_dParticleHash[0]+(Nlh1+Nlh+Nlb1)*2,(Nlb)*2*sizeof(uint),cudaMemcpyDeviceToHost)); //left data
MPI_Send(m_hParticleHashsbuf+Nlh*2,Nlb*2,MPI_UNSIGNED,left,13,MPI_COMM_WORLD);
}

if(Nrb1 !=0)
{
MPI_Recv(m_hParticleHashrbuf+(Nlh1+Nlb1)*2,Nrb1*2,MPI_UNSIGNED,right,12,MPI_COMM_WORLD,&status);
cpuErrchk(cudaMemcpy(m_dParticleHash[0]+(m_hPidrbs+Nrb)*2,m_hParticleHashrbuf+(Nlh1+Nlb1)*2,Nrb1*2*sizeof(uint),cudaMemcpyHostToDevice));
}

if(Nrh1 !=0){
MPI_Recv(m_hParticleHashrbuf+(Nlh1+Nlb1+Nrb1)*2,Nrh1*2,MPI_UNSIGNED,right,13,MPI_COMM_WORLD,&status);
cpuErrchk(cudaMemcpy(m_dParticleHash[0]+(m_hPidrb+Nrh)*2,m_hParticleHashrbuf+(Nlh1+Nlb1+Nrb1)*2,Nrh1*2*sizeof(uint),cudaMemcpyHostToDevice));
}

   m_numParticles=max(m_numParticles-(Nlh+Nrh)+(Nlb1+Nrb1),0);
  /* printf("max=%8d %8d %8d %8d %8d %8d\n",
	   max(m_numParticles-(Nlh+Nrh)+(Nlb1+Nrb1),0),m_numParticles,Nlh,Nrh,Nlb1,Nrb1);*/
   NParticlesrank=m_numParticles+Nlh1+Nlh+Nrh+Nrh1;
   m_hprehead=m_hPidlh;

   if(NParticlesrank>allocatenum)
   {
   printf("not enough memory for arrays, need to allocate!,NParticlesrank=%8d,allocatenum=%8d\n",NParticlesrank,allocatenum);
   }

   MPI_Allreduce(&m_numParticles,&m_numParticles_total,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);

   //printf("it=%8d,rank=%4d,NParticlesrank=%8d,allocatenum=%8d,m_numParticles=%8d,m_numParticles_total=%8d,m_hprehead=%4d,Nlh=%8d,Nlb=%8d,Nrb=%8d,Nrh=%8d,Nlh1=%8d,Nlb1=%8d,Nrb1=%8d,Nrh1=%8d\n",
	//   it,rank,NParticlesrank,allocatenum,m_numParticles,m_numParticles_total,m_hprehead,Nlh,Nlb,Nrb,Nrh,Nlh1,Nlb1,Nrb1,Nrh1);

// pos data move to right 
if(Nrb!=0)
{
cpuErrchk(cudaMemcpy(m_hpossbuf +(Nlh+Nlb)*3,m_dpos +(m_hPidrbs)*3,(Nrb)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hpdissbuf+(Nlh+Nlb)*3,m_dpdis+(m_hPidrbs)*3,(Nrb)*3*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hangvsbuf+(Nlh+Nlb)*3,m_dangv+(m_hPidrbs)*3,(Nrb)*3*sizeof(double),cudaMemcpyDeviceToHost));  

cpuErrchk(cudaMemcpy(m_hradsbuf  +(Nlh+Nlb),m_drad  +(m_hPidrbs),(Nrb)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hrmasssbuf+(Nlh+Nlb),m_drmass+(m_hPidrbs),(Nrb)*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hinertsbuf+(Nlh+Nlb),m_dinert+(m_hPidrbs),(Nrb)*sizeof(double),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hmatIdsbuf+(Nlh+Nlb),m_dmatId+(m_hPidrbs),(Nrb)*sizeof(uint),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hcontactEndsbuf+(Nlh+Nlb)*3,m_dcontactEnd+(m_hPidrbs)*3,(Nrb)*3*sizeof(uint),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hEngdispsbuf+(Nlh+Nlb)*3,m_dEngdisp+(m_hPidrbs)*3,(Nrb)*3*sizeof(double),cudaMemcpyDeviceToHost)); 

MPI_CHECK(MPI_Send(m_hpossbuf +(Nlh+Nlb)*3,Nrb*3,MPI_DOUBLE,right,14,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hpdissbuf+(Nlh+Nlb)*3,Nrb*3,MPI_DOUBLE,right,15,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hangvsbuf+(Nlh+Nlb)*3,Nrb*3,MPI_DOUBLE,right,16,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hradsbuf  +(Nlh+Nlb),Nrb,MPI_DOUBLE,right,17,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hrmasssbuf+(Nlh+Nlb),Nrb,MPI_DOUBLE,right,18,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hinertsbuf+(Nlh+Nlb),Nrb,MPI_DOUBLE,right,19,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hmatIdsbuf+(Nlh+Nlb),Nrb,MPI_UNSIGNED,right,20,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hcontactEndsbuf+(Nlh+Nlb)*3,Nrb*3,MPI_UNSIGNED,right,120,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispsbuf+(Nlh+Nlb)*3,Nrb,MPI_DOUBLE,right,121,MPI_COMM_WORLD)); 
}

if(Nlh1 != 0)
{
MPI_CHECK(MPI_Recv(m_hposrbuf, Nlh1*3,MPI_DOUBLE,left,14,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hpdisrbuf,Nlh1*3,MPI_DOUBLE,left,15,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hangvrbuf,Nlh1*3,MPI_DOUBLE,left,16,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hradrbuf,  Nlh1,MPI_DOUBLE,left,17,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hrmassrbuf,Nlh1,MPI_DOUBLE,left,18,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hinertrbuf,Nlh1,MPI_DOUBLE,left,19,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hmatIdrbuf,Nlh1,MPI_UNSIGNED,left,20,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hcontactEndrbuf,Nlh1*3,MPI_UNSIGNED,left,120,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdisprbuf,Nlh1*3,MPI_DOUBLE,left,121,MPI_COMM_WORLD,&status));

cpuErrchk(cudaMemcpy(m_dpos, m_hposrbuf, Nlh1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dpdis,m_hpdisrbuf,Nlh1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dangv,m_hangvrbuf,Nlh1*3*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_drad,  m_hradrbuf,  Nlh1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_drmass,m_hrmassrbuf,Nlh1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dinert,m_hinertrbuf,Nlh1*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dmatId,m_hmatIdrbuf,Nlh1*sizeof(uint),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dcontactEnd,m_hcontactEndrbuf,Nlh1*3*sizeof(uint),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdisp,m_hEngdisprbuf,Nlh1*3*sizeof(double),cudaMemcpyHostToDevice));
}

if(Nrh!=0)
{
cpuErrchk(cudaMemcpy(m_hpossbuf +(Nlh+Nlb+Nrb)*3,m_dpos+(m_hPidrb)*3,(Nrh)*3*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hpdissbuf+(Nlh+Nlb+Nrb)*3,m_dpdis+(m_hPidrb)*3,(Nrh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hangvsbuf+(Nlh+Nlb+Nrb)*3,m_dangv+(m_hPidrb)*3,(Nrh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hradsbuf  +(Nlh+Nlb+Nrb),m_drad+(m_hPidrb),(Nrh)*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hrmasssbuf+(Nlh+Nlb+Nrb),m_drmass+(m_hPidrb),(Nrh)*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hinertsbuf+(Nlh+Nlb+Nrb),m_dinert+(m_hPidrb),(Nrh)*sizeof(double),cudaMemcpyDeviceToHost));  

/*if(rank==0)
for(i=0;i<Nrh;i++)
	printf("Nrh m_hradsbuf[%4d]=%10.6f,Nlh+Nlb+Nrb=%8d,m_hPidrb=%8d,Nrh=%8d\n",
	i,m_hradsbuf[i+Nlh+Nlb+Nrb],Nlh+Nlb+Nrb,m_hPidrb,Nrh);*/

cpuErrchk(cudaMemcpy(m_hmatIdsbuf+(Nlh+Nlb+Nrb),m_dmatId+(m_hPidrb),(Nrh)*sizeof(uint),cudaMemcpyDeviceToHost));  

cpuErrchk(cudaMemcpy(m_hcontactEndsbuf+(Nlh+Nlb+Nrb)*3,m_dcontactEnd+(m_hPidrb)*3,(Nrh)*3*sizeof(uint),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hEngdispsbuf+(Nlh+Nlb+Nrb)*3,m_dEngdisp+(m_hPidrb)*3,(Nrh)*3*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hEngdispwsbuf+(Nlh+Nlb+Nrb),m_dEngdispw+(m_hPidrb),(Nrh)*sizeof(double),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hEngdispVarwsbuf+(Nlh+Nlb+Nrb)*4,m_dEngdispVarw+(m_hPidrb)*4,(Nrh)*4*sizeof(double),cudaMemcpyDeviceToHost));  

//copyArrayFromDevice(m_hsjgisbuf+(Nlh),m_dsjgi+(m_hPidrb),(Nrh)*sizeof(uint),cudaMemcpyDeviceToHost));  
cpuErrchk(cudaMemcpy(m_hdisptwsbuf+(Nlh)*3,m_ddisptw+(m_hPidrb)*3,(Nrh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hfricpwsbuf+(Nlh),m_dfricpw+(m_hPidrb),(Nrh)*sizeof(double),cudaMemcpyDeviceToHost));  

MPI_CHECK(MPI_Send(m_hpossbuf +(Nlh+Nlb+Nrb)*3,Nrh*3,MPI_DOUBLE,right,22,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hpdissbuf+(Nlh+Nlb+Nrb)*3,Nrh*3,MPI_DOUBLE,right,23,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hangvsbuf+(Nlh+Nlb+Nrb)*3,Nrh*3,MPI_DOUBLE,right,24,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hradsbuf+(Nlh+Nlb+Nrb),Nrh,MPI_DOUBLE,right,25,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hrmasssbuf+(Nlh+Nlb+Nrb),Nrh,MPI_DOUBLE,right,26,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hinertsbuf+(Nlh+Nlb+Nrb),Nrh,MPI_DOUBLE,right,27,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hmatIdsbuf+(Nlh+Nlb+Nrb),Nrh,MPI_UNSIGNED,right,28,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hcontactEndsbuf+(Nlh+Nlb+Nrb)*3,Nrh*3,MPI_UNSIGNED,right,129,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispsbuf+(Nlh+Nlb+Nrb)*3,Nrh*3,MPI_DOUBLE,right,130,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispwsbuf+(Nlh+Nlb+Nrb),Nrh,MPI_DOUBLE,right,131,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispVarwsbuf+(Nlh+Nlb+Nrb)*4,Nrh*4,MPI_DOUBLE,right,132,MPI_COMM_WORLD)); 


//MPI_Send(m_hsjgisbuf+Nlh,Nrh,MPI_UNSIGNED,right,30,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hdisptwsbuf+Nlh*3,Nrh*3,MPI_DOUBLE,right,31,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hfricpwsbuf+Nlh,Nrh,MPI_DOUBLE,right,32,MPI_COMM_WORLD)); 
}

if(Nlb1 != 0)
{
MPI_CHECK(MPI_Recv(m_hposrbuf+Nlh1*3, Nlb1*3,MPI_DOUBLE,left,22,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hpdisrbuf+Nlh1*3,Nlb1*3,MPI_DOUBLE,left,23,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hangvrbuf+Nlh1*3,Nlb1*3,MPI_DOUBLE,left,24,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hradrbuf+Nlh1,  Nlb1,MPI_DOUBLE,left,25,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hrmassrbuf+Nlh1,Nlb1,MPI_DOUBLE,left,26,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hinertrbuf+Nlh1,Nlb1,MPI_DOUBLE,left,27,MPI_COMM_WORLD,&status));

/*if(rank==1)
for(i=0;i<Nlb1;i++)
	printf("m_hradrbuf[%4d]=%10.6f,Nlh1=%8d,Nlb1=%4d\n",i,m_hradrbuf[i+Nlh1],Nlh1,Nlb1);*/

MPI_CHECK(MPI_Recv(m_hmatIdrbuf+Nlh1,Nlb1,MPI_UNSIGNED,left,28,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hcontactEndrbuf+Nlh1*3,Nlb1*3,MPI_UNSIGNED,left,129,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdisprbuf+Nlh1*3,Nlb1*3,MPI_DOUBLE,left,130,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdispwrbuf+Nlh1,Nlb1,MPI_DOUBLE,left,131,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdispVarwrbuf+Nlh1*4,Nlb1*4,MPI_DOUBLE,left,132,MPI_COMM_WORLD,&status));

//MPI_Recv(m_hsjgirbuf,Nlb1,MPI_UNSIGNED,left,30,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hdisptwrbuf,Nlb1*3,MPI_DOUBLE,left,31,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hfricpwrbuf,Nlb1,MPI_DOUBLE,left,32,MPI_COMM_WORLD,&status));

cpuErrchk(cudaMemcpy(m_dpos +(Nlh1+Nlh)*3,m_hposrbuf +Nlh1*3,Nlb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dpdis+(Nlh1+Nlh)*3,m_hpdisrbuf+Nlh1*3,Nlb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dangv+(Nlh1+Nlh)*3,m_hangvrbuf+Nlh1*3,Nlb1*3*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_drad +(Nlh1+Nlh), m_hradrbuf +Nlh1,Nlb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_drmass+(Nlh1+Nlh),m_hrmassrbuf+Nlh1,Nlb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dinert+(Nlh1+Nlh),m_hinertrbuf+Nlh1,Nlb1*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dmatId+(Nlh1+Nlh),m_hmatIdrbuf+Nlh1,Nlb1*sizeof(uint),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dcontactEnd+(Nlh1+Nlh)*3,m_hcontactEndrbuf+Nlh1*3,Nlb1*3*sizeof(uint),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdisp+(Nlh1+Nlh)*3,m_hEngdisprbuf+Nlh1*3,Nlb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdispw+(Nlh1+Nlh),m_hEngdispwrbuf+Nlh1,Nlb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdispVarw+(Nlh1+Nlh)*4,m_hEngdispVarwrbuf+Nlh1*4,Nlb1*4*sizeof(double),cudaMemcpyHostToDevice));

//copyArrayToDevice(m_dsjgi+(Nlh1+Nlh),m_hsjgirbuf,Nlb1*sizeof(uint)); //check here
cpuErrchk(cudaMemcpy(m_ddisptw+(Nlh1+Nlh)*3,m_hdisptwrbuf,Nlb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dfricpw+(Nlh1+Nlh),m_hfricpwrbuf,Nlb1*sizeof(double),cudaMemcpyHostToDevice));

}

// pos data move to left: --------------------------------------------------------------------

if(Nlh!=0)
{
cpuErrchk(cudaMemcpy(m_hpossbuf, m_dpos+Nlh1*3, (Nlh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hpdissbuf,m_dpdis+Nlh1*3,(Nlh)*3*sizeof(double),cudaMemcpyDeviceToHost));
cpuErrchk(cudaMemcpy(m_hangvsbuf,m_dangv+Nlh1*3,(Nlh)*3*sizeof(double),cudaMemcpyDeviceToHost));

cpuErrchk(cudaMemcpy(m_hradsbuf,  m_drad+Nlh1,  (Nlh)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hrmasssbuf,m_drmass+Nlh1,(Nlh)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hinertsbuf,m_dinert+Nlh1,(Nlh)*sizeof(double),cudaMemcpyDeviceToHost)); 

/*if(rank==1)
for(i=0;i<Nlh;i++)
	printf("Nlh m_hradsbuf[%4d]=%10.6f,Nlh1=%8d,Nlh=%4d\n",i,m_hradsbuf[i],Nlh1,Nlh);*/

cpuErrchk(cudaMemcpy(m_hmatIdsbuf,m_dmatId+Nlh1,(Nlh)*sizeof(uint),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hcontactEndsbuf,m_dcontactEnd+Nlh1*3,(Nlh)*3*sizeof(uint),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hEngdispsbuf,m_dEngdisp+Nlh1*3,(Nlh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hEngdispwsbuf,m_dEngdispw+Nlh1,(Nlh)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hEngdispVarwsbuf,m_dEngdispVarw+Nlh1*4,(Nlh)*4*sizeof(double),cudaMemcpyDeviceToHost)); 

//copyArrayFromDevice(m_hsjgisbuf,m_dsjgi+Nlh1,(Nlh)*sizeof(uint)); 
cpuErrchk(cudaMemcpy(m_hdisptwsbuf,m_ddisptw+Nlh1*3,(Nlh)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hfricpwsbuf,m_dfricpw+Nlh1,(Nlh)*sizeof(double),cudaMemcpyDeviceToHost)); 

MPI_CHECK(MPI_Send(m_hpossbuf, Nlh*3,MPI_DOUBLE,left,33,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hpdissbuf,Nlh*3,MPI_DOUBLE,left,34,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hangvsbuf,Nlh*3,MPI_DOUBLE,left,35,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hradsbuf,  Nlh,MPI_DOUBLE,left,36,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hrmasssbuf,Nlh,MPI_DOUBLE,left,37,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hinertsbuf,Nlh,MPI_DOUBLE,left,38,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hmatIdsbuf,Nlh,MPI_UNSIGNED,left,39,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hcontactEndsbuf,Nlh*3,MPI_UNSIGNED,left,140,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispsbuf,Nlh*3,MPI_DOUBLE,left,141,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispwsbuf,Nlh,MPI_DOUBLE,left,142,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispVarwsbuf,Nlh*4,MPI_DOUBLE,left,143,MPI_COMM_WORLD)); 

//MPI_Send(m_hsjgisbuf,Nlh,MPI_UNSIGNED,left,41,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hdisptwsbuf,Nlh*3,MPI_DOUBLE,left,42,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hfricpwsbuf,Nlh,MPI_DOUBLE,left,43,MPI_COMM_WORLD)); 
}

if(Nrb1 != 0)
{
MPI_CHECK(MPI_Recv(m_hposrbuf+(Nlh1+Nlb1)*3, Nrb1*3,MPI_DOUBLE,right,33,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hpdisrbuf+(Nlh1+Nlb1)*3,Nrb1*3,MPI_DOUBLE,right,34,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hangvrbuf+(Nlh1+Nlb1)*3,Nrb1*3,MPI_DOUBLE,right,35,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hradrbuf+(Nlh1+Nlb1),  Nrb1,MPI_DOUBLE,right,36,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hrmassrbuf+(Nlh1+Nlb1),Nrb1,MPI_DOUBLE,right,37,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hinertrbuf+(Nlh1+Nlb1),Nrb1,MPI_DOUBLE,right,38,MPI_COMM_WORLD,&status));

/*if(rank==0)
for(i=0;i<Nrb1;i++)
	printf("Nrb1 m_hradrbuf[%4d]=%10.6f,m_hPidrbs=%8d,Nrb=%8d,m_hPidrbs+Nrb=%8d,Nlh1+Nlb1=%8d\n",
	i,m_hradrbuf[Nlh1+Nlb1+i],m_hPidrbs,Nrb,m_hPidrbs+Nrb,Nlh1+Nlb1);*/

MPI_CHECK(MPI_Recv(m_hmatIdrbuf+(Nlh1+Nlb1),Nrb1,MPI_UNSIGNED,right,39,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hcontactEndrbuf+(Nlh1+Nlb1)*3,Nrb1*3,MPI_UNSIGNED,right,140,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdisprbuf+(Nlh1+Nlb1)*3,Nrb1*3,MPI_DOUBLE,right,141,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdispwrbuf+(Nlh1+Nlb1),Nrb1,MPI_DOUBLE,right,142,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdispVarwrbuf+(Nlh1+Nlb1)*4,Nrb1*4,MPI_DOUBLE,right,143,MPI_COMM_WORLD,&status));

//MPI_Recv(m_hsjgirbuf+(Nlb1),Nrb1,MPI_UNSIGNED,right,41,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hdisptwrbuf+(Nlb1)*3,Nrb1*3,MPI_DOUBLE,right,42,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hfricpwrbuf+(Nlb1),Nrb1,MPI_DOUBLE,right,43,MPI_COMM_WORLD,&status));

cpuErrchk(cudaMemcpy(m_dpos +(m_hPidrbs+Nrb)*3,m_hposrbuf+(Nlh1+Nlb1)*3,Nrb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dpdis+(m_hPidrbs+Nrb)*3,m_hpdisrbuf+(Nlh1+Nlb1)*3,Nrb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dangv+(m_hPidrbs+Nrb)*3,m_hangvrbuf+(Nlh1+Nlb1)*3,Nrb1*3*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_drad  +(m_hPidrbs+Nrb),m_hradrbuf+(Nlh1+Nlb1),Nrb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_drmass+(m_hPidrbs+Nrb),m_hrmassrbuf+(Nlh1+Nlb1),Nrb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dinert+(m_hPidrbs+Nrb),m_hinertrbuf+(Nlh1+Nlb1),Nrb1*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dmatId+(m_hPidrbs+Nrb),m_hmatIdrbuf+(Nlh1+Nlb1),Nrb1*sizeof(uint),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dcontactEnd+(m_hPidrbs+Nrb)*3,m_hcontactEndrbuf+(Nlh1+Nlb1)*3,Nrb1*3*sizeof(uint),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdisp+(m_hPidrbs+Nrb)*3,m_hEngdisprbuf+(Nlh1+Nlb1)*3,Nrb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdispw+(m_hPidrbs+Nrb),m_hEngdispwrbuf+(Nlh1+Nlb1),Nrb1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdispVarw+(m_hPidrbs+Nrb)*4,m_hEngdispVarwrbuf+(Nlh1+Nlb1)*4,Nrb1*4*sizeof(double),cudaMemcpyHostToDevice));

//copyArrayToDevice(m_dsjgi+(m_hPidrbs+Nrb),  m_hsjgirbuf+(Nlb1),Nrb1*sizeof(uint),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_ddisptw+(m_hPidrbs+Nrb)*3,m_hdisptwrbuf+(Nlb1)*3,Nrb1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dfricpw+(m_hPidrbs+Nrb),  m_hfricpwrbuf+(Nlb1),Nrb1*sizeof(double),cudaMemcpyHostToDevice));
}

if(Nlb!=0)
{
cpuErrchk(cudaMemcpy(m_hpossbuf+(Nlh)*3, m_dpos +(Nlh1+Nlh+Nlb1)*3,(Nlb)*3*sizeof(double),cudaMemcpyDeviceToHost));
cpuErrchk(cudaMemcpy(m_hpdissbuf+(Nlh)*3,m_dpdis+(Nlh1+Nlh+Nlb1)*3,(Nlb)*3*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hangvsbuf+(Nlh)*3,m_dangv+(Nlh1+Nlh+Nlb1)*3,(Nlb)*3*sizeof(double),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hradsbuf+(Nlh),  m_drad  +(Nlh1+Nlh+Nlb1),(Nlb)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hrmasssbuf+(Nlh),m_drmass+(Nlh1+Nlh+Nlb1),(Nlb)*sizeof(double),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hinertsbuf+(Nlh),m_dinert+(Nlh1+Nlh+Nlb1),(Nlb)*sizeof(double),cudaMemcpyDeviceToHost)); 

/*if(rank==1)
for(i=0;i<Nlb;i++)
	printf("Nlb m_hradsbuf[%4d]=%10.6f,Nlb=%8d,Nlh=%8d\n",i,m_hradsbuf[Nlh+i],Nlh1+Nlh+Nlb1,Nlb,Nlh);*/


cpuErrchk(cudaMemcpy(m_hmatIdsbuf+(Nlh),m_dmatId+(Nlh1+Nlh+Nlb1),(Nlb)*sizeof(uint),cudaMemcpyDeviceToHost)); 

cpuErrchk(cudaMemcpy(m_hcontactEndsbuf+(Nlh)*3,m_dcontactEnd+(Nlh1+Nlh+Nlb1)*3,(Nlb)*3*sizeof(uint),cudaMemcpyDeviceToHost)); 
cpuErrchk(cudaMemcpy(m_hEngdispsbuf+(Nlh)*3,m_dEngdisp+(Nlh1+Nlh+Nlb1)*3,(Nlb)*3*sizeof(double),cudaMemcpyDeviceToHost)); 

MPI_CHECK(MPI_Send(m_hpossbuf +Nlh*3,Nlb*3,MPI_DOUBLE,left,44,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hpdissbuf+Nlh*3,Nlb*3,MPI_DOUBLE,left,45,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hangvsbuf+Nlh*3,Nlb*3,MPI_DOUBLE,left,46,MPI_COMM_WORLD));

MPI_CHECK(MPI_Send(m_hradsbuf  +Nlh,Nlb,MPI_DOUBLE,left,47,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hrmasssbuf+Nlh,Nlb,MPI_DOUBLE,left,48,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hinertsbuf+Nlh,Nlb,MPI_DOUBLE,left,49,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hmatIdsbuf+Nlh,Nlb,MPI_UNSIGNED,left,50,MPI_COMM_WORLD)); 

MPI_CHECK(MPI_Send(m_hcontactEndsbuf+Nlh*3,Nlb*3,MPI_UNSIGNED,left,151,MPI_COMM_WORLD)); 
MPI_CHECK(MPI_Send(m_hEngdispsbuf+Nlh*3,Nlb*3,MPI_DOUBLE,left,152,MPI_COMM_WORLD)); 
}

if(Nrh1 != 0)
{
MPI_CHECK(MPI_Recv(m_hposrbuf +(Nlh1+Nlb1+Nrb1)*3,Nrh1*3,MPI_DOUBLE,right,44,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hpdisrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3,MPI_DOUBLE,right,45,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hangvrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3,MPI_DOUBLE,right,46,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hradrbuf  +(Nlh1+Nlb1+Nrb1),Nrh1,MPI_DOUBLE,right,47,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hrmassrbuf+(Nlh1+Nlb1+Nrb1),Nrh1,MPI_DOUBLE,right,48,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hinertrbuf+(Nlh1+Nlb1+Nrb1),Nrh1,MPI_DOUBLE,right,49,MPI_COMM_WORLD,&status));

/*if(rank==0)
for(i=0;i<Nrh1;i++)
	printf("m_hradrbuf[%4d]=%10.6f,Nrh1=%8d,m_hPidrb+Nrh=%8d\n",
	i,m_hradrbuf[Nlh1+Nlb1+Nrb1+i],Nrh1,m_hPidrb+Nrh);*/

MPI_CHECK(MPI_Recv(m_hmatIdrbuf+(Nlh1+Nlb1+Nrb1),Nrh1,MPI_UNSIGNED,right,50,MPI_COMM_WORLD,&status));

MPI_CHECK(MPI_Recv(m_hcontactEndrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3,MPI_UNSIGNED,right,151,MPI_COMM_WORLD,&status));
MPI_CHECK(MPI_Recv(m_hEngdisprbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3,MPI_DOUBLE,right,152,MPI_COMM_WORLD,&status));

cpuErrchk(cudaMemcpy(m_dpos +(m_hPidrb+Nrh)*3,m_hposrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dpdis+(m_hPidrb+Nrh)*3,m_hpdisrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dangv+(m_hPidrb+Nrh)*3,m_hangvrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_drad  +(m_hPidrb+Nrh),m_hradrbuf+(Nlh1+Nlb1+Nrb1),Nrh1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_drmass+(m_hPidrb+Nrh),m_hrmassrbuf+(Nlh1+Nlb1+Nrb1),Nrh1*sizeof(double),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dinert+(m_hPidrb+Nrh),m_hinertrbuf+(Nlh1+Nlb1+Nrb1),Nrh1*sizeof(double),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dmatId+(m_hPidrb+Nrh),m_hmatIdrbuf+(Nlh1+Nlb1+Nrb1),Nrh1*sizeof(uint),cudaMemcpyHostToDevice));

cpuErrchk(cudaMemcpy(m_dcontactEnd+(m_hPidrb+Nrh)*3,m_hcontactEndrbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3*sizeof(uint),cudaMemcpyHostToDevice));
cpuErrchk(cudaMemcpy(m_dEngdisp+(m_hPidrb+Nrh)*3,m_hEngdisprbuf+(Nlh1+Nlb1+Nrb1)*3,Nrh1*3*sizeof(double),cudaMemcpyHostToDevice));
}

}