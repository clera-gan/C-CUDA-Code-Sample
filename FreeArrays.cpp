#include "dempacking.h"
//#include "FreeArrays.h"

void FreeParticleDeviceArray()
{
	cpuErrchk(cudaFree(m_dpos));
	cpuErrchk(cudaFree(m_dpdis));
	cpuErrchk(cudaFree(m_dangv));
	cpuErrchk(cudaFree(m_doldpos));
	cpuErrchk(cudaFree(m_doldpdis));
	cpuErrchk(cudaFree(m_doldangv));

	cpuErrchk(cudaFree(m_drad));
	cpuErrchk(cudaFree(m_drmass));
	cpuErrchk(cudaFree(m_dinert));
	cpuErrchk(cudaFree(m_dmatId));
	cpuErrchk(cudaFree(m_doldrad)); 
	cpuErrchk(cudaFree(m_doldrmass));
	cpuErrchk(cudaFree(m_doldinert));
	cpuErrchk(cudaFree(m_doldmatId));

    cpuErrchk(cudaFree(m_dParticleHash[0]));
	cpuErrchk(cudaFree(m_dParticleHash[1]));

	cpuErrchk(cudaFree(m_dAnei));
	cpuErrchk(cudaFree(m_dnjgi));
	cpuErrchk(cudaFree(m_dnjli));
	cpuErrchk(cudaFree(m_dsjgi));
	cpuErrchk(cudaFree(m_doldsjgi));
        cpuErrchk(cudaFree(m_dCN));
        cpuErrchk(cudaFree(m_dsumCN));

   cpuErrchk(cudaFree(m_ddisptw));
   cpuErrchk(cudaFree(m_dfricpw));
   cpuErrchk(cudaFree(m_dolddisptw));
   cpuErrchk(cudaFree(m_doldfricpw));

   cpuErrchk(cudaFree(m_dforcei));
   cpuErrchk(cudaFree(m_dtorquei));
   cpuErrchk(cudaFree(m_dcontacti));
   cpuErrchk(cudaFree(m_dqi));
   cpuErrchk(cudaFree(m_dkeng));
   cpuErrchk(cudaFree(m_doldIP)); //newly added

// add energy dissipition variables
   	cpuErrchk(cudaFree(m_dcontactEnd));
	cpuErrchk(cudaFree(m_dEngdisp));
	cpuErrchk(cudaFree(m_dEngdispw));
	cpuErrchk(cudaFree(m_dEngdispVarw));
	cpuErrchk(cudaFree(m_dBriProb));

	cpuErrchk(cudaFree(m_doldcontactEnd));
	cpuErrchk(cudaFree(m_doldEngdisp));
	cpuErrchk(cudaFree(m_doldEngdispw));
	cpuErrchk(cudaFree(m_doldEngdispVarw));
}

void FreeNeigListDeviceArrays()
{
	cpuErrchk(cudaFree(m_dLp));
	cpuErrchk(cudaFree(m_dLn));
	cpuErrchk(cudaFree(m_doldLn));

	cpuErrchk(cudaFree(m_ddispt));
	cpuErrchk(cudaFree(m_dfricp));
	cpuErrchk(cudaFree(m_dolddispt));
	cpuErrchk(cudaFree(m_doldfricp));

	cpuErrchk(cudaFree(m_dforcepair));
	cpuErrchk(cudaFree(m_dtorqpairi));
	cpuErrchk(cudaFree(m_dtorqpairj));
	cpuErrchk(cudaFree(m_dcontactpair));
	cpuErrchk(cudaFree(m_dpospairi));
	cpuErrchk(cudaFree(m_dpospairj));

	cpuErrchk(cudaFree(m_doldEngdisppair));
	cpuErrchk(cudaFree(m_doldcontactEndpair));

	cpuErrchk(cudaFree(m_dcontactEndpair));
	cpuErrchk(cudaFree(m_dEngdisppair));
	cpuErrchk(cudaFree(m_dEngdispVarpair));
}

void FreeParticleInCellDeviceArray()
{
   if(NallocateH>0)
   cpuErrchk(cudaFree(NumParticleInHCell));

   if(NallocateM>0)
	cpuErrchk(cudaFree(NumParticleInMCell));
	if(NallocateL>0)
	cpuErrchk(cudaFree(NumParticleInLCell));

}

void FreeGPUMeshData()
{
// free wall mesh data
cudaFree(m_dPosnode);
cudaFree(m_dFnodeid);

cudaFree(m_dsharePE);
cudaFree(m_dEdgeHead);
cudaFree(m_dsharePV);
cudaFree(m_dVertHead);

cudaFree(m_dFaceHash);
cudaFree(m_dFaceHashStart);
}


void FreeGPUData()
{
//-----release GPU memory -------------------------	
   FreeParticleDeviceArray();

   FreeParticleInCellDeviceArray();

   FreeGPUMeshData();

   cudaFree(m_dCellStart);
   cudaFree(Addcount);
} 


void FreeCPUParticleArrays()
{
//-------release CPU memory -------------------------	
	 free(m_hpos);
	 free(m_hpdis);
     free(m_hangv);

	 free(m_hrad);
     free(m_hinert);
	 free(m_hrmass);

	 free(m_hmatId);

	 free(m_holdsjgi);

        free(m_holddisptw);
        free(m_holdfricpw);

		free(m_hforcei);
		free(m_hkeng);
        free(m_hcontacti);


		free(m_hParticleHash);
		free(m_holdIP);

		free(m_hcontactEnd);
		free(m_hEngdisp);
		free(m_hEngdispw);
		free(m_hEngdispVarw);
//----------------------------------
} 

void FreeCPUNeigListArrays()
{
	 free(m_holdLn);
     free(m_holddispt);
     free(m_holdfricp);

	free(m_hforcepair);
	free(m_hpospairi);
	free(m_hpospairj);

	free(m_holdcontactEndpair);
	free(m_hEngdisppair);
	free(m_hEngdispVarpair);
}

void FreeCPUMeshData()
{
        free(m_hPosnode);
        free(m_hFnodeid);
}

void FreeCPUParticleBatchData()
{
        free(m_hradbatch);
        free(m_hrmassbatch);
        free(m_hinertbatch);
	    free(m_hpdisbatch);
        free(m_hangvbatch);
		free(m_hmatIdbatch);
}

void FreeMPIParticleBuffer()
{
        free(m_hParticleHashsbuf);
        free(m_hParticleHashrbuf);

        free(m_hpossbuf);
	    free(m_hposrbuf);
        free(m_hpdissbuf);
		free(m_hpdisrbuf);		
        free(m_hangvsbuf);
	    free(m_hangvrbuf);

        free(m_hradsbuf);
		free(m_hradrbuf);				
        free(m_hrmasssbuf);
	    free(m_hrmassrbuf);
        free(m_hinertsbuf);
		free(m_hinertrbuf);
				
        free(m_hmatIdsbuf);
	    free(m_hmatIdrbuf);

		free(m_hdisptwsbuf);
		free(m_hdisptwrbuf);
		free(m_hfricpwsbuf);
		free(m_hfricpwrbuf);

        free(m_hcontactEndsbuf);
	    free(m_hcontactEndrbuf);
        free(m_hEngdispsbuf);
		free(m_hEngdisprbuf);

		free(m_hEngdispwsbuf);
		free(m_hEngdispwrbuf);
		free(m_hEngdispVarwsbuf);
		free(m_hEngdispVarwrbuf);
}
/*
void FreeMPINeigListBuffer()
{
		free(m_hLnsbuf);
		free(m_hLnrbuf);
				
        free(m_hdisptsbuf);
	    free(m_hdisptrbuf);
        free(m_hfricpsbuf);
		free(m_hfricprbuf);

		free(m_hEngdispVarpairsbuf);
		free(m_hEngdispVarpairrbuf);
}*/

void FreeCPUFeedData()
{
free(m_radHmax);
free(m_radLmax);
free(m_radHmin);
free(m_radLmin);
}

void FreeCPUData()
{
FreeCPUParticleArrays();

FreeCPUNeigListArrays();

FreeCPUMeshData();

FreeCPUParticleBatchData();

FreeMPIParticleBuffer();

//FreeMPINeigListBuffer();
}

void JobFinish()
{
       if((fp1=fopen("jobfinished","wb"))==NULL)
		   printf("Cannot open file jobfinished!");
       fprintf(fp1,"job finished\n");
       fclose(fp1);
       if((fp1=fopen("restart.dat","wb"))==NULL)
		   printf("Cannot open file restart.dat");
       fprintf(fp1,"program comes to end!\n");
	   fprintf(fp1,"%8d %8d\n",1,1);
       printf("program comes to end!\n");
       fclose(fp1);
	   MPI_Finalize(); 
}
