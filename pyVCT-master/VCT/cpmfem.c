#include "functions.h"
#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "libcpmfem.h"



typedef struct _THREAD_INFO
{
	int incr;
	//int startincr;
	int NRINC;
	VOX* pv;
	int NVX;
	int NVY;
	CM* CMs;
	int NRc;
	short* CCAlabels;
	FIBERS* pf;
	BOX* pb;
	int* attached;
	int* csize;
	double MAX_FOCALS_CM;
	double MAX_FOCALS_FB;
	double TARGETVOLUME_CM;
	double TARGETVOLUME_FB;
	double INELASTICITY_FB;
	double INELASTICITY_CM;
	double UNLEASH_CM;
	double UNLEASH_FB;
	double GN_CM;
	double GN_FB;
	double LMAX_CM;
	double LMAX_FB;
	double DETACH_CM;
	double DETACH_FB;
	double VOXSIZE;
	double JCMCM;
	double JCMMD;
	double JFBFB;
	double JFBMD;
	double JFBCM;
	double acceptance;
} THREAD_INFO;



void* CpmRoutine(void* Arg)
{
	THREAD_INFO* pThInfo = (THREAD_INFO*)Arg;
	
	if (pThInfo->incr % STEP_PRINT == 0){
			if(!silence)
				printf("\nSTART INCREMENT %d", pThInfo->incr);
			// добавлено (before writing was on each iteration)
			if(pThInfo->incr == (pThInfo->NRINC - 1))
			{
				write_cells(pThInfo->pv, pThInfo->incr, pThInfo->NVX, pThInfo->NVY);
				write_contacts(pThInfo->pv, pThInfo->incr, pThInfo->NVX, pThInfo->NVY);
			}
		}

		findCM(pThInfo->pv, pThInfo->CMs, pThInfo->NRc, pThInfo->NVX, pThInfo->NVY);
		
		pThInfo->acceptance = CPM_moves(
						pThInfo->pv,
						pThInfo->CCAlabels,
						pThInfo->pb,
						pThInfo->pf,
						pThInfo->CMs,
						pThInfo->attached,
						pThInfo->csize,
						pThInfo->MAX_FOCALS_CM,
						pThInfo->MAX_FOCALS_FB, 
						pThInfo->TARGETVOLUME_CM, 
						pThInfo->TARGETVOLUME_FB, 
						pThInfo->INELASTICITY_CM, 
						pThInfo->INELASTICITY_FB, 
						pThInfo->LMAX_CM, 
						pThInfo->LMAX_FB, 
						pThInfo->GN_CM, 
						pThInfo->GN_FB, 
						pThInfo->UNLEASH_CM, 
						pThInfo->UNLEASH_FB, 
						pThInfo->DETACH_CM, 
						pThInfo->DETACH_FB, 
						pThInfo->VOXSIZE, 
						pThInfo->NVX, 
						pThInfo->NVY, 
						pThInfo->JCMCM, 
						pThInfo->JCMMD, 
						pThInfo->JFBFB, 
						pThInfo->JFBMD, 
						pThInfo->JFBCM);

		if (pThInfo->incr % STEP_PRINT == 0 && !silence){
			printf("\nAcceptance rate %.4f", pThInfo->acceptance);
		}
		
		pthread_exit(NULL);
}


int cpmfem(
	int NCX, int NCY, 
	float* PART_matrix,
	double VOXSIZE,
	int NVX, int NVY,
	double GN_CM,
	double GN_FB,
	double TARGETVOLUME_CM,
	double TARGETVOLUME_FB,
	double DETACH_CM,
	double DETACH_FB,
	double INELASTICITY_FB,
	double INELASTICITY_CM,
	double JCMMD,
	double JFBMD,
	double JCMCM,
	double JFBFB,
	double JFBCM,
	double UNLEASH_CM,
	double UNLEASH_FB,
	double LMAX_CM,
	double LMAX_FB,
	double MAX_FOCALS_CM,
	double MAX_FOCALS_FB,
	int shifts,
	double distanceF,
	int SEED,
	int NRINC,
	int cyto,
	char* typ,
	int* cont_m,
	int* fibr,
	int* ctag_m
	
)
{

	//struct timeval tv;
    time_t time1;

	int d;
	VOX *pv;
	FIBERS *pf;
	int NRc,c,v;
	int * types;
	int *csize;
	int incr, startincr;
	double acceptance, acceptance_phi;
	
	if (cyto == 1){
		E_bond = 5.0;
		CONT = (char)1;
    	CONT_INHIB = (char)0;}

	if (cyto == 0){
		E_bond = 0.0;
		CONT = (char)0;
    	CONT_INHIB = (char)1;}

	
	if(!silence){
		printf("SEED = %d\n",SEED);
		printf("Sample size = %d x %d\n",NCX,NCY);
		printf("\n");
		printf("GN_CM = %.2f\n",GN_CM*SCALE);
		printf("TARGETVOLUME_CM = %.2f\n",TARGETVOLUME_CM*1000*VOXSIZE*VOXSIZE);
		printf("INELASTICITY_CM = %.2f\n",INELASTICITY_CM/(SCALE*SCALE*SCALE*SCALE));
		printf("DETACH_CM = %.2f\n",DETACH_CM*SCALE);
		printf("\n");
		printf("GN_FB = %.2f\n",GN_FB*SCALE);
		printf("TARGETVOLUME_FB = %.2f\n",TARGETVOLUME_FB*1000*VOXSIZE*VOXSIZE);
		printf("INELASTICITY_FB = %.2f\n",INELASTICITY_FB/(SCALE*SCALE*SCALE*SCALE));
		printf("DETACH_FB = %.2f\n",DETACH_FB*SCALE);
		printf("\n");
		printf("JMDMD = %.3f\n",JMDMD/VOXSIZE);
		printf("JCMMD = %.3f\n",JCMMD/VOXSIZE);
		printf("JFBMD = %.3f\n",JFBMD/VOXSIZE);
		printf("\n");
		printf("JCMCM = %.3f\n",JCMCM/VOXSIZE);
		printf("JFBFB = %.3f\n",JFBFB/VOXSIZE);
		printf("JFBCM = %.3f\n",JFBCM/VOXSIZE);
		printf("\n");
		printf("UNLEASH_CM = %.2f\n",UNLEASH_CM*SCALE);
		printf("UNLEASH_FB = %.2f\n",UNLEASH_FB*SCALE);
		printf("\n");
		printf("LMAX_CM = %.2f px\n",LMAX_CM);
		printf("LMAX_FB = %.2f px\n",LMAX_FB);
		printf("\n");
		printf("MAX_FOCALS_CM = %.2f\n",MAX_FOCALS_CM);
		printf("MAX_FOCALS_FB = %.2f\n",MAX_FOCALS_FB);
	}

	/// INITIALIZE ///
   	srand(SEED); mt_init();
   	pv = init_voxels(NVX, NVY);
	pf = set_fibers(distanceF, VOXSIZE, NVX, NVY);
	BOX * pb = allocBOX(NCX*NCY+1);

	write_fibers(pf, NVX, NVY);

	startincr = 0;
	types = calloc((NCX*NCY+1), sizeof(int));
	NRc = init_cells(pv,types,pb,NCX,NCY, PART_matrix, shifts, TARGETVOLUME_FB, VOXSIZE, NVX, NVY);
	write_cells(pv,0, NVX, NVY);
	csize = calloc(NRc, sizeof(int)); for(c=0;c<NRc;c++) {csize[c]=0;}
	for(v=0;v<NV;v++) {if(pv[v].ctag) {csize[pv[v].ctag-1]++;}}

	CM* CMs = allocCM(NRc);
	int* attached = alloc_attach(NRc);

	short * CCAlabels = malloc(NV * sizeof(short));

	time1 = time(NULL);
	//time = tv.tv_sec;  

	write_types(types,NRc);		//save types into file

	// START SIMULATION ///
	// Naive multithreading
	pthread_t* threadsCPM = (pthread_t*)malloc((NRINC - startincr)*sizeof(pthread_t));
	
	THREAD_INFO thInfo = {0};
	
	thInfo.incr = incr;
	//thInfo.startincr = startincr;
	thInfo.NRINC = NRINC;
	thInfo.pv = pv;
	thInfo.NVX = NVX;
	thInfo.NVY = NVY;
	thInfo.CMs = CMs;
	thInfo.NRc = NRc;
	thInfo.CCAlabels = CCAlabels;
	thInfo.pf = pf;
	thInfo.pb = pb;
	thInfo.attached = attached;
	thInfo.csize = csize;
	thInfo.MAX_FOCALS_CM = MAX_FOCALS_CM;
	thInfo.MAX_FOCALS_FB = MAX_FOCALS_FB;
	thInfo.TARGETVOLUME_CM = TARGETVOLUME_CM;
	thInfo.TARGETVOLUME_FB = TARGETVOLUME_FB;
	thInfo.INELASTICITY_FB = INELASTICITY_FB;
	thInfo.INELASTICITY_CM = INELASTICITY_CM;
	thInfo.UNLEASH_CM = UNLEASH_CM;
	thInfo.UNLEASH_FB = UNLEASH_FB;
	thInfo.GN_CM = GN_CM;
	thInfo.GN_FB = GN_FB;
	thInfo.LMAX_CM = LMAX_CM;
	thInfo.LMAX_FB = LMAX_FB;
	thInfo.DETACH_CM = DETACH_CM;
	thInfo.DETACH_FB = DETACH_FB;
	thInfo.VOXSIZE = VOXSIZE;
	thInfo.JCMCM = JCMCM;
	thInfo.JCMMD = JCMMD;
	thInfo.JFBFB = JFBFB;
	thInfo.JFBMD = JFBMD;
	thInfo.JFBCM = JFBCM;
	thInfo.acceptance = acceptance;
	
	for(incr=startincr; incr<NRINC; incr++)
	{
		pthread_create(&threadsCPM[incr], NULL, CpmRoutine, &thInfo);	
	}

	/// END ///
	if(!silence)
	printf("\nSIMULATION FINISHED!\n");
	
	
	for(incr=startincr; incr<NRINC; incr++)
	{
		pthread_join(threadsCPM[incr], NULL);
	}
	
	free(threadsCPM);

	write_contacts(pv,0, NVX, NVY);

	/*pv = init_voxels();
	read_cells(pv,types, NRc, "./output/ctags1.sout","./output/conts1.sout","./output/types.sout");*/

	/// START DISTRIBUTION ///
	findCM(pv,CMs,NRc, NVX, NVY);
	for(incr=startincr; incr<NRINC_CH; incr++)
	{
		if (incr % 100 == 0){
			if(!silence)
				printf("\nSTART CHANNEL DISTRIBUTION %d",incr);
			write_contacts(pv,incr+1, NVX, NVY);
		}

		acceptance = CH_moves(pv, CMs, 0.5 + 0.5*incr/NRINC, VOXSIZE, NVX, NVY);

		if (incr % 100 == 0 && !silence){
			printf("\nAcceptance rate %.4f",acceptance);
		}
	}

	write_cells(pv,1, NVX, NVY);
	write_contacts(pv,1, NVX, NVY);

	/// END ///
	if(!silence)
	printf("\nSIMULATION FINISHED!\n");

	time_t time2 = time(NULL);
	if(!silence)
	printf("Took %lds\n", time2 - time1);
	
	//fill our outer massives//
	int i,j,k;
	for(i=0; i<NVX; i++) {
    	for (j=0; j<NVY; j++) {
        	k = i + j * NVX;
   			ctag_m[k]=pv[k].ctag;
			cont_m[k]=pv[k].contact;
			fibr[k]=(char)pf[k].Q;
		}
    }
	int l;
	for(l=0; l<NRc; l++) {
      typ[l] = types[1+l];
    }

	printf("\nmassives done!\n");

	free(pv); 
	free(pf);
	free(CCAlabels);
	
	return 0;
}

