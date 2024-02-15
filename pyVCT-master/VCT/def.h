// def.h
#ifndef _DEF
#define _DEF

// #define NULL 0
#define FALSE 0
#define TRUE 1
typedef int BOOL;
//#define SEED 2
extern int SEED;

#define PI 3.14159265359 

#define rounder(a) (((a) + ((a) < 0 ? 0.5 : -0.5) < (int)(a))? (int)(a): (int)(a) + 1)

//sample size
#define MULT 1
extern double VOXSIZE;		// [mm]
extern double sizeX;
extern double sizeY;
#define SCALE (VOXSIZE/.0025)	//		// [mm]
	// [mm]
#define sizeMarginX 0.100 		// [mm] from one side
#define sizeMarginY 0.100 		// [mm] from one side
#define MARGINX rounder(sizeMarginX/VOXSIZE)
#define MARGINY rounder(sizeMarginY/VOXSIZE)
extern int NVX; 
extern int NVY;
#define NV  (NVX*NVY)
//#define NRINC 901
extern int NRINC;
#define NRINC_CH 500

#define STEP_PRINT 100

#define MAXNRITER 1000
#define ACCURACY .00001

//parameters for channels distribution
#define IMMOTILITY_CH 1.0*SCALE*SCALE
#define JB	10.0
#define JH	2.0
#define G_NCH 50.0

extern double E_bond;  /* 5.0 */

// cells
#define IMMOTILITY 1.0*SCALE*SCALE						// 1/T
extern int NCX;
extern int NCY;

//fibers
extern double distanceF;
#define fiberD	0.0025 									// 0.0025 [mm] fibers diameter
#define F_DISTANCE rounder(distanceF/VOXSIZE)			// [pixels]
#define F_ANGLE 0										// angle of fibers with horizont

//spreading
extern char CONT;
extern char CONT_INHIB;
extern double GN_CM;
extern double GN_FB;
extern double PART; 			//% FBs

//elasticity
extern double TARGETVOLUME_CM;
extern double TARGETVOLUME_FB;
#define STARTVOLUME (TARGETVOLUME_FB/10)
extern double INELASTICITY_CM;
extern double INELASTICITY_FB;
extern double LMAX_CM;
extern double LMAX_FB;
#define INF 10000000.0

//nucleus protection
#define NUCLEI_R .007/VOXSIZE			// nucleus radius [pixels]
#define NUCL 2.0						// penalty for nucleus penetration (NUCL * DETACH)

//Js
extern double DETACH_CM;
extern double DETACH_FB;
#define JMDMD 0
extern double JCMMD;
extern double JFBMD;
extern double JCMCM;
extern double JFBFB;
extern double JFBCM;
extern double JCMCMc;
extern double JFBFBc;
extern double JFBCMc;

extern double UNLEASH_CM;
extern double UNLEASH_FB;

#define SQ05 .707107 					//sqrt(.5), used often enough to make this convenient

//number of focal adhesions
extern double MAX_FOCALS_CM;
extern double MAX_FOCALS_FB;

extern int silence;
extern int shifts;

#endif
