/***
 *
 * Version:                             2.0
 *
 * Package:                             DANSS SiPm Signal Processing and Calibration
 *
 * Description:                         Example user functions to process digitized SiPm data
 *
 *	Calculate different event parameters and put to root file
 *
 ***/

#include <stdio.h>
#include <string.h>

#include "Riostream.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TTreeCacheUnzip.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TKey.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
//#include "TSpectrum.h"

#include "readDigiData.h"
#include "danssScmGlobals.h"
#include "evtbuilder.h"

/***********************	Definitions	****************************/
//	Veto
#define VETO2MEV	0.002			// Translate corrected VETO signal to MeV
//	Initial clean parameters
#define MINSIPMPIXELS	3			// Minimum number of pixels to consider SiPM hit
#define MINSIPMPIXELS2	2			// Minimum number of pixels to consider SiPM hit without confirmation (method 2)
#define MINPMTENERGY	0.1			// Minimum PMT energy for a hit
#define MINVETOENERGY	0.1			// Minimum VETO energy for a hit
#define SIPMEARLYTIME	45			// ns - shift from fine time
#define SOMEEARLYTIME	130			// ns - absolute if fineTime is not defined
//	fine time
#define MINENERGY4TIME	0.25			// Minimum energy to use for fine time averaging
#define TCUT		15			// fine time cut, ns
#define NOFINETIME	10000			// something out of range
//	Positron
#define ATTENUATION	0.00342			// Signal attenuation for positron energy correction
#define MCATTENUATION	0.005			// Signal attenuation for positron energy correction MC
//	Flags
#define FLG_PRINTALL		1
#define FLG_NOCLEANNOISE 	0x10000
#define FLG_NOTIMECUT		0x20000
#define FLG_NOCONFIRM		0x40000
#define FLG_NOCONFIRM2		0x80000

using namespace std;

// Globals:

long long                               iNevtTotal;
long long				upTime;
long long				fileFirstTime;
long long				fileLastTime;
int                                     progStartTime;
char *					chTimeCalibration;
char *					chOutputFile;
int					iFlags;
int					MaxEvents;
int					IsMc;				// MC run flag


TFile *					OutputFile;
TTree *					OutputTree;
TTree *					InfoTree;
struct DanssEventStruct2		DanssEvent;
struct DanssInfoStruct			DanssInfo;
struct DanssMcStruct			DanssMc;
int 					HitFlag[iMaxDataElements];	// array to flag out SiPM hits

/************************	Analysis functions					*****************************/
//				Auxillary
float VetoEnergy(int chan, float amp)
{
//	We use here on-line coefficiencies
        const float TrigCoef[64] = {
		0.417, 1.429, 0.714, 1.786, 0.714, 0.526, 0.714, 1.316, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
		0.152, 0.294, 0.185, 0.217, 0.227, 0.083, 0.192, 0.263, 0.152, 0.161, 0.139, 0.227, 0.185, 0.156, 0.200, 0.156,
		0.091, 0.179, 0.150, 0.083, 0.172, 0.152, 0.156, 0.200, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
		0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
	return amp * TrigCoef[chan] * VETO2MEV;
}

int IsInModule(int SiPm, int Pmt, ReadDigiDataUser *user)
{
	int SiPmXY, PmtXY;
	int SiPmZ, PmtZ;

	if (user->side(SiPm) != user->side(Pmt)) return false;
	if (user->type(SiPm) != SiPmHit || user->type(Pmt) != PmtHit) return false;	// wrong request
	SiPmXY = user->firstCoord(SiPm);
	PmtXY  = user->firstCoord(Pmt);
	SiPmZ  = user->zCoord(SiPm);
	PmtZ   = user->zCoord(Pmt);
	if (SiPmXY / 5 != PmtXY || SiPmZ / 20 != PmtZ) return false;
	return true;
}

int IsNeighbor(int hitA, int hitB, ReadDigiDataUser *user)
{
	if (user->zCoord(hitA) == user->zCoord(hitB) && abs(user->firstCoord(hitA) - user->firstCoord(hitB)) <= 1) return 1;
	if (abs(user->zCoord(hitA) - user->zCoord(hitB)) == 1) return 1;
	return 0;
}

float acorr(float energy, float dist, int side)
{
	float C;

	if (dist >= 0) {
		if (IsMc) {
			C = exp(MCATTENUATION * (dist - 50.0));	// 50 cm is the middle
		 	if (side) C = 1/C;			// side = 1 - X
		} else {
			C = exp(ATTENUATION * (dist - 50.0));	// 50 cm is the middle
		}
	} else {
		C = 1;
	}
	return C * energy;
}

//				Main analysis

void CalculatePositron(ReadDigiDataUser *user)
{
	int i, j, N;
	float A;
	float x, y, z;
	float nx, ny;
	int maxHit;
	int clHits;
	int repeat;

	N = user->nhits();
//		Find the maximum hit
	A = 0;
	maxHit = -1;
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && user->type(i) == SiPmHit && user->e(i) > A) {
		A = user->e(i);
		maxHit = i;
	}
	if (maxHit < 0) return;	// nothing to do - no usable SiPM hits
	HitFlag[maxHit] = 10;
	DanssEvent.MaxHitEnergy = A;
//		Find cluster
	for (;;) {
		repeat = 0;
		for (i=0; i<N; i++) if (HitFlag[i] >= 10) for (j=0; j<N; j++) if (HitFlag[j] >= 0 && HitFlag[j] < 10 && user->type(j) == SiPmHit && IsNeighbor(i, j, user)) {
			HitFlag[j] = 20;
			repeat = 1;
		}
		if (!repeat) break;
	}
//		Find cluster position
	x = y = z = 0;
	nx = ny = 0;
	for (i=0; i<N; i++) if (HitFlag[i] >= 10) {
		DanssEvent.PositronHits++;
		if (user->side(i) == 'X') {
			x += user->firstCoord(i) * fStripWidth * user->e(i);
			z += user->zCoord(i) * fStripHeight * user->e(i);
			nx += user->e(i);
		} else {
			y += user->firstCoord(i) * fStripWidth * user->e(i);
			z += user->zCoord(i) * fStripHeight * user->e(i);
			ny += user->e(i);
		}
	}
	DanssEvent.PositronX[0] = (nx > 0) ? x / nx : -1;		// 50 cm is DANSS center
	DanssEvent.PositronX[1] = (ny > 0) ? y / ny : -1;		// 50 cm is DANSS center
	DanssEvent.PositronX[2] = (nx + ny > 0) ? z / (nx + ny) : -1;	// 50 cm is DANSS center
//		Find corrected energy
	for (i=0; i<N; i++) if (HitFlag[i] >= 10) {
		if (user->side(i) == 'X') {
			DanssEvent.PositronSiPmEnergy += acorr(user->e(i), DanssEvent.PositronX[1], 1);
		} else {
			DanssEvent.PositronSiPmEnergy += acorr(user->e(i), DanssEvent.PositronX[0], 0);
		}
	}
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && user->type(i) == PmtHit) {
		for (j=0; j<N; j++) if (IsInModule(j, i, user) && HitFlag[j] >= 10) break;
		if (j >= N) continue;
		if (user->side(i) == 'X') {
			DanssEvent.PositronPmtEnergy += acorr(user->e(i), DanssEvent.PositronX[1], 1);
		} else {
			DanssEvent.PositronPmtEnergy += acorr(user->e(i), DanssEvent.PositronX[0], 0);
		}
	}
//		Count possible gammas
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && HitFlag[i] < 10 && user->type(i) == SiPmHit) {
		DanssEvent.AnnihilationGammas++;
		DanssEvent.AnnihilationEnergy += user->e(i);
	}
}

void CalculateNeutron(ReadDigiDataUser *user)
{
	float x, y, z, r;
	int nx, ny;
	int i, N;

	N = user->nhits();
//	Find the center (1st approximation)
	x = y = z = 0;
	nx = ny = 0;
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && user->type(i) == SiPmHit) {
		if (user->side(i) == 'X') {
			x += user->firstCoord(i) * fStripWidth;
			z += user->zCoord(i) * fStripHeight;
			nx++;
		} else {
			y += user->firstCoord(i) * fStripWidth;
			z += user->zCoord(i) * fStripHeight;
			ny++;
		}
	}
	DanssEvent.NeutronX[0] = (nx) ? x / nx : -1;			// 50 cm is DANSS center
	DanssEvent.NeutronX[1] = (ny) ? y / ny : -1;			// 50 cm is DANSS center
	DanssEvent.NeutronX[2] = (nx + ny) ? z / (nx + ny) : -1;	// 50 cm is DANSS center
//	Average distance
	nx = 0;
	for (i=0; i<N; i++) if (user->type(i) == SiPmHit && HitFlag[i] >= 0) {
		r = DanssEvent.NeutronX[(user->side(i) == 'X') ? 0 : 1]; 
		r = (user->firstCoord(i) * fStripWidth - r) * (user->firstCoord(i) * fStripWidth - r);
		r += (user->zCoord(i) * fStripHeight - DanssEvent.NeutronX[2]) * (user->zCoord(i) * fStripHeight - DanssEvent.NeutronX[2]);
		r = sqrt(r);
		DanssEvent.NeutronRadius += r;
		nx++;
	}
	DanssEvent.NeutronRadius = (nx) ? DanssEvent.NeutronRadius / nx : -1;
}

void CleanZeroes(ReadDigiDataUser *user)
{
	int i, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) if (user->type(i) == SiPmHit && user->npix(i) <= 0) HitFlag[i] = -1;
}

void CleanNoise(ReadDigiDataUser *user)
{
	int i, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) switch (user->type(i)) {
	case SiPmHit:
		if (user->npix(i) < MINSIPMPIXELS) HitFlag[i] = -1;
		break;
	case PmtHit:
		if (user->e(i) < MINPMTENERGY) HitFlag[i] = -1;
		break;
	case VetoHit:
		if (VetoEnergy(user->adcChan(i), user->signal(i)) < MINVETOENERGY) HitFlag[i] = -1;
		break;
	}
}

void CleanByConfirmation(ReadDigiDataUser *user)
{
	int i, j, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) if (HitFlag[i] >= 0) switch (user->type(i)) {
	case SiPmHit:
		for (j=0; j<N; j++) if (HitFlag[j] >= 0 && user->type(j) == PmtHit && IsInModule(i, j, user)) break;
		if (j == N) HitFlag[i] = -1;
		break;
	case PmtHit:
		for (j=0; j<N; j++) if (HitFlag[j] >= 0 && user->type(j) == SiPmHit && IsInModule(j, i, user)) break;
		if (j == N) HitFlag[i] = -1;
		break;
	}
}


/*	Clean SiPM only if npix == 1 and no PMT confirmation	*/
void CleanByConfirmation2(ReadDigiDataUser *user)
{
	int i, j, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && user->type(i) == SiPmHit) {
		if (user->npix(i) >= MINSIPMPIXELS2) continue;		// that's enough
		for (j=0; j<N; j++) if (HitFlag[j] >= 0 && user->type(j) == PmtHit && IsInModule(i, j, user)) break;
		if (j < N) continue;
		HitFlag[i] = -1;
	}
//		"early" hits
	for (i=0; i<N; i++) if (HitFlag[i] == -100)
	{
		if (user->npix(i) >= MINSIPMPIXELS2) continue;		// that's enough
		for (j=0; j<N; j++) if (HitFlag[j] >= 0 && user->type(j) == PmtHit && IsInModule(i, j, user)) break;
		if (j < N) continue;
		HitFlag[i] = -1;
	}
}

void CleanByTime(ReadDigiDataUser *user)
{
	int i, N;
	float tearly;
	
	N = user->nhits();
	if (DanssEvent.fineTime != NOFINETIME) {
		for (i=0; i<N; i++) if (fabs(user->t(i) - DanssEvent.fineTime) > TCUT) HitFlag[i] = -1;
		tearly = DanssEvent.fineTime - SIPMEARLYTIME;
	} else {
		tearly = SOMEEARLYTIME;
	}
	for (i=0; i<N; i++) if (HitFlag[i] >= 0 && user->type(i) == SiPmHit && fabs(user->t(i) - tearly) <= TCUT) HitFlag[i] = -100;	// mark early hit candidates
}

void DebugFullPrint(ReadDigiDataUser *user)
{
	int i, N;
	time_t tm;

	N = user->nhits();
	tm = DanssEvent.unixTime;
	printf("******************************************************************************************************************\n");
	printf("Event: %Ld globalTime: %Ld fineTime: %6.1f ns   linux time: %s", 
		DanssEvent.number, DanssEvent.globalTime, DanssEvent.fineTime, ctime(&tm));
	printf("Total %d hits: %d SiPM %d PMT %d Veto; Clean: %d SiPM %d PMT %d Veto\n", 
		N, DanssEvent.SiPmHits, DanssEvent.PmtHits, DanssEvent.VetoHits, 
		DanssEvent.SiPmCleanHits, DanssEvent.PmtCleanHits, DanssEvent.VetoCleanHits);
	printf("Energy: %6.1f SiPM %6.1f PMT %6.1f Veto; Clean: %6.1f SiPM %6.1f PMT %6.1f Veto\n", 
		DanssEvent.SiPmEnergy, DanssEvent.PmtEnergy, DanssEvent.VetoEnergy, 
		DanssEvent.SiPmCleanEnergy, DanssEvent.PmtCleanEnergy, DanssEvent.VetoCleanEnergy);
	if (N) {
		printf("N    Type  N  S       E    time  ADC.Ch side XY Z  Flag\n");
//			1234512345123412345678123451234561231123123451231231234
		for(i=0; i<N; i++) switch(user->type(i)) {
		case SiPmHit:
			printf("%4d SiPM %3.0f %7.1f %4.1f %5.1f %2d.%2.2d    %c  %2d %2d  %c\n", i+1, user->npix(i), user->signal(i),
				user->e(i), user->adc(i), user->t(i), user->adcChan(i), user->side(i), user->firstCoord(i), user->zCoord(i),
				(HitFlag[i]<0) ? 'X' : ' ');
			break;
		case PmtHit:
			printf("%4d PMT      %7.1f %4.1f %5.1f %2d.%2.2d    %c  %2d %2d  %c\n", i+1, user->signal(i),
				user->e(i), user->t(i), user->adc(i), user->adcChan(i), user->side(i), user->firstCoord(i), user->zCoord(i),
				(HitFlag[i]<0) ? 'X' : ' ');
			break;
		case VetoHit:
			printf("%4d VETO     %7.1f %4.1f %5.1f %2d.%2.2d    -  xx xx  %c\n", i+1, user->signal(i),
				VetoEnergy(user->adcChan(i), user->signal(i)), user->t(i), user->adc(i), user->adcChan(i),
				(HitFlag[i]<0) ? 'X' : ' ');
			break;
		}
	}
}

void FindFineTime(ReadDigiDataUser *user)
{
	float tsum;
	float asum;
	float e;
	int i, N;
	
	tsum = asum = 0;
	N = user->nhits();
	for (i=0; i<N; i++) if (HitFlag[i] >= 0) {
		switch(user->type(i)) {
		case SiPmHit:
			e = user->e(i);
			if (user->npix(i) < MINSIPMPIXELS) e = 0;
			break;
		case PmtHit:
			e = user->e(i);
			break;
		case VetoHit:
			e = VetoEnergy(user->adcChan(i), user->signal(i));
			break;
		}
		if (e > MINENERGY4TIME) {
			tsum += user->t(i) * e;
			asum += e;
		}
	}
	DanssEvent.fineTime = (asum > 0) ? tsum / asum : NOFINETIME;	// some large number if not usable hits found
}

void SumClean(ReadDigiDataUser *user)
{
	int i, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) if (HitFlag[i] >= 0) switch (user->type(i)) {
	case SiPmHit:
		DanssEvent.SiPmCleanHits++;
		DanssEvent.SiPmCleanEnergy += user->e(i);
		break;
	case PmtHit:
		DanssEvent.PmtCleanHits++;
		DanssEvent.PmtCleanEnergy += user->e(i);
		break;
	case VetoHit:
		DanssEvent.VetoCleanHits++;
		DanssEvent.VetoCleanEnergy += VetoEnergy(user->adcChan(i), user->signal(i));
		break;
	}

	for (i=0; i<N; i++) if (HitFlag[i] == -100) {
		DanssEvent.SiPmEarlyHits++;
		DanssEvent.SiPmEarlyEnergy += user->e(i);		
	}
}

void SumEverything(ReadDigiDataUser *user)
{
	int i, N;
	
	N = user->nhits();
	for (i=0; i<N; i++) switch (user->type(i)) {
	case SiPmHit:
		DanssEvent.SiPmHits++;
		DanssEvent.SiPmEnergy += user->e(i);
		break;
	case PmtHit:
		DanssEvent.PmtHits++;
		DanssEvent.PmtEnergy += user->e(i);
		break;
	case VetoHit:
		DanssEvent.VetoHits++;
		DanssEvent.VetoEnergy += VetoEnergy(user->adcChan(i), user->signal(i));
		break;
	}
}

/************************	class ReadDigiDataUser user functions			*****************************/

/***
 *
 * A function to initialize calibration delays
 *
 * Input parameters: -
 *
 * Return value: -
 *
 ***/

void ReadDigiDataUser::init_Tds()
{
	FILE *f;
	char str[1024];
	char *ptr;
	int i, k;
	int iAdcNum, iAdcChan;
	

//	Set all zeroes
	for(i = 0; i < iNElements; i++) {
    		iAdcNum = i / 100;
    		iAdcChan = i % 100;
    		if(!isAdcChannelExist(iAdcNum, iAdcChan)) continue;
    		setTd(i, 0); // set all td = 0
  	}
	if (!chTimeCalibration) return;

//	Read and implement tcalib file
	f = fopen(chTimeCalibration, "rt");
	if (!f) {
		printf("Can not open file %s: %m\n", chTimeCalibration);
		return;
	}
	k = 0;
	for(;;) {
		ptr = fgets(str, sizeof(str), f);
		if (!ptr) break;		// EOF or error
		if (str[0] != 'C') {
			printf("Time calibration: %s", str);
			continue;	// Comment ?
		}
		ptr = strstr(str, "Channel=");
		if (!ptr) {
			printf("Time calibration no Channel=: %s", str);
			continue;		// strange string
		}
		ptr += strlen("Channel=");
		i = 100 * (strtod(ptr, NULL) + 0.002);
    		iAdcNum = i / 100;
    		iAdcChan = i % 100;
    		if(!isAdcChannelExist(iAdcNum, iAdcChan)) {	// non-existing channel - strange
			printf("Time calibration wrong channel i=%d (%d.%d): %s", i, iAdcNum, iAdcChan, str);
			continue;
		}
		ptr = strstr(str, "DT=");
		if (!ptr) {		// strange string	
			printf("Time calibration no value DT=: %s", str);
			continue;
		}
		ptr += strlen("DT=");
		setTd(i, strtod(ptr, NULL));
		k++;
	}

	printf("Time calibration used: %s. %d channels found.\n", chTimeCalibration, k);

	fclose(f);	
}

//------------------------------->


/***
 *
 * A function which is called to initialized user data
 *
 * Input parameters: standard command line parameters argc, argv
 *
 * Return value: -
 *
 ***/

void ReadDigiDataUser::initUserData(int argc, const char **argv)
{
	int i, j;
	char strs[128];
	char strl[1024];

	progStartTime = time(NULL);
	chOutputFile = NULL;
	chTimeCalibration = NULL;
	iFlags = 0;
	MaxEvents = -1;
	IsMc = 0;

	for (i=1; i<argc; i++) {
		if (!strcmp(argv[i], "-output")) {
			i++;
			chOutputFile = (char *)argv[i];
		} else if (!strcmp(argv[i], "-tcalib")) {
			i++;
			chTimeCalibration = (char *)argv[i];
		} else if (!strcmp(argv[i], "-flag")) {
			i++;
			iFlags = strtol(argv[i], NULL, 0);
		} else if (!strcmp(argv[i], "-events")) {
			i++;
			MaxEvents = strtol(argv[i], NULL, 0);
		} else if (!strcmp(argv[i], "-mcdata")) {
			IsMc = 1;
		}
	}

	if (chTimeCalibration) init_Tds();

	if (!chOutputFile) {
		chOutputFile = (char *) malloc(strlen(argv[0]) + 6);
		sprintf(chOutputFile, "%s%s", argv[0], ".root");
	}
	
	OutputFile = new TFile(chOutputFile, "RECREATE");
	if (!OutputFile->IsOpen()) throw "Panic - can not open output file!";
	OutputTree = new TTree("DanssEvent", "Danss event tree");
	OutputTree->Branch("Data", &DanssEvent, 
//		Common parameters
		"globalTime/L:"		// time in terms of 125 MHz
		"number/L:"		// event number in the file
		"unixTime/I:"		// linux time, seconds
		"fineTime/F:"		// fine time of the event (for hit selection)
//		Veto parameters
		"VetoHits/I:"		// hits in the event record
		"VetoEnergy/F:"		// Energy Sum of all hits
		"VetoCleanHits/I:"	// hits above threshold and in time window
		"VetoCleanEnergy/F:"	// Energy Sum of clean hits
//		PMT parameters
		"PmtHits/I:"		// the same as above for PMT
		"PmtEnergy/F:"
		"PmtCleanHits/I:"
		"PmtCleanEnergy/F:"
//		SiPM parameters
		"SiPmHits/I:"		// the same as above for PMT
		"SiPmEnergy/F:"
		"SiPmCleanHits/I:"
		"SiPmCleanEnergy/F:"
		"SiPmEarlyHits/I:"
		"SiPmEarlyEnergy/F:"		
//		"positron cluster" parameters
		"PositronHits/I:"	// hits in the cluster
		"PositronSiPmEnergy/F:"	// Energy sum of the cluster, corrected (SiPM)
		"PositronPmtEnergy/F:"	// Energy sum of the cluster, corrected (PMT)
		"MaxHitEnergy/F:"	// Energy of the maximum hit
		"PositronX[3]/F:"	// cluster position
		"AnnihilationGammas/I:"	// number of possible annihilation gammas
		"AnnihilationEnergy/F:"	// Energy in annihilation gammas
//		"neutron" parameters
		"NeutronX[3]/F:"	// center of gammas position
		"NeutronRadius/F"	// average distance between hits and the center
	);
	if (IsMc) OutputTree->Branch("MC", &DanssMc,
		"McEnergy/F:"		// MC true energy
		"McX[3]/F"		// MC vertex position
	);

	InfoTree = new TTree("DanssInfo", "Run info tree");	
	InfoTree->Branch("Info", &DanssInfo,  
		"gTime/L:"		// running time in terms of 125 MHz
		"startTime/I:"		// linux start time, seconds
		"stopTime/I:"		// linux stop time, seconds
		"events/I"		// number of events
	);
	iNevtTotal = 0;
	upTime = 0;
	fileFirstTime = -1;
	fileLastTime = -1;
	memset(&DanssInfo, 0, sizeof(struct DanssInfoStruct));
}

//------------------------------->

/***
 *
 * A function which is called once per each event
 *
 * Input parameters: -
 *
 * Return value: - 0 - OK, -1 - stop
 *
 ***/

int ReadDigiDataUser::processUserEvent()
{
	float fineTime;

  	if( ttype() != 1 ) return 0;
	
	memset(HitFlag, 0, nhits() * sizeof(int));
	memset(&DanssEvent, 0, sizeof(struct DanssEventStruct2));

	fileLastTime = globalTime();
	DanssInfo.stopTime = absTime();
	if (fileFirstTime < 0) {
		fileFirstTime = fileLastTime;
		DanssInfo.startTime = absTime();
		DanssInfo.events = 0;
	}
  	iNevtTotal++;
	DanssInfo.events++;

	DanssEvent.globalTime = globalTime();
	DanssEvent.number     = nevt();
	DanssEvent.unixTime   = absTime();

	CleanZeroes(this);
	SumEverything(this);
	if (!(iFlags & FLG_NOCLEANNOISE)) CleanNoise(this);
	FindFineTime(this);
	if (!(iFlags & FLG_NOTIMECUT)) CleanByTime(this);
	if (!(iFlags & FLG_NOCONFIRM)) CleanByConfirmation(this);
	if (!(iFlags & FLG_NOCONFIRM2)) CleanByConfirmation2(this);
	SumClean(this);
	CalculateNeutron(this);
	CalculatePositron(this);
	if (iFlags & FLG_PRINTALL) DebugFullPrint(this);

	if (IsMc) mcTruth(DanssMc.Energy, DanssMc.X[0], DanssMc.X[1], DanssMc.X[2]);

	OutputTree->Fill();

	if (MaxEvents > 0 && iNevtTotal >= MaxEvents) return -1;
  	return 0;
}

//------------------------------->

/***
 *
 * A function which is called to finish user data processing
 *
 * Input parameters: -
 *
 * Return value: -
 *
 ***/

void ReadDigiDataUser::finishUserProc()
{
	if (fileFirstTime > 0) {
		DanssInfo.upTime = fileLastTime - fileFirstTime;
		InfoTree->Fill();
	}
	OutputTree->Write();
	InfoTree->Write();	
  	OutputFile->Close();
  
	printf("Total up time        %Ld seconds\n", upTime / (long long) GLOBALFREQ);
	printf("Total physics events %Ld\n", iNevtTotal);
	printf("Run completed in     %d seconds\n", time(NULL) - progStartTime);
  	return;
}

//--------------------->


int ReadDigiDataUser::userActionAtFileChange()
{
	DanssInfo.upTime = fileLastTime - fileFirstTime;
	InfoTree->Fill();
	upTime += DanssInfo.upTime;
	fileFirstTime = -1;
	return 0;
}

