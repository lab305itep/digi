/***
 *
 * Version:                             2.0
 *
 * Package:                             DANSS SiPm Signal Processing and Calibration
 *
 * Description:                         Example user functions to process digitized SiPm data
 *
 *	Get times for analysis
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
#include "TSpectrum.h"

#include "readDigiData.h"
#include "danssScmGlobals.h"

#define EHITMIN		1.0
#define TOTALMIN	20.0
#define NSIPMMIN	10
#define GLOBALFREQ	125000000

using namespace std;

// Globals:

long long                               iNevtTotal;
long long				upTime;
long long				fileFirstTime;
long long				fileLastTime;
int                                     progStartTime;
char *					chTimeCalibration;
char *					chOutputFile;

// Here one can list histograms:

TH1D *                                  hDt[iN_AdcBoards][iNChannels_AdcBoard];

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

	chOutputFile = NULL;
	chTimeCalibration = NULL;

	for (i=1; i<argc; i++) {
		if (!strcmp(argv[i], "-output")) {
			i++;
			chOutputFile = (char *)argv[i];
		} else if (!strcmp(argv[i], "-tcalib")) {
			i++;
			chTimeCalibration = (char *)argv[i];
		}
	}

	if (chTimeCalibration) init_Tds();

	if (!chOutputFile) {
		chOutputFile = (char *) malloc(strlen(argv[0]) + 6);
		sprintf(chOutputFile, "%s%s", argv[0], ".root");
	}

	progStartTime = time(NULL);
	for (i=0; i<iN_AdcBoards; i++) for (j=0; j<iNChannels_AdcBoard; j++) {
		sprintf(strs, "hDt%2.2dx%2.2d", i, j);
		sprintf(strl, "Difference from SiPM average time. Channel %d.%d;ns", i, j);
		hDt[i][j] = new TH1D(strs, strl, 400, -200, 200);
	}
	iNevtTotal = 0;
	upTime = 0;
	fileFirstTime = -1;
	fileLastTime = -1;
}

//------------------------------->

/***
 *
 * A function which is called once per each event
 *
 * Input parameters: -
 *
 * Return value: -
 *
 ***/

void ReadDigiDataUser::processUserEvent()
{
	float eSiPM;
	float ePMT;
	float eVeto;
	float tSiPM;
	int nSiPM;
	int nPMT;
	int nVeto;
	int i;

  	if( ttype() != 1 ) return;

  	iNevtTotal++;
	fileLastTime = globalTime();
	if (fileFirstTime < 0) fileFirstTime = fileLastTime;

	eSiPM = ePMT = eVeto = tSiPM = 0;
	nSiPM = nPMT = nVeto = 0;

  	for(i = 0; i < nhits(); i++ ) {

		switch(type(i)) {
		case SiPmHit:
			if (e(i) < EHITMIN) break;
			eSiPM += e(i);
			tSiPM += t(i);
			nSiPM++;
			break;
		case PmtHit:
			ePMT += e(i);
			nPMT++;
			break;
		case VetoHit:
			eVeto += signal(i) / 500;	// no calibration yet
			nVeto++;
			break;
		}
  	} // for i

	if (nSiPM >= NSIPMMIN && ePMT >= TOTALMIN) {
		tSiPM /= nSiPM;
  		for(i = 0; i < nhits(); i++ ) {
			switch(type(i)) {
			case SiPmHit:
			case PmtHit:
				if (e(i) < EHITMIN) break;
				hDt[adc(i)][adcChan(i)]->Fill(t(i) - tSiPM);
				break;
			case VetoHit:
				if (signal(i) / 500 < EHITMIN) break;	// no calibration yet
				hDt[adc(i)][adcChan(i)]->Fill(t(i) - tSiPM);
				break;
			}
		}
	}

  	return;
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
	int i, j;

  	TFile *fout = TFile::Open(chOutputFile, "recreate");
	for (i=0; i<iN_AdcBoards; i++) for (j=0; j<iNChannels_AdcBoard; j++) hDt[i][j]->Write();
  	fout->Close();
  
	printf("Total up time        %Ld seconds\n", upTime / GLOBALFREQ);
	printf("Total physics events %Ld\n", iNevtTotal);
	printf("Run completed in     %d seconds\n", time(NULL) - progStartTime);
  	return;
}

//--------------------->


void ReadDigiDataUser::userActionAtFileChange()
{
	upTime += fileLastTime - fileFirstTime;
	fileFirstTime = -1;
}

