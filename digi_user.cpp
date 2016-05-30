/***
 *
 * Version:                             2.0
 *
 * Package:                             DANSS SiPm Signal Processing and Calibration
 *
 * Description:                         Example user functions to process digitized SiPm data
 *
 ***/

#include <stdio.h>

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

#define SIPMEMIN	0.15
#define VETOEMIN	6.0
#define GLOBALFREQ	125000000

using namespace std;

// Globals:

long long                               iNevtTotal;
long long                               iNevtNoVeto;
long long				upTime;
long long				fileFirstTime;
long long				fileLastTime;
int                                     progStartTime;

// Here one can list histograms:

TH1D *                                  hNSiPM;
TH1D *                                  hNPMT;
TH1D *                                  hESiPM;
TH1D *                                  hEPMT;

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
	for( int i = 0; i < iNElements; i++ ) {
    		int iAdcNum = static_cast<int> ( static_cast<float>(i) / 100.) ;
    		int iAdcChan = i - 100 * iAdcNum;
    		if(!isAdcChannelExist(iAdcNum, iAdcChan)) continue;
    		setTd(i, 1); // set all td = 1
  	}
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
	progStartTime = time(NULL);
	hNSiPM = new TH1D("hNSiPM", "Number of SiPM hits", 50, 0, 50);
	hNPMT  = new TH1D("hNPMT" , "Number of PMT hits" , 10, 0, 10);
	hESiPM = new TH1D("hESiPM", "SiPM energy", 100, 0, 10);
	hEPMT  = new TH1D("hEPMT" , "PMT energy" , 100, 0, 10);
	iNevtTotal = iNevtNoVeto = 0;
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

	eSiPM = ePMT = eVeto = 0;
	nSiPM = nPMT = nVeto = 0;

  	for(i = 0; i < nhits(); i++ ) {

		switch(type(i)) {
		case SiPmHit:
			if (e(i) < SIPMEMIN) break;
			eSiPM += e(i);
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

	if (eVeto > VETOEMIN || nVeto > 1) return;
	iNevtNoVeto++;    

  	hNSiPM->Fill(nSiPM);
  	hESiPM->Fill(eSiPM);
  	hNPMT->Fill(nPMT);
  	hEPMT->Fill(ePMT);

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
  	TFile *fout = TFile::Open("digi_anal.root", "recreate");
  	hNSiPM->Write();
  	hESiPM->Write();
  	hNPMT->Write();
  	hEPMT->Write();
  	fout->Close();
  
	printf("Total up time        %Ld seconds\n", upTime / GLOBALFREQ);
	printf("Total physics events %Ld\n", iNevtTotal);
	printf("No veto events       %Ld\n", iNevtNoVeto);
	printf("Run completed in     %d seconds\n", time(NULL) - progStartTime);
  	return;
}

//--------------------->


void ReadDigiDataUser::userActionAtFileChange()
{
	upTime += fileLastTime - fileFirstTime;
	fileFirstTime = -1;
}

