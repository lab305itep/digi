#include <math.h>
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

#include "evtbuilder.h"

#define GFREQ2US	(GLOBALFREQ / 1000000)
#define MAXTDIFF	50.0	// us
#define MINPOSE		1.0	// MeV
#define MAXPOSE		8.0	// MeV
#define MINNEUTE	3.0	// MeV
#define MAXNEUTE	10.0	// MeV
#define NEUTN		3	// number of gammas
#define MINVETOE	4.0	// MeV
#define VETON		2	// number of hits
#define VETOBLK		100.0	// us

int IsNeutron(struct DanssEventStruct *DanssEvent)
{
	if (DanssEvent->NeutronSiPmEnergy < MINNEUTE || DanssEvent->NeutronSiPmEnergy > MAXNEUTE) return 0;
	if (DanssEvent->PmtCleanEnergy < MINNEUTE || DanssEvent->PmtCleanEnergy > MAXNEUTE) return 0;
	if (DanssEvent->SiPmCleanEnergy < MINNEUTE || DanssEvent->SiPmCleanEnergy > MAXNEUTE) return 0;
	
	if (DanssEvent->SiPmCleanHits < NEUTN) return 0;
	return 1;
}

int IsPositron(struct DanssEventStruct *DanssEvent)
{
	if (DanssEvent->PositronSiPmEnergy < MINPOSE || DanssEvent->PmtCleanEnergy < MINPOSE || DanssEvent->SiPmCleanEnergy < MINPOSE) return 0;
	if (DanssEvent->PositronSiPmEnergy > MAXPOSE) return 0;
	return 1;
}

int IsVeto(struct DanssEventStruct *Event)
{
	if (Event->VetoCleanEnergy > MINVETOE || Event->VetoCleanHits > VETON) return 1;
	return 0;
}

void MakePair(struct DanssEventStruct *DanssEvent, struct DanssEventStruct *SavedEvent, struct DanssPairStruct *DanssPair)
{
	DanssPair->number[0] = SavedEvent->number;
	DanssPair->number[1] = DanssEvent->number;
	DanssPair->unixTime = DanssEvent->unixTime;
	DanssPair->SiPmCleanEnergy[0] = SavedEvent->SiPmCleanEnergy;
	DanssPair->PmtCleanEnergy[0] = SavedEvent->PmtCleanEnergy;
	DanssPair->SiPmCleanEnergy[1] = DanssEvent->SiPmCleanEnergy;
	DanssPair->PmtCleanEnergy[1] = DanssEvent->PmtCleanEnergy;
	
	DanssPair->PositronHits = SavedEvent->PositronHits;
	DanssPair->PositronSiPmEnergy = SavedEvent->PositronSiPmEnergy;
	memcpy(DanssPair->PositronX, SavedEvent->PositronX, sizeof(SavedEvent->PositronX));
	DanssPair->MaxHitEnergy = SavedEvent->MaxHitEnergy;
	DanssPair->AnnihilationGammas = SavedEvent->AnnihilationGammas;
	DanssPair->AnnihilationEnergy = SavedEvent->AnnihilationEnergy;
	
	DanssPair->NeutronHits = DanssEvent->NeutronHits;
	DanssPair->NeutronSiPmEnergy = DanssEvent->NeutronSiPmEnergy;
	memcpy(DanssPair->NeutronX, DanssEvent->NeutronX, sizeof(DanssEvent->NeutronX));
	memcpy(DanssPair->NeutronGammaEnergy, DanssEvent->NeutronGammaEnergy, sizeof(DanssEvent->NeutronGammaEnergy));
	memcpy(DanssPair->NeutronGammaDistance, DanssEvent->NeutronGammaDistance, sizeof(DanssEvent->NeutronGammaDistance));
	DanssPair->NeutronRadius = DanssEvent->NeutronRadius;
	
	DanssPair->gtDiff = DanssEvent->globalTime - SavedEvent->globalTime;
	DanssPair->Distance = sqrt(
		(DanssEvent->NeutronX[0] - SavedEvent->PositronX[0]) * (DanssEvent->NeutronX[0] - SavedEvent->PositronX[0]) +
		(DanssEvent->NeutronX[1] - SavedEvent->PositronX[1]) * (DanssEvent->NeutronX[1] - SavedEvent->PositronX[1]) +
		(DanssEvent->NeutronX[2] - SavedEvent->PositronX[2]) * (DanssEvent->NeutronX[2] - SavedEvent->PositronX[2])
	);
	DanssPair->DistanceZ = DanssEvent->NeutronX[2] - SavedEvent->PositronX[2];
}

int main(int argc, char **argv)
{
	struct DanssEventStruct			DanssEvent;
	struct DanssPairStruct			DanssPair;
	struct DanssEventStruct			SavedEvent;
	struct DanssInfoStruct			DanssInfo;
	struct DanssInfoStruct			SumInfo;

	TChain *EventChain;
	TChain *InfoChain;
	TTree *tOut;
	TTree *InfoOut;
	TFile *fOut;
	FILE *fList;
	char str[1024];
	long long iEvt, nEvt;
	long long lastgTime, lastVeto;
	int PairCnt;
	int i;
	char *ptr;
	
	if (argc < 3) {
		printf("Usage: ./pairbuilder list_file.txt output_file.root\n");
		printf("Will process files in the list_file and create root-file\n");
		return 10;
	}

	fOut = new TFile(argv[2], "RECREATE");
	if (!fOut->IsOpen()) {
		printf("Can not open the output file %s: %m\n", argv[2]);
		return -10;
	}

	tOut = new TTree("DanssPair", "Time Correlated events");
	tOut->Branch("Pair", &DanssPair,
		"number[2]/L:"		// event numbers in the file
		"unixTime/I:"		// linux time, seconds
		"SiPmCleanEnergy[2]/F:"	// Full Clean energy SiPm
		"PmtCleanEnergy[2]/F:"	// Full Clean energy Pmt
//		"positron cluster" parameters
		"PositronHits/I:"	// hits in the cluster
		"PositronSiPmEnergy/F:"	// Energy sum of the cluster (SiPM)
		"MaxHitEnergy/F:"	// Energy of the maximum hit (SiPM)
		"PositronX[3]/F:"	// cluster position
		"AnnihilationGammas/I:"	// number of possible annihilation gammas
		"AnnihilationEnergy/F:"	// Energy in annihilation gammas
//		"neutron" parameters
		"NeutronHits/I:"	// number of hits considered as neutron capture gammas
		"NeutronSiPmEnergy/F:"	// Energy sum of above (SiPM)
		"NeutronX[3]/F:"	// center of gammas position
		"NeutronGammaEnergy[5]/F:"	// sorted list of the 5 most energetic gammas
		"NeutronGammaDistance[5]/F:"	// distances for the gammas above to the "neutron" center
		"NeutronRadius/F:"	// average distance between hits and the center
//		Pair parameters
		"gtDiff/F:"		// time difference in us (from 125 MHz clock)
		"Distance/F:"		// distance between neutron and positron, cm
		"DistanceZ/F"		// in Z, cm
	);
	
	InfoOut = new TTree("SumInfo", "Summary information");
	InfoOut->Branch("Info", &SumInfo,  
		"gTime/L:"		// running time in terms of 125 MHz
		"startTime/I:"		// linux start time, seconds
		"stopTime/I:"		// linux stop time, seconds
		"events/I"		// number of events
	);
	memset(&SumInfo, 0, sizeof(struct DanssInfoStruct));

	EventChain = new TChain("DanssEvent");
	EventChain->SetBranchAddress("Data", &DanssEvent);
	InfoChain = new TChain("DanssInfo");
	InfoChain->SetBranchAddress("Info", &DanssInfo);

	fList = fopen(argv[1], "rt");
	if (!fList) {
		printf("Can not open list of files %s: %m\n", argv[1]);
		goto fin;
	}
	
	for(;;) {
		if (!fgets(str, sizeof(str), fList)) break;
		ptr = strchr(str, '\n');
		if (ptr) *ptr = '\0';
		EventChain->Add(str);
		InfoChain->Add(str);
	}
	fclose(fList);

	nEvt = EventChain->GetEntries();
	PairCnt = 0;
	lastVeto = lastgTime = -GLOBALFREQ;
	for (iEvt =0; iEvt < nEvt; iEvt++) {
		EventChain->GetEntry(iEvt);
		if (IsVeto(&DanssEvent)) lastVeto = DanssEvent.globalTime;
		if (DanssEvent.globalTime - lastVeto < VETOBLK * GFREQ2US) {
			continue;	// Veto is active
		}
		if (DanssEvent.globalTime - lastgTime < MAXTDIFF * GFREQ2US && IsNeutron(&DanssEvent)) {
			MakePair(&DanssEvent, &SavedEvent, &DanssPair);
			tOut->Fill();
			lastgTime = -GLOBALFREQ;
			PairCnt++;
			continue;
		} else if (IsPositron(&DanssEvent)) {
			memcpy(&SavedEvent, &DanssEvent, sizeof(struct DanssEventStruct));
			lastgTime = DanssEvent.globalTime;
		}
	}
	
	for(i=0; i<InfoChain->GetEntries(); i++) {
		InfoChain->GetEntry(i);
		SumInfo.upTime += DanssInfo.upTime;
		SumInfo.stopTime = DanssInfo.stopTime;
		SumInfo.events += DanssInfo.events;
		if (!i) SumInfo.startTime = DanssInfo.startTime;
	}
	InfoOut->Fill();	

	printf("%Ld events processed - %d pairs found. Aquired time %f7.0 s\n", iEvt, PairCnt, SumInfo.upTime / GLOBALFREQ);
fin:
	delete EventChain;
	delete InfoChain;

	InfoOut->Write();
	tOut->Write();
	fOut->Close();
	return 0;
}

