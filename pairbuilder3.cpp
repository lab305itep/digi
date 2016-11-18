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

#include "evtbuilder.h"

#define GFREQ2US	(GLOBALFREQ / 1000000.0)
#define MAXTDIFF	50.0	// us
#define MINPOSE		1.0	// MeV
#define MAXPOSE		8.0	// MeV
#define AGAMMAN		0	// number of annihilation gamma hits (0 no requirement)
#define MINNEUTE	3.0	// MeV
#define MAXNEUTE	10.0	// MeV
#define NEUTN		3	// number of hits
#define MINVETOE	4.0	// MeV
#define VETON		2	// number of hits
#define DANSSVETOE	20.0	// Make veto if VETO counters are silent from Pmt or SiPM
#define VETOBLK		100.0	// us
#define ATTENUATION	0.00342	// Signal attenuation for positron energy correction

//	Correction based on the neutron position if this was not done before based on the positron position
void acorr(struct DanssPairStruct3 *DanssPair) {
	if (DanssPair->PositronX[0] < 0 && DanssPair->NeutronX[0] >= 0) {
		DanssPair->PositronEnergy *= exp(ATTENUATION * (DanssPair->NeutronX[0] - 50.0));
	} else if (DanssPair->PositronX[1] < 0 && DanssPair->NeutronX[1] >= 0) {
		DanssPair->PositronEnergy *= exp(ATTENUATION * (DanssPair->NeutronX[1] - 50.0));
	}
}

int IsNeutron(struct DanssEventStruct3 *DanssEvent)
{
	float E;
	int rc;
	
	E = (DanssEvent->SiPmCleanEnergy + DanssEvent->PmtCleanEnergy) / 2;
	rc = (E >= MINNEUTE && E < MAXNEUTE && DanssEvent->SiPmCleanHits >= NEUTN);

	return rc;
}

int IsPositron(struct DanssEventStruct3 *DanssEvent)
{
	float E;
	int rc;

	E = DanssEvent->PositronEnergy;
	rc = (E >= MINPOSE && E < MAXPOSE && DanssEvent->AnnihilationGammas >= AGAMMAN);

	return rc;
}

int IsVeto(struct DanssEventStruct3 *Event)
{
	if (Event->VetoCleanEnergy > MINVETOE || Event->VetoCleanHits >= VETON || Event->PmtCleanEnergy + Event->SiPmCleanEnergy > 2*DANSSVETOE) return 1;
	return 0;
}

void MakePair(struct DanssEventStruct3 *DanssEvent, struct DanssEventStruct3 *SavedEvent, struct DanssEventStruct3 *VetoEvent, struct DanssPairStruct3 *DanssPair)
{
	memset(DanssPair, 0, sizeof(struct DanssPairStruct3));
	DanssPair->number[0] = SavedEvent->number;
	DanssPair->number[1] = DanssEvent->number;
	DanssPair->unixTime = DanssEvent->unixTime;
	DanssPair->SiPmCleanEnergy[0] = SavedEvent->SiPmCleanEnergy;
	DanssPair->PmtCleanEnergy[0] = SavedEvent->PmtCleanEnergy;
	DanssPair->SiPmCleanEnergy[1] = DanssEvent->SiPmCleanEnergy;
	DanssPair->PmtCleanEnergy[1] = DanssEvent->PmtCleanEnergy;
	
	DanssPair->PositronHits = SavedEvent->PositronHits;
	DanssPair->PositronEnergy = SavedEvent->PositronEnergy;
	memcpy(DanssPair->PositronX, SavedEvent->PositronX, sizeof(SavedEvent->PositronX));
	DanssPair->MaxHitEnergy = SavedEvent->MaxHitEnergy;
	DanssPair->AnnihilationGammas = SavedEvent->AnnihilationGammas;
	DanssPair->AnnihilationEnergy = SavedEvent->AnnihilationEnergy;
	
	DanssPair->NeutronHits = DanssEvent->SiPmCleanHits;
	DanssPair->NeutronEnergy = (DanssEvent->SiPmCleanEnergy + DanssEvent->PmtCleanEnergy) / 2;
	memcpy(DanssPair->NeutronX, DanssEvent->NeutronX, sizeof(DanssEvent->NeutronX));
	DanssPair->NeutronRadius = DanssEvent->NeutronRadius;
	
	DanssPair->gtDiff = (DanssEvent->globalTime - SavedEvent->globalTime) / GFREQ2US;
	DanssPair->Distance = sqrt(
		(DanssEvent->NeutronX[0] - SavedEvent->PositronX[0]) * (DanssEvent->NeutronX[0] - SavedEvent->PositronX[0]) +
		(DanssEvent->NeutronX[1] - SavedEvent->PositronX[1]) * (DanssEvent->NeutronX[1] - SavedEvent->PositronX[1]) +
		(DanssEvent->NeutronX[2] - SavedEvent->PositronX[2]) * (DanssEvent->NeutronX[2] - SavedEvent->PositronX[2])
	);
	DanssPair->DistanceZ = DanssEvent->NeutronX[2] - SavedEvent->PositronX[2];

	DanssPair->gtFromVeto = (SavedEvent->globalTime - VetoEvent->globalTime) / GFREQ2US;
	DanssPair->VetoHits = VetoEvent->VetoCleanHits;
	DanssPair->VetoEnergy = VetoEvent->VetoCleanEnergy;
	DanssPair->DanssEnergy = (VetoEvent->SiPmCleanEnergy + VetoEvent->PmtCleanEnergy) / 2;
	
	acorr(DanssPair);
}

int main(int argc, char **argv)
{
	struct DanssEventStruct3		DanssEvent;
	struct DanssPairStruct3			DanssPair;
	struct DanssEventStruct3		Neutron;
	struct DanssEventStruct3		Positron;
	struct DanssEventStruct3		Veto;
	struct DanssInfoStruct3			DanssInfo;
	struct DanssInfoStruct			SumInfo;

	TChain *EventChain;
	TChain *InfoChain;
	TTree *tOut;
	TTree *InfoOut;
	TFile *fOut;
	FILE *fList;
	char str[1024];
	long long iEvt, nEvt;
	int PairReady;
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
		"PositronEnergy/F:"	// Energy sum of the cluster (SiPM)
		"MaxHitEnergy/F:"	// Energy of the maximum hit (SiPM)
		"PositronX[3]/F:"	// cluster position
		"AnnihilationGammas/I:"	// number of possible annihilation gammas
		"AnnihilationEnergy/F:"	// Energy in annihilation gammas
//		"neutron" parameters
		"NeutronHits/I:"	// number of hits considered as neutron capture gammas
		"NeutronEnergy/F:"	// Energy sum of above (SiPM)
		"NeutronX[3]/F:"	// center of gammas position
		"NeutronRadius/F:"	// average distance between hits and the center
//		Pair parameters
		"gtDiff/F:"		// time difference in us (from 125 MHz clock)
		"Distance/F:"		// distance between neutron and positron, cm
		"DistanceZ/F:"		// in Z, cm
//		Environment
		"gtFromPrevious/F:"	// time from the previous hit before positron, us
		"PreviousEnergy/F:"	// energy of the previous event
		"gtToNext/F:"		// time to the next hit after neutron, counted from positron, us
		"NextEnergy/F:"		// energy of the next event
		"EventsBetween/I:"	// Events between positron and neutron
//		Veto
		"gtFromVeto/F:"		// time from the last Veto event
		"VetoHits/I:"		// hits in Veto counters
		"VetoEnergy/F:"		// Energy in Veto counters
		"DanssEnergy/F"		// Veto Energy in Danss (Pmt + SiPm)/2
	);
	
	InfoOut = new TTree("SumInfo", "Summary information");
	InfoOut->Branch("Info", &SumInfo,  
		"gTime/L:"		// running time in terms of 125 MHz
		"startTime/I:"		// linux start time, seconds
		"stopTime/I:"		// linux stop time, seconds
		"events/L"		// number of events
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
	for (iEvt =0; iEvt < nEvt; iEvt++) {
		EventChain->GetEntry(iEvt);
//	Veto
		if (IsVeto(&DanssEvent)) {
			memcpy(&Veto, &DanssEvent, sizeof(struct DanssEventStruct3));
			continue;
		}
//	Get Neutron
		if (IsNeutron(&DanssEvent)) {
			memcpy(&Neutron, &DanssEvent, sizeof(struct DanssEventStruct3));
//	Now look backward for positron
			for (i=iEvt-1; i>=0; i--) {
				EventChain->GetEntry(i);
				if (Neutron.globalTime - DanssEvent.globalTime >= MAXTDIFF * GFREQ2US) break;
				if (Neutron.globalTime - DanssEvent.globalTime < 0) break;
				if (IsPositron(&DanssEvent)) break;
			}
			if (Neutron.globalTime - DanssEvent.globalTime < 0) continue;	// 125 MHz counter was restarted
//	less than 50 us from neutron
			if (Neutron.globalTime - DanssEvent.globalTime < MAXTDIFF * GFREQ2US && i >= 0) {
				memcpy(&Positron, &DanssEvent, sizeof(struct DanssEventStruct3));
				MakePair(&Neutron, &Positron, &Veto, &DanssPair);
				DanssPair.EventsBetween = iEvt - i - 1;
//	look backward
				if (i>0) {
					EventChain->GetEntry(i-1);
					DanssPair.gtFromPrevious = (Positron.globalTime - DanssEvent.globalTime) / GFREQ2US;
					DanssPair.PreviousEnergy = (DanssEvent.SiPmCleanEnergy + DanssEvent.PmtCleanEnergy) / 2;
				}
//	look forward
				if (iEvt + 1 < nEvt) {
					EventChain->GetEntry(iEvt+1);
					DanssPair.gtToNext = (DanssEvent.globalTime - Positron.globalTime) / GFREQ2US;
					DanssPair.NextEnergy = (DanssEvent.SiPmCleanEnergy + DanssEvent.PmtCleanEnergy) / 2;
				}
				tOut->Fill();
				PairCnt++;
			}
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

