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

#define GFREQ2MS	(GLOBALFREQ / 1000)
#define MAXT		100	// ms
#define MUONENERGY	800	// MeV
#define POSITRONMIN	1	// MeV
#define POSITRONMAX	12	// MeV
#define DAY		(24.0*60*60)


void Add2Chain(TChain *ch, int run_from, int run_to, int mask, const char *f_template)
{
	const char stat_file_name[] = "stat_all.txt";
	FILE *f;
	char str[1024];
	char fname[1024];
	char *ptr;
	int i, num, cond, N;

	f = fopen(stat_file_name, "rt");
	if (!f) {
		printf("Stat file %s not found: %m\n", stat_file_name);
		return;
	}

	N = 0;
	for(i=0;;i++) {
		ptr = fgets(str, sizeof(str), f);
		if (!ptr || feof(f)) break;
		ptr = strtok(str, " \t");
		if (!isdigit(ptr[0])) continue;
		num = strtol(ptr, NULL, 10);
		ptr = strtok(NULL, " \t");
		if (!ptr) {
			printf("Strange record at line %d\n", i);
			continue;
		}
		cond = strtol(ptr, NULL, 10);
		if (num < run_from) continue;
		if (num > run_to) break;
		if (cond <= 0) continue;
		if (mask & (1 << (cond - 1))) {
			sprintf(fname, f_template, num);
			ch->AddFile(fname);
			N++;
		}
	}
	fclose(f);
	printf("%d runs found.\n", N);
}

int main(int argc, char **argv)
{
	struct DanssEventStruct4	DanssEvent;
	struct DanssInfoStruct4		DanssInfo;

	TChain *EventChain;
	TChain *InfoChain;
	char str[1024];
	long long iEvt, nEvt;
	long long lastMuon;
	long long SelCnt, MuonCnt;
	long long upTime;
	double fUpTime;
	int i;
	int iFirst, iLast;
	double t, e;
	
	TH1D *hE = new TH1D("hE", "Cluster Energy;MeV", 4 * (POSITRONMAX - POSITRONMIN), POSITRONMIN, POSITRONMAX);
	TH1D *hT = new TH1D("hT", "Time difference;ms", MAXT, 0, MAXT);
	
	if (argc < 3) {
		printf("Usage: %s first_run last_run\n", argv[0]);
		return 10;
	}
	iFirst = strtol(argv[1], NULL, 10);
	iLast = strtol(argv[2], NULL, 10);

	EventChain = new TChain("DanssEvent");
	EventChain->SetBranchAddress("Data", &DanssEvent);
	Add2Chain(EventChain, iFirst, iLast, 0x801E, "/mnt/root0/danss_root4/danss_%6.6d.root");
	InfoChain = new TChain("DanssInfo");
	InfoChain->SetBranchAddress("Info", &DanssInfo);
	Add2Chain(InfoChain, iFirst, iLast, 0x801E, "/mnt/root0/danss_root4/danss_%6.6d.root");

	nEvt = EventChain->GetEntries();
	SelCnt = MuonCnt = 0;
	lastMuon = -GLOBALFREQ;
	printf("Starting for %Ld events.\n", nEvt);
	for (iEvt = 0; iEvt < nEvt; iEvt++) {
		EventChain->GetEntry(iEvt);
		if (!(iEvt % 1000000)) {
			printf(">");
			fflush(stdout);
		}
		if (lastMuon > DanssEvent.globalTime) {
			lastMuon = -GLOBALFREQ;		// time break			
		} else if (DanssEvent.SiPmCleanEnergy + DanssEvent.PmtCleanEnergy > 2 * MUONENERGY) {
			MuonCnt++;
			lastMuon = DanssEvent.globalTime;
		} else if (DanssEvent.globalTime - lastMuon < GFREQ2MS * MAXT && 
			DanssEvent.PositronEnergy > POSITRONMIN &&
			DanssEvent.PositronEnergy < POSITRONMAX) {
			SelCnt++;
			t = 1.0 * (DanssEvent.globalTime - lastMuon) / GFREQ2MS;
			e = DanssEvent.PositronEnergy;
			if (t > 5 && t < 50) hE->Fill(e);
			if (e > 6) hT->Fill(t);
		}
	}
	
	upTime = 0;
	for(i=0; i<InfoChain->GetEntries(); i++) {
		InfoChain->GetEntry(i);
		upTime += DanssInfo.upTime;
	}
	upTime /= GLOBALFREQ;
	fUpTime = upTime / DAY;

	hT->Scale(1/fUpTime);
	hE->Scale(1/fUpTime);

	printf("\n%Ld events (%5.0f days) processed - %Ld muons, %Ld events.\n", nEvt, fUpTime, MuonCnt, SelCnt);

	TCanvas *cv = new TCanvas("CV", "CV", 1200, 800);
	
	cv->Divide(2, 1);
	cv->cd(1)->SetLogy();
	hE->Draw();
	cv->cd(2);
	hT->Draw();
	cv->SaveAs("mu12B.pdf");

	delete cv;
	delete EventChain;
	delete InfoChain;
	delete hE;
	delete hT;

	return 0;
}

