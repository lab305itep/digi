#include <stdio.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>

#include "HPainter.h"

void spectr(char *base1, float bgnd1 = 0.05, char *base2 = NULL, float bgnd2 = 0.05)
{
	TH1D *h[2][2];	// spectrum
	TH1D *hSum[2];
	int i, j;
	char str[64], strl[1024];
	HPainter *p[2];
	TCanvas *cv;

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);

	p[0] = new HPainter(base1);
	if (!p[0]->IsOpen()) {
		printf("Something wrong with base=%s.\n", base1);
		return;
	}
	if (base2 && strlen(base2)) {
		p[1] = new HPainter(base2);
		if (!p[1]->IsOpen()) {
			printf("Something wrong with base=%s.\n", base2);
			return;
		}
	} else {
		p[1] = NULL;
	}

	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
	// cut edges - more muon background there
	TCut cex("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cey("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cez("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut chot("!(PositronX[2] >= 80 && PositronX[2] < 81 && PositronX[1] >= 24 && PositronX[1] < 28)");
	// Distance cut - nothing beyond 100 cm in R and +-40cm in RZ
	TCut cr("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cg20("gtDiff > 2");
	TCut cgamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
	TCut cpe("PositronEnergy>=1");

	cf = cg20 && cr && cgamma && cex && cey && cez && chot && cpe;

	for (i=0; i<2; i++) {
		sprintf(str, "hS%s%d", base1, i);
		sprintf(strl, "Positron spectrum, DANSS position = %s;E, MeV;mHz", base1);
		h[0][i] = new TH1D(str, strl, 40, 0, 8);
		if (p[1]) {
			sprintf(str, "hS%s%d", base2, i);
			sprintf(strl, "Positron spectrum, DANSS position = %s;E, MeV;mHz", base2);
			h[1][i] = new TH1D(str, strl, 40, 0, 8);			
		}
	}

	cv = new TCanvas("CV", "Spectrum", 1500, 600);
	cv->Divide(3, 1);
	
	for (i=0; i<2; i++) if (p[i]) {
		cv->cd(i+1);
		p[i]->Project(h[i][0], "PositronEnergy", cs && cf);
		p[i]->Project(h[i][1], "PositronEnergy", (!cs) && cf);
		h[i][1]->Scale((i) ? bgnd2 : bgnd1);
		h[i][0]->Add(h[i][1], -1);	// background
		h[i][0]->SetLineColor(kGreen);
		h[i][1]->SetLineColor(kRed);
		h[i][0]->DrawCopy();
		h[i][1]->DrawCopy("sames");
	}
	
	sprintf(str, "hs%s", base1);
	hSum[0] = (TH1D *) h[0][0]->Clone(str);
	sprintf(str, "hb%s", base1);
	hSum[1] = (TH1D *) h[0][1]->Clone(str);
	sprintf(strl, "Positron spectrum, DANSS position = %s+%s;E, MeV;mHz", base1, base2);
	if (p[1]) for (i=0; i<2; i++) {
		sprintf(strl, "Positron spectrum, DANSS position = %s+%s;E, MeV;mHz", base1, base2);
		hSum[i]->SetTitle(strl);
		hSum[i]->SetBit(TH1::kIsAverage);
		h[1][i]->SetBit(TH1::kIsAverage);
		hSum[i]->Add(h[1][i]);
	}
	
	cv->cd(3);
	hSum[0]->DrawCopy();
	hSum[1]->DrawCopy("sames");
	cv->Update();
	for (i=1; i<=40; i++) printf("%3d %8.5f +- %8.5f    %8.5f +- %8.5f\n", i, hSum[0]->GetBinContent(i), hSum[0]->GetBinError(i), hSum[1]->GetBinContent(i), hSum[1]->GetBinError(i));
}

