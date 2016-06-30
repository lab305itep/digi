#include <stdio.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TH1D.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>

class HPainter;

void posi_spectr(char *base)
{
	TH1D *h[6][2];
	int i, j;
	char str[64];
	TCanvas *cv;
	TPaveStats *st;
	float y, dy;
	TText *txt;
	TCut cf[6];
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);
	gStyle->SetOptFit();

	HPainter *p = new HPainter(base);
	if (!p->IsOpen()) {
		printf("Something wrong with base=%s.\n", base);
		return;
	}
	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
	TCut cg06("gtDiff > 0.6");									// dead time
	// cut edges - more muon background there
	TCut cex("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cey("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cez("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut chot("!(PositronX[2] >= 80 && PositronX[2] < 81 && PositronX[1] >= 24 && PositronX[1] < 28)");
	// Distance cut - nothing beyond 100 cm in R and +-40cm in RZ
	TCut cr("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cr40("Distance < 40");									// strong cut
	TCut cg20("gtDiff > 2");
	TCut cg250("gtDiff < 25");									// strong cut
	TCut cxy("PositronX[0]>=0 && PositronX[1]>=0");							// strong cut
	TCut cn5("NeutronEnergy > 5 && NeutronHits >= 5");						// strong cut
	TCut cgamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
	TCut cgamma0("AnnihilationGammas > 0");								// strong cut

	cf[0] = cg06;
	cf[1] = cg20 && cr && cgamma;
	cf[2] = cf[1] && cex && cey && cez && chot;
	cf[3] = cf[2] && cgamma0;
	cf[4] = cf[3] && cr40 && cg250 && cn5;
	cf[5] = cf[4] && cxy;

	for (i=0; i<2; i++) {
		sprintf(str, "hEA%d", i);
		h[0][i] = new TH1D(str, "Positron spectrum, only time isolation cuts;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hEB%d", i);
		h[1][i] = new TH1D(str, "Positron spectrum, weak cuts;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hEC%d", i);
		h[2][i] = new TH1D(str, "Positron spectrum, weak cuts + fiducial volume;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hED%d", i);
		h[3][i] = new TH1D(str, "Positron spectrum, weak cuts + fiducial volume + gamma > 0;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hEE%d", i);
		h[4][i] = new TH1D(str, "Positron spectrum, strong cuts;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hEF%d", i);
		h[5][i] = new TH1D(str, "Positron spectrum, strong + xy cuts;E, MeV;mHz", 40, 0, 8);
	}
	
	printf("Histograms are created\n");
	
	for (i=0; i<6; i++) {
		p->Project(h[i][0], "PositronEnergy", cs && cf[i]);
		p->Project(h[i][1], "PositronEnergy", (!cs) && cf[i]);
	}
	
	printf("Projections are done\n");
	
	for (i=0; i<6; i++) {
		for (j=0; j<2; j++) h[i][j]->SetLineWidth(4);
		h[i][0]->SetLineColor(kGreen);
		h[i][1]->SetLineColor(kRed);
		h[i][1]->SetMinimum(0);
	}

	cv = new TCanvas(str, "Plots", 1200, 1200);
	txt = new TText;
	cv->Divide(3, 2);
	for (i=0; i<6; i++) {
		cv->cd(i+1);
		h[i][1]->Draw();
		h[i][0]->Draw("sames");
		gPad->Update();
		st = (TPaveStats *) h[i][0]->FindObject("stats");
		st->SetLineColor(kGreen);
		st->SetTextColor(kGreen);
		y = st->GetY1NDC();
		dy = st->GetY2NDC() - y;
		st->SetX1NDC(0.72);
		st->SetY1NDC(y + dy/2);
		st = (TPaveStats *) h[i][1]->FindObject("stats");
		st->SetLineColor(kRed);
		st->SetTextColor(kRed);
		st->SetX1NDC(0.72);
		st->SetY2NDC(y + dy/2);
		if (h[i][0]->GetMaximum() > h[i][1]->GetMaximum()) {
			h[i][0]->Draw();
			h[i][1]->Draw("sames");
		} else {
			h[i][1]->Draw();
			h[i][0]->Draw("sames");
		}
		txt->DrawTextNDC(0.7, 0.7, base);
		gPad->Update();
	}
	cv->Update();
}
