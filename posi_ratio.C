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

void posi_ratio(void)
{
	const char *base[] = {"up", "stuck", "down"};
	TH1D *hs[3][3];	// spectrum
	TH1D *hr[3][2];	// ratios
	int i, j;
	char str[64], strl[1024];
	TCanvas *cv;
	TPaveStats *st;
	float y, dy;
	TText *txt;
	TCut cf;
	HPainter *p[3];
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);
	gStyle->SetOptFit();

	for (i=0; i<3; i++) {
		p[i] = new HPainter(base[i]);
		if (!p[i]->IsOpen()) {
			printf("Something wrong with base=%s.\n", base[i]);
			return;
		}
	}

	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
//	TCut cg06("gtDiff > 0.6");									// dead time
	// cut edges - more muon background there
	TCut cex("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cey("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cez("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut chot("!(PositronX[2] >= 80 && PositronX[2] < 81 && PositronX[1] >= 24 && PositronX[1] < 28)");
	// Distance cut - nothing beyond 100 cm in R and +-40cm in RZ
	TCut cr("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
//	TCut cr40("Distance < 40");									// strong cut
	TCut cg20("gtDiff > 2");
//	TCut cg250("gtDiff < 25");									// strong cut
//	TCut cxy("PositronX[0]>=0 && PositronX[1]>=0");							// strong cut
//	TCut cn5("NeutronEnergy > 5 && NeutronHits >= 5");						// strong cut
	TCut cgamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
//	TCut cgamma0("AnnihilationGammas > 0");								// strong cut

	cf = cg20 && cr && cgamma && cex && cey && cez && chot;

	for (i=0; i<3; i++) for (j=0; j<3; j++) {
		sprintf(str, "hS%c%d", 'A'+i, j);
		sprintf(strl, "Positron spectrum, DANSS position = %s;E, MeV;mHz", base[i]);
		hs[i][j] = new TH1D(str, strl, 40, 0, 8);
	}
	for (i=0; i<3; i++) for (j=0; j<2; j++) {
		sprintf(str, "hR%c%d", 'A'+i, j);
		hr[i][j] = new TH1D(str, "Positron spectrum ratio;E, MeV;mHz", 40, 0, 8);
	}
	
	printf("Histograms are created\n");
	
	for (i=0; i<3; i++) {
		p[i]->Project(hs[i][0], "PositronEnergy", cs && cf);
		p[i]->Project(hs[i][1], "PositronEnergy", (!cs) && cf);
		hs[i][2]->Add(hs[i][0], hs[i][1], 1.0, -0.05);	// 5% background
	}
	
	printf("Projections are done\n");
	for (i=0; i<3; i++) {
		j = (i+2)%3;
		hr[i][0]->Divide(hs[i][0], hs[j][0]);
		hr[i][1]->Divide(hs[i][2], hs[j][2]);
	}	
	printf("Ratios are done\n");
	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			hs[i][j]->SetLineWidth(4);
			hs[i][j]->SetMinimum(0);
		}
		hs[i][0]->SetLineColor(kGreen);
		hs[i][1]->SetLineColor(kRed);
	}

	hs[0][2]->SetLineColor(kRed);
	hs[1][2]->SetLineColor(kGreen);
	hs[2][2]->SetLineColor(kBlue);
	for (i=0; i<2; i++) {
		hr[0][i]->SetLineColor(kRed);
		hr[1][i]->SetLineColor(kGreen);
		hr[2][i]->SetLineColor(kBlue);
	}

	for (i=0; i<3; i++) for (j=0; j<2; j++) {
		hr[i][j]->SetMinimum(0.4);
		hr[i][j]->SetMaximum(1.8);
	}

	cv = new TCanvas(str, "Plots", 1200, 1200);
	cv->Divide(3, 2);
	txt = new TText;
//		Spectrum
	for (i=0; i<3; i++) {
		cv->cd(i+1);
		hs[i][1]->Draw();
		hs[i][0]->Draw("sames");
		gPad->Update();
		st = (TPaveStats *) hs[i][0]->FindObject("stats");
		st->SetLineColor(kGreen);
		st->SetTextColor(kGreen);
		y = st->GetY1NDC();
		dy = st->GetY2NDC() - y;
		st->SetX1NDC(0.72);
		st->SetY1NDC(y + dy/2);
		st = (TPaveStats *) hs[i][1]->FindObject("stats");
		st->SetLineColor(kRed);
		st->SetTextColor(kRed);
		st->SetX1NDC(0.72);
		st->SetY2NDC(y + dy/2);
		if (hs[i][0]->GetMaximum() > hs[i][1]->GetMaximum()) {
			hs[i][0]->Draw();
			hs[i][1]->Draw("sames");
		} else {
			hs[i][1]->Draw();
			hs[i][0]->Draw("sames");
		}
		txt->DrawTextNDC(0.7, 0.7, base[i]);
		gPad->Update();
	}
//		Spectrum, background subtracted
	cv->cd(4);
	hs[0][2]->SetTitle("Positron spectrum");
	hs[0][2]->Draw();
	hs[1][2]->Draw("sames");
	hs[2][2]->Draw("sames");
	gPad->Update();
	st = (TPaveStats *) hs[0][2]->FindObject("stats");
	st->SetLineColor(kRed);
	st->SetTextColor(kRed);
	y = st->GetY1NDC();
	dy = st->GetY2NDC() - y;
	st->SetX1NDC(0.72);
	st->SetY1NDC(y + dy/2);
	st = (TPaveStats *) hs[1][2]->FindObject("stats");
	st->SetLineColor(kGreen);
	st->SetTextColor(kGreen);
	st->SetX1NDC(0.72);
	st->SetY2NDC(y + dy/2);
	st = (TPaveStats *) hs[2][2]->FindObject("stats");
	st->SetLineColor(kBlue);
	st->SetTextColor(kBlue);
	st->SetX1NDC(0.72);
	st->SetY2NDC(y);
	st->SetY1NDC(y - dy/2);
	hs[0][2]->Draw();
	hs[1][2]->Draw("sames");
	hs[2][2]->Draw("sames");
	TLegend *lg = new TLegend(0.7, 0.5, 0.9, 0.65);
	for (i=0; i<3; i++) lg->AddEntry(hs[i][2], base[i], "l");
	lg->Draw();
//		Ratios
	TLegend *lgr = new TLegend(0.7, 0.1, 0.9, 0.25);
	for (i=0; i<3; i++) {
		sprintf(str, "%s/%s", base[i], base[(i+2)%3]);
		lgr->AddEntry(hs[i][2], str, "l");
	}
	for (i=0; i<2; i++) {
		cv->cd(i+5);
		hr[0][i]->Fit("pol0");
		hr[1][i]->Fit("pol0", "", "sames");
		hr[2][i]->Fit("pol0", "", "sames");
		gPad->Update();
		st = (TPaveStats *) hr[0][i]->FindObject("stats");
		st->SetOptStat(0);
		st->SetOptFit(1);
		st->SetLineColor(kRed);
		st->SetTextColor(kRed);
		y = st->GetY1NDC();
		dy = st->GetY2NDC() - y;
		st->SetY1NDC(y + dy/2);
		st->SetX1NDC(0.72);
		st = (TPaveStats *) hr[1][i]->FindObject("stats");
		st->SetOptStat(0);
		st->SetOptFit(1);
		st->SetLineColor(kGreen);
		st->SetTextColor(kGreen);
		st->SetX1NDC(0.72);
		st->SetY2NDC(y + dy/2);
		st = (TPaveStats *) hr[2][i]->FindObject("stats");
		st->SetOptStat(0);
		st->SetOptFit(1);
		st->SetLineColor(kBlue);
		st->SetTextColor(kBlue);
		st->SetX1NDC(0.72);
		st->SetY2NDC(y);
		st->SetY1NDC(y - dy/2);
		lgr->Draw();
	}
//
	cv->Update();
}

