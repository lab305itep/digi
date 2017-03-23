#include "HPainter.h"
void random_bgnd_Dist(void)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);

	TFile *fRoot = new TFile("random_bgnd_Dist.root");
	if (!fRoot->IsOpen()) {
		printf("Call random_bgnd_Dist_calc() first\n");
		return;
	}
	
	TH1D * hSig1 = (TH1D *) fRoot->Get("hR1-sig");
	TH1D * hRand1 = (TH1D *) fRoot->Get("hR1-rand");
	TH1D * hDiff1 = (TH1D *) fRoot->Get("hR1-diff");
	TH1D * hSig2 = (TH1D *) fRoot->Get("hR2-sig");
	TH1D * hRand2 = (TH1D *) fRoot->Get("hR2-rand");
	TH1D * hDiff2 = (TH1D *) fRoot->Get("hR2-diff");
	if (!hSig1 || !hRand1 || !hDiff1 || !hSig2 || !hRand2 || !hDiff2) {
		printf("Call random_bgnd_gtDiff_calc() first\n");
		return;
	}

	hDiff1->SetMarkerStyle(kFullCircle);
	hDiff1->SetMarkerSize(2);
	hDiff1->SetLineColor(kGreen);
	hDiff1->SetMarkerColor(kGreen);

	hSig1->SetMarkerStyle(22);
	hSig1->SetMarkerSize(2);
	hSig1->SetLineColor(kYellow);
	hSig1->SetMarkerColor(kYellow);
	hSig1->GetYaxis()->SetLabelSize(0.05);
	
	hRand1->SetMarkerStyle(kOpenCircle);
	hRand1->SetMarkerSize(2);
	hRand1->SetLineColor(kBlue);
	hRand1->SetMarkerColor(kBlue);

	hDiff2->SetMarkerStyle(kFullSquare);
	hDiff2->SetMarkerSize(2);
	hDiff2->SetLineColor(kRed);
	hDiff2->SetMarkerColor(kRed);

	hSig2->SetMarkerStyle(27);
	hSig2->SetMarkerSize(2);
	hSig2->SetLineColor(kOrange);
	hSig2->SetMarkerColor(kOrange);
	hSig2->GetYaxis()->SetLabelSize(0.05);
	
	hRand2->SetMarkerStyle(kOpenSquare);
	hRand2->SetMarkerSize(2);
	hRand2->SetLineColor(kViolet);
	hRand2->SetMarkerColor(kViolet);

	TCanvas *c1 = new TCanvas("CV", "Random bgnd: Dist", 1000, 1000);
	TLegend *lg = new TLegend(0.45, 0.75, 0.9, 0.9);
	lg->AddEntry(hDiff1, "Neutrino events, X or Y", "P");
	lg->AddEntry(hRand1, "Random background, X or Y", "P");
	lg->AddEntry(hDiff2, "Neutrino events, X and Y", "P");
	lg->AddEntry(hRand2, "Random background, X and Y", "P");
	hRand1->SetTitle(";us;");
	hRand1->SetMinimum(0);
	hRand1->DrawCopy();
	hDiff2->DrawCopy("same");
	hRand2->DrawCopy("same");
	hDiff1->DrawCopy("same");

	lg->Draw();
	c1->Update();
	fRoot->Close();
}

void random_bgnd_Dist_calc(void)
{
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR2("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cSel1 = cX && cY && cZ && c20 && cGamma && cPe && cVeto && cIso && !cR2;
	TCut cSel2 = cX && cY && cZ && c20 && cGamma && cPe && cVeto && cIso && cR2;

	TH1D *hDist1 = new TH1D("hR1", ";cm;", 40, 0, 160);
	TH1D *hDist2 = new TH1D("hR2", ";cm;", 40, 0, 160);
	HPainter *hp = new HPainter(0x801E, 5808, 11688);
	TFile *fRoot = new TFile("random_bgnd_Dist.root", "RECREATE");
	hp->SetFile(fRoot);
	hp->Project(hDist1, "Distance", cSel1);
	hp->Project(hDist2, "Distance", cSel2);

	fRoot->cd();
	hDist1->Write();
	hDist2->Write();
	delete hp;
	fRoot->Close();
}
