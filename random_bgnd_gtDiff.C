#include "HPainter.h"
void random_bgnd_gtDiff(void)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
//	gStyle->SetLineWidth(4);

	TFile *fRoot = new TFile("random_bgnd_gtDiff.root");
	if (!fRoot->IsOpen()) {
		printf("Call random_bgnd_gtDiff_calc() first\n");
		return;
	}
	
	TH1D * hSig = (TH1D *) fRoot->Get("hGt-sig");
	TH1D * hRand = (TH1D *) fRoot->Get("hGt-rand");
	TH1D * hDiff = (TH1D *) fRoot->Get("hGt-diff");
	if (!hSig || !hRand || !hDiff) {
		printf("Call random_bgnd_gtDiff_calc() first\n");
		return;
	}

	TF1 *fDec = new TF1("FDEC", "[0]*(exp(-x/[1]) - exp(-x/[2]))", 0.0, 50.0);
	fDec->SetParNames("Const.", "#tau_{capture}", "#tau_{therm.}");
	fDec->SetLineWidth(4);
	fDec->SetLineColor(kOrange);

	TF1 *fPol = new TF1("FPOL", "pol1", 0.0, 50.0);
	fPol->SetLineWidth(4);
	fPol->SetLineColor(kOrange);

	hDiff->SetMarkerStyle(21);
	hDiff->SetMarkerSize(2);
	hDiff->SetLineColor(kGreen);
	hDiff->SetMarkerColor(kGreen);

	hSig->SetMarkerStyle(22);
	hSig->SetMarkerSize(2);
	hSig->SetLineColor(kRed);
	hSig->SetMarkerColor(kRed);
	hSig->GetYaxis()->SetLabelSize(0.05);
	
	hRand->SetMarkerStyle(23);
	hRand->SetMarkerSize(2);
	hRand->SetLineColor(kBlue);
	hRand->SetMarkerColor(kBlue);

	fDec->SetParameters(hDiff->Integral(), 15, 5);
//	hDiff->Fit(fDec, "0", "", 1, 50);
//	hRand->Fit(fPol, "0", "", 1, 50);


	TCanvas *c1 = new TCanvas("CV", "Random bgnd: gtDiff", 1000, 1000);
	TLegend *lg = new TLegend(0.45, 0.75, 0.9, 0.9);
	lg->AddEntry(hSig, "Events found", "LP");
	lg->AddEntry(hRand, "Random background", "LP");
//	lg->AddEntry(hDiff, "Neutrino events", "LP");
	hSig->SetTitle(";us;");
	hSig->SetMinimum(0);
	hSig->DrawCopy();
//	hDiff->DrawCopy("sames");
//	fDec->Draw("same");
	hRand->DrawCopy("sames");
//	fPol->Draw("same");

	lg->Draw();
	c1->Update();
	fRoot->Close();
}

void random_bgnd_gtDiff_calc(void)
{
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR("Distance < 40 && DistanceZ > -40 && DistanceZ < 40");
	TCut c10("gtDiff > 1");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cSel = cX && cY && cZ && cR && c10 && cGamma && cPe;
	TCut cSig = cSel && cVeto && cIso;

	TH1D *hGt = new TH1D("hGt", ";us;", 50, 0, 50);
	HPainter *hp = new HPainter(0x801E, 5808, 11688);
	TFile *fRoot = new TFile("random_bgnd_gtDiff.root", "RECREATE");
	hp->SetFile(fRoot);
	hp->Project(hGt, "gtDiff", cSig);

	TH1D *hR1 = new TH1D("hR1", ";cm;", 150, 0, 150);
	TH1D *hR2 = new TH1D("hR2", ";cm;", 150, 0, 150);
	TCut cR0 = cX && cY && cZ && c20 && cGamma && cPe && cVeto && cIso;
	TCut cXY = ("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	hp->Project(hR1, "Distance", cR0 && !cXY);
	hp->Project(hR2, "Distance", cR0 && cXY);
	
	fRoot->cd();
	hGt->Write();
	hR1->Write();
	hR2->Write();
	delete hp;
	fRoot->Close();
}
