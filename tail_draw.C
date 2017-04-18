#include "HPainter.h"

void tail_draw(int run_first = 5808, int run_last = 15028)
{
	char strs[128];
	char strl[1024];
	TH1D *h[3][2];
	int i, j;
	
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut cXY("PositronX[0] >= 0 && PositronX[1] >= 0");
	TCut cT20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe9("PositronEnergy > 9");
        TCut cR60("Distance < 60");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct;

	TCanvas *cv;
	TVirtualPad *pd;
	TChain *chain;
	TH2D *hXY[2];
	TH1D *hX[2];
	TH1D *hY[2];
	TH1D *hZ[2];
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

	HPainter *hp = new HPainter(0x801E, run_first, run_last);
	chain = hp->GetPairChain();
	
	for (j=0; j<2; j++) {
		hXY[j] = new TH2D((j) ? "hXYC" : "hXYN", (j) ? "Cosmic;X, cm;Y, cm" : "Neutrino;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
		hX[j]  = new TH1D((j) ? "hXC" : "hXN", "Positron X distribution;X, cm", 25, 0, 100);
		hX[j]->GetYaxis()->SetLabelSize(0.06);
		hX[j]->SetLineColor((j) ? kRed : kBlue);
		hY[j]  = new TH1D((j) ? "hYC" : "hYN", "Positron Y distribution;Y, cm", 25, 0, 100);
		hY[j]->GetYaxis()->SetLabelSize(0.06);
		hY[j]->SetLineColor((j) ? kRed : kBlue);
		hZ[j]  = new TH1D((j) ? "hZC" : "hZN", "Positron Z distribution;Z, cm", 100, 0, 100);
		hZ[j]->GetYaxis()->SetLabelSize(0.06);
		hZ[j]->SetLineColor((j) ? kRed : kBlue);
	}

	for (j=0; j<2; j++) {
		ct = cIso && cXY && cZ && cT20 && cR && cPe9 && cGamma && cN;
		chain->Project(hXY[j]->GetName(), "PositronX[1]+2:PositronX[0]+2", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cY && cZ && cT20 && cR && cPe9 && cGamma && cN && "PositronX[0] >= 0";
		chain->Project(hX[j]->GetName(), "PositronX[0]+2", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cZ && cT20 && cR && cPe9 && cGamma && cN && "PositronX[1] >= 0";
		chain->Project(hY[j]->GetName(), "PositronX[1]+2", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cT20 && cR && cPe9 && cGamma && cN;
		chain->Project(hZ[j]->GetName(), "PositronX[2]+0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
	}
	
	cv = new TCanvas("CVXY", "XY", 1200, 700);
	cv->Divide(2, 1);
	for (j=0; j<2; j++) {
		pd = cv->cd(j+1);
		pd->SetRightMargin(0.12);
		hXY[j]->Draw("colz");
	}

	cv = new TCanvas("CVXYZ", "XYZ", 1200, 800);
	cv->Divide(3, 1);
	cv->cd(1);
	hX[1]->SetMinimum(0);
	hX[1]->Draw("");
	hX[0]->Draw("same");
	cv->cd(2);
	hY[1]->SetMinimum(0);
	hY[1]->Draw("");
	hY[0]->Draw("same");
	cv->cd(3);
	hZ[1]->SetMinimum(0);
	hZ[1]->Draw("");
	hZ[0]->Draw("same");
	TLegend *lg = new TLegend(0.2, 0.7, 0.7, 0.8);
	lg->SetTextSize(0.04);
	lg->AddEntry(hZ[0], "Neutrino", "l");
	lg->AddEntry(hZ[1], "Cosmic", "l");
	lg->Draw();
}
