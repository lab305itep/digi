#include "HPainter.h"
void LongTimeScale(float ShowerThreshold = 800, int run_first = 2306, int run_last = 15028)
{
	char str[128];
	double upTime;
	TChain *ch;
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct;

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	TH1D *h = new TH1D("hLongTime", ";s;Events / day", 50, 0, 1);
	HPainter *hp = new HPainter(0x801E, run_first, run_last, "/space/danss_root2b/");
	ct = cIso && cX && cY && cZ && c20 && cR && cPe && cGamma && cN && cVeto;
	sprintf(str, "ShowerEnergy > %8.1f", ShowerThreshold);
	hp->Project(h, "gtFromShower/1E6", ct && str);
	upTime = hp->GetUpTime();
	h->Scale(86.4);
	h->SetLineWidth(2);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
//		800 MeV: 1319 events in 25859 s
//	TF1 *f2Exp = new TF1("f2Exp", "[0]*exp(-x/19.6)+[1]*exp(-x/0.2572)", 0, 10);
//		2500 MeV 105 events in 25859 s
	TF1 *f2Exp = new TF1("f2Exp", "[0]*exp(-x/246)+[1]*exp(-x/0.2572)", 0, 10);
	f2Exp->SetParameters(4, 0);
	f2Exp->SetParNames("Const.", "^{9}Li");
	new TCanvas("CV", "LongScale", 1200, 800);
	h->SetTitle("Fit with sum of two exponents");
	h->Fit(f2Exp, "", "", 0.02, 1);
}
