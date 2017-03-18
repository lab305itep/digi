#include "HPainter.h"

TFile *fRoot;

TH1D *spectr5(const char *prefix, int mask, int run_from, int run_to, double bgnd, TCut cAux)
{
	char str[256];
//	TFile *fRoot;
	
	gStyle->SetOptStat(1001100);
//		Set cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR("Distance < 60 && DistanceZ > -40 && DistanceZ < 40");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cPe && cAux;
	TCut cSig = cSel && cVeto && cIso;
	TCut cBgnd = cSel && (!cVeto) && (cIso || "gtFromPrevious == gtFromVeto");
//#ifdef STRONG_CUTS
//        TCut cStrong("NeutronEnergy > 4 && gtDiff < 20");
//        cSig = cSig && cStrong;
//        cBgnd = cBgnd && cStrong;
//	fRoot = new TFile ("danss_report_strong.root", "UPDATE");
//#else
//	fRoot = new TFile ("danss_report.root", "UPDATE");
//#endif

	sprintf(str, "%s_hSig", prefix);
	TH1D *hSig  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.2 MeV", 55, 1, 12);
	sprintf(str, "%s_hCosm", prefix);
	TH1D *hBgnd = new TH1D(str, "Positron spectrum;MeV;mHz/0.2 MeV", 55, 1, 12);
	sprintf(str, "%s_hRes", prefix);
	TH1D *hRes  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.2 MeV", 55, 1, 12);
	
	HPainter *ptr = new HPainter(mask, run_from, run_to);
	ptr->SetFile(fRoot);
	ptr->Project(hSig,  "PositronEnergy", cSig);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd);

	hRes->Add(hSig, hBgnd, 1, -bgnd);
	fRoot->cd();
	hSig->Write();
	hBgnd->Write();
	hRes->Write();
//	fRoot->Close();
	delete ptr;
	return hRes;
}

void spectr_all(const char *fname = "danss_report.root", TCut cAux = (TCut)"")
{
#include "positions.h"
	const int mask = 0x801E;
	int i;
	int N;
	
	fRoot = new TFile (fname, "UPDATE");
	N = sizeof(positions) / sizeof(positions[0]);
	for (i=0; i<N; i++) 
		spectr5(positions[i].name, mask, positions[i].first, positions[i].last, positions[i].bgnd, cAux);
	fRoot->Close();
}
