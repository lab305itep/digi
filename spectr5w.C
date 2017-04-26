#include "HPainter.h"

TFile *fRoot;

TH1D *spectr5(const char *prefix, int mask, int run_from, int run_to, double bgnd, TCut cAux)
{
	char str[256];
	
	gStyle->SetOptStat(1001100);
//		Set cuts
	TCut cVeto("gtFromVeto > 60");
//	TCut cVeto("gtFromVeto > 200");
	TCut cIso("(gtFromPrevious > 45 || gtFromPrevious == gtFromVeto) && gtToNext > 80 && EventsBetween == 0");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe("PositronEnergy > 1");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
//        TCut cN("NeutronEnergy > 3");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cN && cPe && cIso && cAux;
	TCut cSig = cSel && cVeto;
	TCut cBgnd = cSel && (!cVeto);
//		Background tail correction
	TF1 fBgndN("fBgndN", "0.0163-0.0007*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.0722-0.0034*x", 0, 100);

	sprintf(str, "%s_hSig", prefix);
	TH1D *hSig  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	sprintf(str, "%s_hCosm", prefix);
	TH1D *hBgnd = new TH1D(str, "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	sprintf(str, "%s_hRes", prefix);
	TH1D *hRes  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	
	HPainter *ptr = new HPainter(mask, run_from, run_to);
	ptr->SetFile(fRoot);
	ptr->Project(hSig,  "PositronEnergy", cSig);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd);

	hSig->Add(&fBgndN, -1);
	hBgnd->Add(&fBgndC, -1);
	hRes->Add(hSig, hBgnd, 1, -bgnd);
	fRoot->cd();
	hSig->Write();
	hBgnd->Write();
	hRes->Write();
	delete ptr;
	return hRes;
}

void spectr_all(const char *fname = "danss_report.root", TCut cAux = (TCut)"", double bgScale = 2.24)	// 5.6% from reactor OFF data
{
#include "positions.h"
	const int mask = 0x801E;
	int i;
	int N;
	
	fRoot = new TFile (fname, "RECREATE");
	N = sizeof(positions) / sizeof(positions[0]);
	for (i=0; i<N; i++) 
		spectr5(positions[i].name, mask, positions[i].first, positions[i].last, positions[i].bgnd * bgScale, cAux);
	fRoot->Close();
}

