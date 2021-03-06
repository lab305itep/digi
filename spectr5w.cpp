#include <stdio.h>
#include <stdlib.h>

#include <TCut.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>

#include "HPainter2.h"


void make_cuts(TCut &cSig, TCut &cBgnd, TCut cAux)
{
//		Set cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cMuonA("gtFromVeto == 0");
	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
//        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1);
        TCut cN("NeutronEnergy > 3.5");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cGammaMax && cN && cPe && cIso && cShower && cAux;
	cSig = cSel && cVeto;
	cBgnd = cSel && (!cVeto);
}

void spectr5(const char *name, int mask, int run_from, int run_to, double bgnd)
{
	char str[256];
	TCut cSig;
	TCut cBgnd;
	TCut cAux;
	int nSect;
	double bgScale;
	char *ptr;
	char fname[1024];
	char pair_dir[] = "/mnt/root1/danss_pair6";
	
//	Environment
	nSect = 1;
	ptr = getenv("SPECTR_NSECT");
	if (ptr) nSect = strtol(ptr, NULL, 10);
	
	bgScale = 2.24;
	ptr = getenv("SPECTR_BGSCALE");
	if (ptr) bgScale = strtod(ptr, NULL);
	
	cAux = (TCut)"";
	ptr = getenv("SPECTR_AUXCUT");
	if (ptr) cAux = (TCut)ptr;
	
	sprintf(fname, "period/%s.root", name);
	ptr = getenv("SPECTR_OUTDIR");
	if (ptr) sprintf(fname, "%s/%s.root", ptr, name);
	TFile *fRoot = new TFile (fname, "RECREATE");
	
	ptr = getenv("SPECTR_PAIRDIR");
	if (!ptr) ptr = pair_dir;
	
	make_cuts(cSig, cBgnd, cAux);
//		Background tail correction (mHz)
//	TF1 fBgndN("fBgndN", "0.01370-0.00057*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.06217-0.00288*x", 0, 100);
	TF1 fBgndN("fBgndN", "0.01426-0.000613*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07142-0.003486*x", 0, 100);

	TH1D *hSig  = new TH1D("hSig",  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	TH1D *hBgnd = new TH1D("hCosm", "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	TH1D *hRes  = new TH1D("hRes",  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	TH1D *hConst = new TH1D("hConst", "Various period parameters", 10, 0, 10);
	
	HPainter2 *ptr2 = new HPainter2(mask, run_from, run_to, ptr);
	ptr2->SetFile(fRoot);
	ptr2->Project(hSig,  "PositronEnergy", cSig);
	ptr2->Project(hBgnd, "PositronEnergy", cBgnd);
	hConst->Fill("UpTime", ptr2->GetUpTime());

	hSig->Add(&fBgndN, -1.0/nSect);
	hBgnd->Add(&fBgndC, -1.0/nSect);
	hRes->Add(hSig, hBgnd, 1, -bgnd * bgScale);
	fRoot->cd();
	hSig->Write();
	hBgnd->Write();
	hRes->Write();
	hConst->Write();
	fRoot->Close();
	delete ptr2;
}

int main(int argc, char **argv)
// void spectr_all(int nSect, const char *fname = "danss_report_v4.root", TCut cAux = (TCut)"", double bgScale = 2.24)	// 5.6% from reactor OFF data
{
#include "positions.h"
	const int mask = 0x801E;
	int i;
	int N;
	

	if (argc < 2) {
		printf("Usage: %s period_number (from 1)\n", argv[0]);
		return 100;
	}
	i = strtol(argv[1], NULL, 10);
	N = sizeof(positions) / sizeof(positions[0]);
	if (i<1 || i > N) {
		printf("Illegal period number %d (max = %d)\n", i, N);
		return 110;
	}
	
	spectr5(positions[i-1].name, mask, positions[i-1].first, positions[i-1].last, positions[i-1].bgnd);
	return 0;
}
