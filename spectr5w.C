#include "HPainter2.h"

TFile *fRoot;

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
        TCut cRZ("fabs(DistanceZ) < 40");
//        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cR = cR2 && (cRXY || cR1);
        TCut cN("NeutronEnergy > 3.5");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cGammaMax && cN && cPe && cIso && cShower && cAux;
	cSig = cSel && cVeto;
	cBgnd = cSel && (!cVeto);
}

TH1D *spectr5(int nSect, const char *prefix, int mask, int run_from, int run_to, double bgnd, TCut cAux)
{
	char str[256];
	TCut cSig;
	TCut cBgnd;
	
	gStyle->SetOptStat(1001100);
	make_cuts(cSig, cBgnd, cAux);
//		Background tail correction
//	TF1 fBgndN("fBgndN", "0.01370-0.00057*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.06217-0.00288*x", 0, 100);
	TF1 fBgndN("fBgndN", "0.01370-0.00057*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.06217-0.00288*x", 0, 100);

	sprintf(str, "%s_hSig", prefix);
	TH1D *hSig  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	sprintf(str, "%s_hCosm", prefix);
	TH1D *hBgnd = new TH1D(str, "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	sprintf(str, "%s_hRes", prefix);
	TH1D *hRes  = new TH1D(str,  "Positron spectrum;MeV;mHz/0.25 MeV", 60, 1, 16);
	
	HPainter2 *ptr = new HPainter2(mask, run_from, run_to);
	ptr->SetFile(fRoot);
	ptr->Project(hSig,  "PositronEnergy", cSig);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd);

	hSig->Add(&fBgndN, -1.0/nSect);
	hBgnd->Add(&fBgndC, -1.0/nSect);
	hRes->Add(hSig, hBgnd, 1, -bgnd);
	fRoot->cd();
	hSig->Write();
	hBgnd->Write();
	hRes->Write();
	delete ptr;
	return hRes;
}

void spectr_all(int nSect, const char *fname = "danss_report_v4.root", TCut cAux = (TCut)"", double bgScale = 2.24)	// 5.6% from reactor OFF data
{
#include "positions.h"
	const int mask = 0x801E;
	int i;
	int N;
	
	fRoot = new TFile (fname, "RECREATE");
	N = sizeof(positions) / sizeof(positions[0]);
	for (i=0; i<N; i++) {
		spectr5(nSect, positions[i].name, mask, positions[i].first, positions[i].last, positions[i].bgnd * bgScale, cAux);
		printf(".");
		fflush(stdout);
	}
	printf("\n");
	fRoot->Close();
}

void make_sect_cut(TCut &cSect, int num, int nSect)
{
	double s3_cuts[3][2] = {{3.5, 34.5}, {34.5, 64.5}, {64.5, 95.5}};	// special case for old calculations compatibility
	double zmin, zmax;
	char str[1024];
	
	if (nSect < 1 || num < 0 || num >= nSect) {
		cSect = (TCut) "0";	// forbidden
		return;
	} else if (nSect == 3) {	// special case
		zmin = s3_cuts[num][0];
		zmax = s3_cuts[num][1];
	} else {
		zmin = 3.5 + 92.0 * num / nSect;
		zmax = 3.5 + 92.0 * (num + 1) / nSect;
	}
	sprintf(str, "PositronX[2] > %7.2f && PositronX[2] <= %7.2f", zmin, zmax);
	cSect = (TCut) str;
}

void spectr_sect(int nSect = 3, const char *fname = "danss_report_v4.root", TCut cAux = (TCut)"", double bgScale = 2.24)	// 5.6% from reactor OFF data
{
	char str[1024];
	char *ptr;
	int i;
	TCut cSect;
	
	for (i=0; i<nSect; i++) {
		strcpy(str, fname);
		ptr = strstr(str, ".root");
		if (!ptr) ptr = &str[strlen(str)];
		sprintf(ptr, "-sect%d_of_%d.root", i+1, nSect);
		make_sect_cut(cSect, i, nSect);
		spectr_all(nSect, str, cAux && cSect, bgScale);
	}
}

TH1D *count_z5(int mask, int run_from, int run_to, double bgnd, TCut cAux)
{
	char str[256];
	TCut cSig;
	TCut cBgnd;
	
	make_cuts(cSig, cBgnd, cAux);
	TH1D *hSig  = new TH1D("hSigz5", ";Z, cm", 1000, 0, 100);
	TH1D *hBgnd = new TH1D("hBgndz5", ";Z, cm", 1000, 0, 100);
	TH1D *hRes  = new TH1D("hResz5", ";Z, cm", 1000, 0, 100);
	
	HPainter2 *ptr = new HPainter2(mask, run_from, run_to);
	ptr->Project(hSig,  "PositronX[2]+0.5", cSig);
	ptr->Project(hBgnd, "PositronX[2]+0.5", cBgnd);

	hRes->Add(hSig, hBgnd, 1, -bgnd);
	delete ptr;
	printf("MeanZ=%7.3f +- %5.3f\n", hRes->GetMean(), hRes->GetMeanError());
	return hRes;
}

void z5_sect(int nSect, int mask, int run_from, int run_to, TCut cAux = (TCut)"", double bgnd = 0.056)
{
	int i;
	TCut cSect;
	for (i=0; i<nSect; i++) {
		printf("Sect %d: ", i + 1);
		make_sect_cut(cSect, i, nSect);
		count_z5(mask, run_from, run_to, bgnd, cSect && cAux);
	}
}
