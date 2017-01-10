#include "HPainter.h"

TH1D *spectr3(const char *name, const char *namer, const char *resname, unsigned int t0 = 0, unsigned int t1 = 0, double bgnd = 0.05)
{
//	gROOT->ProcessLine(".L HPainter.cpp+");
//	gROOT->ProcessLine(".x cuts-nov16.C");
	gStyle->SetOptStat(1001100);

	TH1D *hSig  = new TH1D("hSig",  "Positron Energy;MeV;mHz", 35, 1, 8);
	TH1D *hBgnd = new TH1D("hBgnd", "Positron Energy;MeV;mHz", 35, 1, 8);
	TH1D *hRes  = new TH1D(resname,  "Positron Energy;MeV;mHz", 35, 1, 8);
	
	HPainter *ptr = new HPainter(name, namer);
	if (t1 > t0) ptr->SetUpTime(t0, t1);
	
	ptr->Project(hSig,  "PositronEnergy", cSig);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd);

	hRes->Add(hSig, hBgnd, 1, -bgnd);
	hRes->Draw();
	printf("Running time = %g s\n", ptr->GetUpTime());
	delete hSig;
	delete hBgnd;
	delete ptr;
	return hRes;
}

