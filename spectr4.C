#include "HPainter.h"

TH1D *spectr4(const char *name, const char *namer, TCut cP = "")
{
//	gROOT->ProcessLine(".L HPainter.cpp+");
//	gROOT->ProcessLine(".x cuts-nov16.C");
	gStyle->SetOptStat(1001100);

	TH1D *hSig  = new TH1D("hSig",  "Positron Energy;MeV;mHz", 30, 2, 8);
	TH1D *hBgnd = new TH1D("hBgnd", "Positron Energy;MeV;mHz", 30, 2, 8);
	TH1D *hRes  = new TH1D("hRes",  "Positron Energy;MeV;mHz", 30, 2, 8);
	
	HPainter *ptr = new HPainter(name, namer);
	
	ptr->Project(hSig,  "PositronEnergy", cSig && cP);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd && cP);

	hRes->Add(hSig, hBgnd, 1, -0.05);
	hRes->Draw();
	return hRes;
}

