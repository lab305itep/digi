#include <stdio.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TTree.h>
#include "HPainter.h"

#define NRANDOM	16

HPainter::HPainter(const char *base)
{
	char strs[2048], strr[2048];

	snprintf(strs, sizeof(strs), "sel_%s.root", base);
	snprintf(strr, sizeof(strr), "sel_%s_random.root", base);

	Init(strs, strr);
}

HPainter::HPainter(const char *sname, const char *rname)
{
	Init(sname, rname);
}

void HPainter::Init(const char *sname, const char *rname)
{
	TTree *info;

	tSig = tRand = NULL;
	upTime = 0;
	
	fSig = new TFile(sname);
	if (!fSig->IsOpen()) {
		printf("No file %s\n", sname);
		return;
	}
	
	fRand = new TFile(rname);
	if (!fRand->IsOpen()) {
		printf("No file %s\n", rname);
		return;
	}
	
	tSig = (TTree *) fSig->Get("DanssPair");
	tRand = (TTree *) fRand->Get("DanssPair");
	
	info = (TTree *) fSig->Get("SumInfo");
	info->GetEntry(0);
	upTime = info->GetLeaf("gTime")->GetValue() / 125E6;
	gROOT->cd();
}

HPainter::~HPainter(void)
{
	fSig->Close();
	fRand->Close();
}

void HPainter::Project(TH1 *hist, const char *what, TCut cut)
{
	TH1 *hSig;
	TH1 *hRand;
	
	if (!IsOpen()) return;
	
	hSig = (TH1 *) hist->Clone("_sigtmp");
	hRand = (TH1 *) hist->Clone("_randtmp");
	
	tSig->Project("_sigtmp", what, cut.GetTitle());
	tRand->Project("_randtmp", what, cut.GetTitle());
	
	hSig->Sumw2();
	hRand->Sumw2();
	hRand->Scale(1.0/NRANDOM);
	
	hist->Add(hSig, hRand, 1.0, -1.0);
	hist->Scale(1000.0 / upTime);	// mHz
	
	delete hSig;
	delete hRand;
}

