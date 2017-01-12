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
	tBegin = 0;
	tEnd = 0;
	
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
	TCut ct;
	char str[256];
	
	if (!IsOpen()) return;
	
	hSig = (TH1 *) hist->Clone("_sigtmp");
	hRand = (TH1 *) hist->Clone("_randtmp");
	if (tBegin < tEnd) {
		sprintf(str, "unixTime > %d && unixTime < %d", tBegin, tEnd);
		ct = cut && str;
		tSig->Project("_sigtmp", what, ct);
		tRand->Project("_randtmp", what, ct);
	} else {
		tSig->Project("_sigtmp", what, cut.GetTitle());
		tRand->Project("_randtmp", what, cut.GetTitle());
	}
	
	hSig->Sumw2();
	hRand->Sumw2();
	hRand->Scale(1.0/NRANDOM);
	
	hist->Add(hSig, hRand, 1.0, -1.0);
	hist->Scale(1000.0 / upTime);	// mHz
	
	delete hSig;
	delete hRand;
}

void HPainter::SetUpTime(unsigned int t0, unsigned int t1)
{
	int i, N, S;
	TH1S *h;
	tBegin = t0;
	tEnd = t1;
	N = (t1 - t0) / 100;
	h = (TH1S*) gROOT->FindObject("_tmp_utime");
	if (h) delete h;
	h = new TH1S("_tmp_utime", "unix time", N, tBegin, tBegin + 100*N);
	tSig->Project("_tmp_utime", "unixTime");
	S = 0;
	for (i=0; i<N; i++) if (h->GetBinContent(i+1)) S++;
	upTime = S*100;
}

