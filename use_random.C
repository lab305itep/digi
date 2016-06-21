#include <TCanvas.h>
#include <TColor.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include "evtbuilder.h"

TH1D *hrat[6][4];

struct RTabStruct {
	double min;
	double max;
	int n;
	double *r;
};

struct RTabStruct RTab[6];

double ksingle(double x, struct RTabStruct *Tab)
{
	int k;
	double val;
	
	if (x < Tab->min || x > Tab->max) return 1E-10;		// something small
	k = Tab->n * (x - Tab->min) / (Tab->max - Tab->min);
	val = Tab->r[k+1];
	if (val <= 0) val = 1E-10;
	return val;
}

double kfunction(struct DanssPairStruct2 *pair) 
{
	double x[6];
	int i;
	double sum;
	
	x[0] = pair->Distance;
	x[1] = pair->gtDiff;
	x[2] = pair->NeutronHits;
	x[3] = pair->NeutronEnergy;
	x[4] = pair->AnnihilationGammas;
	x[5] = pair->AnnihilationEnergy;
	
	sum = 0;
	for(i=0; i<6; i++) sum += TMath::Log(ksingle(x[i], &RTab[i]));
	
	return sum;
}

void use_random(char *fname, char *rname)
{
	int i, j, N;
	char str[64];
	struct DanssPairStruct2 pair;
	float ksum;
//		constants
	const char cs[] = "EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100 && gtDiff > 0.6 && gtFromVeto > 100";
	const enum EColor colors[] = {kBlue, kRed, kGreen, kBlack};
//		Root settings
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);
//		Histogramms
	for (i=0; i<4; i++) {
		sprintf(str, "hd%d", i);
		hrat[0][i] = new TH1D(str, "Positron to Neutron distance;cm", 40, 0, 160);
		sprintf(str, "ht%d", i);
		hrat[1][i] = new TH1D(str, "Positron to Neutron time;us", 50, 0, 50);
		sprintf(str, "hnn%d", i);
		hrat[2][i] = new TH1D(str, "Hits in Neutron capture", 20, 0, 20);
		sprintf(str, "hne%d", i);
		hrat[3][i] = new TH1D(str, "Energy in Neutron capture", 50, 0, 10);
		sprintf(str, "hpn%d", i);
		hrat[4][i] = new TH1D(str, "Hits beyond positron cluster", 10, 0, 10);
		sprintf(str, "hpe%d", i);
		hrat[5][i] = new TH1D(str, "Energy beyond positron cluster", 20, 0, 2);

		for (j=0; j<6; j++) hrat[j][i]->SetLineColor(colors[i]);
	}
//		Files and trees
	TFile *fs = new TFile(fname, "UPDATE");
	if (!fs->IsOpen()) {
		printf("Can not open file %s\n", fname);
		return;
	}
	TFile *fr = new TFile(rname, "UPDATE");
	if (!fr->IsOpen()) {
		printf("Can not open file %s\n", rname);
		return;
	}
	TTree *ts = (TTree *) fs->Get("DanssPair");
	if (!ts) {
		printf("No tree DanssPair in file %s\n", fname);
		return;
	}
	TTree *tr = (TTree *) fr->Get("DanssPair");
	if (!tr) {
		printf("No tree DanssPair in file %s\n", rname);
		return;
	}
	gROOT->cd();
//		Projections
	ts->Project("hd0", "Distance", cs);
	tr->Project("hd1", "Distance", cs);
	ts->Project("ht0", "gtDiff", cs);
	tr->Project("ht1", "gtDiff", cs);
	ts->Project("hnn0", "NeutronHits", cs);
	tr->Project("hnn1", "NeutronHits", cs);
	ts->Project("hne0", "NeutronEnergy", cs);
	tr->Project("hne1", "NeutronEnergy", cs);
	ts->Project("hpn0", "AnnihilationGammas", cs);
	tr->Project("hpn1", "AnnihilationGammas", cs);
	ts->Project("hpe0", "AnnihilationEnergy", cs);
	tr->Project("hpe1", "AnnihilationEnergy", cs);
//		Subtract and divide
	for (j=0; j<6; j++) {
		hrat[j][0]->Sumw2();
		hrat[j][1]->Sumw2();
		hrat[j][2]->Add(hrat[j][0], hrat[j][1], 1, -1);
		hrat[j][3]->Divide(hrat[j][2], hrat[j][0]);
	}
//		Draw
	TCanvas *cv = new TCanvas("CV", "Random background", 1200, 800);
	cv->Divide(3, 2);
	for (j=0; j<6; j++) {
		cv->cd(j+1);
		hrat[j][0]->SetMinimum(0);
		hrat[j][0]->Draw("hist");
		hrat[j][1]->Draw("hist,same");
		hrat[j][2]->Draw("same");
	}
//		Crerate tables
	for (j=0; j<6; j++) {
		RTab[j].r = hrat[j][3]->GetArray();
		RTab[j].n = hrat[j][3]->GetXaxis()->GetNbins();
		RTab[j].min = hrat[j][3]->GetXaxis()->GetXmin();
		RTab[j].max = hrat[j][3]->GetXaxis()->GetXmax();
		for(i=0; i<RTab[j].n; i++) if (RTab[j].r[i+1] < 0) RTab[j].r[i+1] = 0;
	}
//		KSUM - signal
	ts->SetBranchAddress("Pair", &pair);
	TTree *tks = (TTree *) fs->Get("KTree");
	if (tks) {
		tks->Reset();
		tks->SetBranchAddress("K", &ksum);
	} else {
		fs->cd();
		tks = new TTree("KTree", "K-values");
		tks->Branch("K", &ksum, "K/F");
	}
	N = ts->GetEntries();
	for(i=0; i<N; i++) {
		ts->GetEntry(i);
		ksum = kfunction(&pair);
		tks->Fill();
	}
	fs->cd();
	tks->Write();
//		KSUM - random
	tr->SetBranchAddress("Pair", &pair);
	TTree *tkr = (TTree *) fr->Get("KTree");
	if (tkr) {
		tkr->Reset();
		tkr->SetBranchAddress("K", &ksum);
	} else {
		fr->cd();
		tkr = new TTree("KTree", "K-values");
		tkr->Branch("K", &ksum, "K/F");
	}
	N = tr->GetEntries();
	for(i=0; i<N; i++) {
		tr->GetEntry(i);
		ksum = kfunction(&pair);
		tkr->Fill();
	}
	fr->cd();
	tkr->Write();	
//		Close files
	fs->Close();
	fr->Close();
//
	printf("Everything is done.\n");
}

void k_draw(char *fname, char *rname)
{
	TH1D *hk[3];
	const char cs[] = "EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100 && gtDiff > 0.6 && gtFromVeto > 100";
//		Hists
	hk[0] = new TH1D("hk0", "K-distribution", 100, -20, 0);
	hk[1] = new TH1D("hk1", "K-distribution", 100, -20, 0);
	hk[2] = new TH1D("hk2", "K-distribution", 100, -20, 0);
	hk[0]->SetLineColor(kBlue);
	hk[1]->SetLineColor(kRed);
	hk[2]->SetLineColor(kGreen);
//		Files and trees
	TFile *fs = new TFile(fname);
	if (!fs->IsOpen()) {
		printf("Can not open file %s\n", fname);
		return;
	}
	TFile *fr = new TFile(rname);
	if (!fr->IsOpen()) {
		printf("Can not open file %s\n", rname);
		return;
	}
	TTree *ts = (TTree *) fs->Get("DanssPair");
	if (!ts) {
		printf("No tree DanssPair in file %s\n", fname);
		return;
	}
	TTree *tr = (TTree *) fr->Get("DanssPair");
	if (!tr) {
		printf("No tree DanssPair in file %s\n", rname);
		return;
	}
	TTree *tks = (TTree *) fs->Get("KTree");
	if (!tks) {
		printf("No tree KTree in file %s\n", fname);
		return;
	}
	TTree *tkr = (TTree *) fr->Get("KTree");
	if (!tkr) {
		printf("No tree KTree in file %s\n", rname);
		return;
	}
	ts->AddFriend(tks);
	tr->AddFriend(tkr);
//		Project	
	gROOT->cd();
	ts->Project("hk0", "K", cs);
	tr->Project("hk1", "K", cs);

	hk[0]->Sumw2();	
	hk[1]->Sumw2();	
	hk[2]->Add(hk[0], hk[1], 1, -1);
	hk[0]->Draw("hist");
	hk[1]->Draw("hist,same");
	hk[2]->Draw("same");
}

void posi_draw(char *fname, char *rname, char *cut)
{
	TH1D *hp[3];
//		Hists
	hp[0] = new TH1D("hp0", "Positron kinetic energy, MeV", 40, 0, 8);
	hp[1] = new TH1D("hp1", "Positron kinetic energy, MeV", 40, 0, 8);
	hp[2] = new TH1D("hp2", "Positron kinetic energy, MeV", 40, 0, 8);
	hp[0]->SetLineColor(kBlue);
	hp[1]->SetLineColor(kRed);
	hp[2]->SetLineColor(kGreen);
	TCut cs("EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100 && gtDiff > 0.6 && gtFromVeto > 100");
//		Files and trees
	TFile *fs = new TFile(fname);
	if (!fs->IsOpen()) {
		printf("Can not open file %s\n", fname);
		return;
	}
	TFile *fr = new TFile(rname);
	if (!fr->IsOpen()) {
		printf("Can not open file %s\n", rname);
		return;
	}
	TTree *ts = (TTree *) fs->Get("DanssPair");
	if (!ts) {
		printf("No tree DanssPair in file %s\n", fname);
		return;
	}
	TTree *tr = (TTree *) fr->Get("DanssPair");
	if (!tr) {
		printf("No tree DanssPair in file %s\n", rname);
		return;
	}
	TTree *tks = (TTree *) fs->Get("KTree");
	if (!tks) {
		printf("No tree KTree in file %s\n", fname);
		return;
	}
	TTree *tkr = (TTree *) fr->Get("KTree");
	if (!tkr) {
		printf("No tree KTree in file %s\n", rname);
		return;
	}
	ts->AddFriend(tks);
	tr->AddFriend(tkr);
//		Project	
	gROOT->cd();
	ts->Project("hp0", "PositronEnergy", cs && cut);
	tr->Project("hp1", "PositronEnergy", cs && cut);

	hp[0]->Sumw2();
	hp[1]->Sumw2();	
	hp[2]->Add(hp[0], hp[1], 1, -1);
	
	TCanvas *cp = new TCanvas("CP", "Positron", 1200, 800);
	cp->Divide(2, 1);
	cp->cd(1);
	hp[0]->SetStats(0);
	hp[1]->SetStats(0);
	hp[2]->SetStats(0);
	hp[0]->DrawCopy("hist");
	hp[1]->DrawCopy("hist,same");
	hp[2]->DrawCopy("same");
	cp->cd(2);
	gStyle->SetOptStat(1000000);
	hp[2]->UseCurrentStyle();
	hp[2]->DrawCopy();	
}

