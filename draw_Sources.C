#include <stdio.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TChain.h>
#include <TStyle.h>
#include <TText.h>
#include <TLine.h>
#include <TLegend.h>
#include <TCut.h>
#include <TGraph.h>
#include <TRandom2.h>
#include <TROOT.h>

TRandom2 rnd;

class MyRandom {
    public:
	inline MyRandom(void) {;};
	inline ~MyRandom(void) {;};
	static inline double Gaus(double mean = 0, double sigma = 1) 
		{
		return rnd.Gaus(mean, sigma);
	};
	static inline double GausAdd(double val, double sigma)
	{
		return rnd.Gaus(val, sqrt(val)*sigma);
	};
	static inline double GausAdd2(double val, double sigma)
	{
		return rnd.Gaus(val, val*sigma);
	};
};

void draw_Exp(TChain *tExpA, TChain *tExpB, TChain *tInfoA, TChain *tInfoB, 
	const char *name, const char *fname, TCut cXY, TCut cZ, double Efit[2], double kSP = 0.5)
{
	char str[256];
	long long gtA, gtB;
	double timeA, timeB;
	int i;
	
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
//	gStyle->SetOptStat(0);
//	gStyle->SetOptFit(0);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
//	gStyle->SetLineWidth(4);
	
	TH2D *hXY = new TH2D("hXY", "XY distribution of gamma flash center;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	sprintf(str, "DANSS energy deposit in %s decay;E, MeV", name);
	TH1D *hExpA = new TH1D("hExpA", str, 70, 0, 7);
	TH1D *hExpB = new TH1D("hExpB", str, 70, 0, 7);
	TH1D *hExpC = new TH1D("hExpC", str, 70, 0, 7);
	sprintf(str, "SiPM energy deposit in %s decay;E, MeV", name);
	TH1D *hExpSiPMA = new TH1D("hExpSiPMA", str, 70, 0, 7);
	TH1D *hExpSiPMB = new TH1D("hExpSiPMB", str, 70, 0, 7);
	TH1D *hExpSiPMC = new TH1D("hExpSiPMC", str, 70, 0, 7);
	sprintf(str, "PMT energy deposit in %s decay;E, MeV", name);
	TH1D *hExpPMTA = new TH1D("hExpPMTA", str, 70, 0, 7);
	TH1D *hExpPMTB = new TH1D("hExpPMTB", str, 70, 0, 7);
	TH1D *hExpPMTC = new TH1D("hExpPMTC", str, 70, 0, 7);
	sprintf(str, "Number of hits from %s decay", name);
	TH1D *hHitsA = new TH1D("hHitsA", str, 20, 0, 20);
	TH1D *hHitsB = new TH1D("hHitsB", str, 20, 0, 20);
	TH1D *hHitsC = new TH1D("hHitsC", str, 20, 0, 20);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
	TCut cn("SiPmCleanHits > 5");
	TCut cSel = cxyz && cVeto;
	
	sprintf(str, "%6.4f*SiPmCleanEnergy+%6.4f*PmtCleanEnergy", kSP, 1-kSP);
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cSel && cZ && cn);
	tExpA->Project("hExpA", str, cSel && cZ && cXY && cn);
	tExpB->Project("hExpB", str, cSel && cZ && cXY && cn);
	tExpA->Project("hExpSiPMA", "SiPmCleanEnergy", cSel && cZ && cXY && cn);
	tExpB->Project("hExpSiPMB", "SiPmCleanEnergy", cSel && cZ && cXY && cn);
	tExpA->Project("hExpPMTA", "PmtCleanEnergy", cSel && cZ && cXY && cn);
	tExpB->Project("hExpPMTB", "PmtCleanEnergy", cSel && cZ && cXY && cn);
	tExpA->Project("hHitsA", "SiPmCleanHits", cSel && cZ && cXY);
	tExpB->Project("hHitsB", "SiPmCleanHits", cSel && cZ && cXY);
	
	gtA = gtB = 0;
	for(i=0; i<tInfoA->GetEntries(); i++) {
		tInfoA->GetEntry(i);
		gtA += tInfoA->GetLeaf("gTime")->GetValueLong64();
	}
	for(i=0; i<tInfoB->GetEntries(); i++) {
		tInfoB->GetEntry(i);
		gtB += tInfoB->GetLeaf("gTime")->GetValueLong64();
	}
	timeA = gtA / 1.25E8;
	timeB = gtB / 1.25E8;
	
	printf("A = %6.1fs    B = %6.1fs\n", timeA, timeB);
	
	hExpA->Sumw2();
	hExpB->Sumw2();
	hExpSiPMA->Sumw2();
	hExpSiPMB->Sumw2();
	hExpPMTA->Sumw2();
	hExpPMTB->Sumw2();
	hHitsA->Sumw2();
	hHitsB->Sumw2();
	
	hExpC->Add(hExpA, hExpB, 1.0/timeA, -1.0/timeB);
	hExpSiPMC->Add(hExpSiPMA, hExpSiPMB, 1.0/timeA, -1.0/timeB);
	hExpPMTC->Add(hExpPMTA, hExpPMTB, 1.0/timeA, -1.0/timeB);
	hHitsC->Add(hHitsA, hHitsB, 1.0/timeA, -1.0/timeB);

	hXY->GetXaxis()->SetLabelSize(0.045);
	hXY->GetYaxis()->SetLabelSize(0.045);
	hXY->GetZaxis()->SetLabelSize(0.045);
	hExpC->GetYaxis()->SetLabelSize(0.05);
	hExpC->SetLineWidth(2);
	hExpC->SetLineColor(kBlue);
	hExpSiPMC->GetYaxis()->SetLabelSize(0.05);
	hExpPMTC->GetYaxis()->SetLabelSize(0.05);
	hHitsC->GetYaxis()->SetLabelSize(0.05);
	hHitsC->SetLineColor(kBlue);
	hHitsC->SetMarkerColor(kBlue);
	hHitsC->SetMarkerStyle(kFullSquare);

	hHitsC->SetStats(0);
	
	TCanvas *cExpA = new TCanvas("cExpA", "Data", 1200, 800);
	cExpA->Divide(3, 1);
	cExpA->cd(1);
	hExpC->Fit("gaus", "", "", Efit[0], Efit[1]);
	cExpA->cd(2);
	hExpSiPMC->Fit("gaus", "", "", Efit[0], Efit[1]);
	cExpA->cd(3);
	hExpPMTC->Fit("gaus", "", "", Efit[0], Efit[1]);

	sprintf(str, "src_v3/%s.pdf(", fname);
	cExpA->SaveAs(str);

	TCanvas *cExpB = new TCanvas("cExpB", "Hits", 1200, 800);
	cExpB->Divide(2, 1);
	cExpB->cd(1);
	hXY->Draw("colz");
	cExpB->cd(2);
	hHitsC->Draw();
	sprintf(str, "src_v3/%s.pdf)", fname);
	cExpB->SaveAs(str);

	sprintf(str, "src_v3/%s.root", fname);
	TFile *f = new TFile(str, "RECREATE");
	if (f->IsOpen()) {
		f->cd();
		hXY->Write();
		hExpA->Write();
		hExpB->Write();
		hExpC->Write();
		hExpSiPMA->Write();
		hExpSiPMB->Write();
		hExpSiPMC->Write();
		hExpPMTA->Write();
		hExpPMTB->Write();
		hExpPMTC->Write();
		hHitsA->Write();
		hHitsB->Write();
		hHitsC->Write();
		f->Close();
	}
}

void draw_MC(TChain *tMc, const char *name, const char *fname, TCut cXY, TCut cZ, double Efit[2], double kSP = 0.5, double kRndm = 0.0)
{
	char str[256];
	double rAB;
	long NA, NB;
	
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
//	gStyle->SetOptStat(0);
//	gStyle->SetOptFit(0);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
//	gStyle->SetLineWidth(4);
	
	sprintf(str, "Monte Carlo energy deposit in %s decay;E, MeV", name);
	TH1D *hMc = new TH1D("hMc", str, 70, 0, 7);
	sprintf(str, "Monte Carlo SiPM energy deposit in %s decay;E, MeV", name);
	TH1D *hMcSiPM = new TH1D("hMcSiPM", str, 35, 0, 7);
	sprintf(str, "Monte Carlo PMT energy deposit in %s decay;E, MeV", name);
	TH1D *hMcPMT = new TH1D("hMcPMT", str, 35, 0, 7);
	sprintf(str, "Monte Carlo number of hits from %s decay", name);
	TH1D *hMcHits = new TH1D("hMcHits", str, 20, 0, 20);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
	TCut cn("SiPmCleanHits > 5");
	TCut cSel = cxyz && cVeto;
	
	sprintf(str, "MyRandom::GausAdd(%6.4f*SiPmCleanEnergy+%6.4f*PmtCleanEnergy, %6.4f)", kSP, 1-kSP, kRndm);
	tMc->Project("hMc", str, cSel && cZ && cXY && cn);
	sprintf(str, "MyRandom::GausAdd(SiPmCleanEnergy, %6.4f)", kRndm);
	tMc->Project("hMcSiPM", str, cSel && cZ && cXY && cn);
	sprintf(str, "MyRandom::GausAdd(PmtCleanEnergy, %6.4f)", kRndm);
	tMc->Project("hMcPMT", str, cSel && cZ && cXY && cn);
	tMc->Project("hMcHits", "SiPmCleanHits", cSel && cZ && cXY);
	
	hMc->Sumw2();
	hMcSiPM->Sumw2();
	hMcPMT->Sumw2();
	hMcHits->Sumw2();
	
	hMc->GetYaxis()->SetLabelSize(0.05);
	hMcSiPM->GetYaxis()->SetLabelSize(0.05);
	hMcPMT->GetYaxis()->SetLabelSize(0.05);
	hMcHits->GetYaxis()->SetLabelSize(0.05);
	hMcHits->SetLineColor(kRed);
	hMcHits->SetMarkerColor(kRed);
	hMcHits->SetMarkerStyle(kFullCircle);
	hMcHits->SetStats(0);
	
	TCanvas *cMc = new TCanvas("cMc", "Monte Carlo", 1200, 800);
	cMc->Divide(2, 2);
	cMc->cd(1);
	hMc->Fit("gaus", "", "", Efit[0], Efit[1]);
	cMc->cd(2);
	hMcSiPM->Fit("gaus", "", "", Efit[0], Efit[1]);
	cMc->cd(3);
	hMcPMT->Fit("gaus", "", "", Efit[0], Efit[1]);
	cMc->cd(4);
	hMcHits->Draw("hist");
	sprintf(str, "src_v3/%s.pdf", fname);
	cMc->SaveAs(str);
	
	sprintf(str, "src_v3/%s.root", fname);
	TFile *f = new TFile(str, "RECREATE");
	if (f->IsOpen()) {
		f->cd();
		hMc->Write();
		hMcSiPM->Write();
		hMcPMT->Write();
		hMcHits->Write();
		f->Close();
	}
}

void draw_Sources(int iser, double kSP = 0.5, double kRndm = 0.0)
{
	const char *name;
	char fname[1024];
	char str[1024];
	double Efit[2] = {1.5, 3.5};
	TCut cXY;
	TCut cZ;
	int code;
	int i;
	TChain *tMc = new TChain("DanssEvent");
	TChain *tExpA = new TChain("DanssEvent");
	TChain *tExpB = new TChain("DanssEvent");
	TChain *tInfoA = new TChain("DanssInfo");
	TChain *tInfoB = new TChain("DanssInfo");
	
	code = iser / 1000;
	switch (iser) {
	case 2:		// Na Feb 17, center, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012380.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012381.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012380.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012381.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		name = "22Na";
		sprintf(fname, "22Na_feb17_short_old_center_ksp_%4.2f", kSP);
		break;
	case 4:		// Na May 17, center, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020243.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020244.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_020243.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_020244.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		name = "22Na";
		sprintf(fname, "22Na_may17_old_center_ksp_%4.2f", kSP);
		break;
	case 8:		// Na Feb 17, center, all files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12376; i<=12407; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "22Na";
		sprintf(fname, "22Na_feb17_old_center_ksp_%4.2f", kSP);
		break;
	case 12:	// Na Feb 17, Y90cm, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012370.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012371.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012370.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012371.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		name = "22Na";
		sprintf(fname, "22Na_feb17_short_old_90cm_ksp_%4.2f", kSP);
		break;
	case 18:	// Na Feb 17, Y90cm, all files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12364; i<=12373; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "22Na";
		sprintf(fname, "22Na_feb17_old_90cm_ksp_%4.2f", kSP);
		break;
	case 22:	// Na Feb 17, center, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012380.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012381.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012380.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012381.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		name = "22Na";
		sprintf(fname, "22Na_feb17_short_new_center_ksp_%4.2f", kSP);
		break;
	case 24:	// Na May 17, center, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020243.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020244.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020252.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020253.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_020243.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_020244.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_020252.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_020253.root");
		name = "22Na";
		sprintf(fname, "22Na_may17_new_center_ksp_%4.2f", kSP);
		break;
	case 28:	// Na Feb 17, center, all files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12376; i<=12407; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "22Na";
		sprintf(fname, "22Na_feb17_new_center_ksp_%4.2f", kSP);
		break;
	case 32:	// Na Feb 17, Y90cm, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012370.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012371.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012370.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012371.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		name = "22Na";
		sprintf(fname, "22Na_feb17_short_new_90cm_ksp_%4.2f", kSP);
		break;
	case 38:	// Na Feb 17, Y90cm, all files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12364; i<=12373; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "22Na";
		sprintf(fname, "22Na_feb17_new_90cm_ksp_%4.2f", kSP);
		break;
	case 102:	// Co Feb 17, center, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012310.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012311.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012310.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012311.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		name = "60Co";
		sprintf(fname, "60Co_feb17_short_old_center_ksp_%4.2f", kSP);
		break;
	case 104:	// Co May 17, center, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020233.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020234.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_020233.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_020234.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		name = "60Co";
		sprintf(fname, "60Co_may17_old_center_ksp_%4.2f", kSP);
		break;
	case 108:	// Co Feb 17, center, all files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12306; i<=12346; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "60Co";
		sprintf(fname, "60Co_feb17_old_center_ksp_%4.2f", kSP);
		break;
	case 112:	// Co Feb 17, Y90cm, 2 files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012350.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012351.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012350.root");
		tInfoA->AddFile("/mnt/root0/danss_root4/danss_012351.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tInfoB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		name = "60Co";
		sprintf(fname, "60Co_feb17_short_old_90cm_ksp_%4.2f", kSP);
		break;
	case 118:	// Co Feb 17, Y90cm, all files, old analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12348; i<=12361; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "60Co";
		sprintf(fname, "60Co_feb17_old_90cm_ksp_%4.2f", kSP);
		break;
	case 122:	// Co Feb 17, center, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012310.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012311.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012310.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012311.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		name = "60Co";
		sprintf(fname, "60Co_feb17_short_new_center_ksp_%4.2f", kSP);
		break;
	case 124:	// Co May 17, center, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020233.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020234.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020252.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020253.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_020233.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_020234.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_020252.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_020253.root");
		name = "60Co";
		sprintf(fname, "60Co_may17_new_center_ksp_%4.2f", kSP);
		break;
	case 128:	// Co Feb 17, center, all files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12306; i<=12346; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "60Co";
		sprintf(fname, "60Co_feb17_new_center_ksp_%4.2f", kSP);
		break;
	case 132:	// Co Feb 17, Y90cm, 2 files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012350.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012351.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012350.root");
		tInfoA->AddFile("/mnt/root1/danss_root5/danss_012351.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tInfoB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		name = "60Co";
		sprintf(fname, "60Co_feb17_short_new_90cm_ksp_%4.2f", kSP);
		break;
	case 138:	// Co Feb 17, Y90cm, all files, new analysys
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		for (i=12348; i<=12361; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
			tInfoA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
			tInfoB->AddFile(str);
		}
		name = "60Co";
		sprintf(fname, "60Co_feb17_new_90cm_ksp_%4.2f", kSP);
		break;
	case 1000:	// Na MC, center
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tMc->AddFile("/mnt/root0/danss_root4/mc_22Na_center_transcode.root");
		name = "22Na";
		sprintf(fname, "22Na_MC_center_ksp_%4.2f_rndm_%4.2f", kSP, kRndm);
		break;
	case 1010:	// Na MC, edge
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tMc->AddFile("/mnt/root0/danss_root4/mc_22Na_90cm_transcode.root");
		name = "22Na";
		sprintf(fname, "22Na_MC_90cm_ksp_%4.2f_rndm_%4.2f", kSP, kRndm);
		Efit[0] = 1.3;
		break;
	case 1100:	// Co MC, center
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tMc->AddFile("/mnt/root0/danss_root4/MC_newGeo_G/mc_60Co_center_transcode.root");
		name = "60Co";
		sprintf(fname, "60Co_MC_center_ksp_%4.2f_rndm_%4.2f", kSP, kRndm);
		break;
	case 1110:	// Co MC, edge
		cXY = (TCut) "(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 88) * (NeutronX[1] - 88) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400";
		cZ = (TCut) "(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100";
		tMc->AddFile("/mnt/root0/danss_root4/mc_60Co_90cm_transcode.root");
		name = "60Co";
		sprintf(fname, "60Co_MC_90cm_ksp_%4.2f_rndm_%4.2f", kSP, kRndm);
		break;
	default:
		printf("%d - unknown\n", iser);
		printf("Available valuse for iser = MIDD:\n");
		printf("M - MC (1) or experiment (0);\n");
		printf("I - isotope: 60Co (1) or 22Na (0);\n");
		printf("DD - serial for experiment:\n");
		printf("2  - February 17, center, old analysis, small dataset\n");
		printf("4  - May 17, center, old analysis\n");
		printf("8  - February 17, center, old analysis, large dataset\n");
		printf("12 - February 17, 90cm, old analysis, small dataset\n");
//		printf("14 - May 17, 90cm, old analysis\n");
		printf("18 - February 17, 90cm, old analysis, large dataset\n");
		printf("22 - February 17, center, new analysis, small dataset\n");
		printf("24 - May 17, center, new analysis\n");
		printf("28 - February 17, center, new analysis, large dataset\n");
		printf("32 - February 17, 90cm, new analysis, small dataset\n");
//		printf("34 - May 17, 90cm, new analysis\n");
		printf("38 - February 17, 90cm, new analysis, large dataset\n");
		printf("DD - for MC:\n");
		printf("0  - center (50, 50, 50) position\n");
		printf("10 - Y90cm (50, 90, 50) position\n");
		code = -1;
		break;
	}

	switch (code) {
	case 0:	// Experiment
		draw_Exp(tExpA, tExpB, tInfoA, tInfoB, name, fname, cXY, cZ, Efit, kSP);
		break;
	case 1:	// MC
		draw_MC(tMc, name, fname, cXY, cZ, Efit, kSP, kRndm);
		break;		
	default:
		break;
	}

	delete tMc;
	delete tExpA;
	delete tExpB;
}

void src2mc(const char *expfile, const char *mcfile, const char *newname)
{
	int i;

	TFile *fExp = new TFile(expfile);
	TFile *fMC = new TFile(mcfile);
	TH1 *hExpE = (TH1 *) fExp->Get("hExpC");
	TH1 *hExpH = (TH1 *) fExp->Get("hHitsC");
	TH1 *hMCE = (TH1 *) fMC->Get("hMc");
	TH1 *hMCH = (TH1 *) fMC->Get("hMcHits");
	if (!(hExpE && hExpH && hMCE && hMCH)) {
		printf("Something is wrong.\n");
		return;
	}
	hExpE->RecursiveRemove(hExpE->GetFunction("gaus"));
	hMCE->RecursiveRemove(hMCE->GetFunction("gaus"));
	TCanvas *cv = new TCanvas("CV", "SRC", 1200, 900);
	cv->Divide(2, 1);
	TVirtualPad * pd = cv->cd(1);
	pd->SetLeftMargin(0.21);
	pd->SetRightMargin(0.03);
	hExpE->SetLineColor(kBlack);
	hExpE->SetLineWidth(3);
	hExpE->SetStats(0);
	hExpE->GetXaxis()->SetTitle("Energy, MeV");
	hExpE->GetYaxis()->SetTitle("Events/100 keV");
	hExpE->GetYaxis()->SetTitleOffset(2.0);
	hExpE->Draw("e");
	hMCE->Scale(hExpE->Integral() / hMCE->Integral());
	hMCE->SetLineColor(kBlue);
	hMCE->SetStats(0);
	hMCE->Draw("same,hist");
	TLegend *lg = new TLegend(0.65, 0.8, 0.97, 0.9);
	lg->AddEntry(hExpE, "Experiment", "LE");
	lg->AddEntry(hMCE, "MC", "L");
	lg->Draw();

	pd = cv->cd(2);
	pd->SetLeftMargin(0.13);
	pd->SetRightMargin(0.03);
	hMCH->Scale(hExpH->Integral() / hMCH->Integral());
	hMCH->SetLineColor(kBlue);
	hMCH->SetStats(0);
	hMCH->GetXaxis()->SetTitle("SiPM hits");
	hMCH->GetYaxis()->SetTitle("Events/100 keV");
	hMCH->GetYaxis()->SetTitleOffset(1.3);
	hMCH->Draw("hist");
	hExpH->SetLineColor(kBlack);
	hExpH->SetLineWidth(3);
	hExpH->SetStats(0);
	hExpH->Draw("e,same");
	lg->Draw();
	cv->Update();
	cv->SaveAs(newname);
	
	printf("%s - experiment energy:\nEnergy     Events\n", newname);
	for (i=0; i<hExpE->GetNbinsX(); i++) printf("%5.3f-%5.3f    %f +- %f\n",
		hExpE->GetBinLowEdge(i+1), hExpE->GetBinLowEdge(i+1) + hExpE->GetBinWidth(i+1),
		hExpE->GetBinContent(i+1), hExpE->GetBinError(i+1));
	printf("%s - experiment SiPM hits:\nN     Events\n", newname);
	for (i=0; i<hExpH->GetNbinsX(); i++) printf("%d    %f +- %f\n", i,
		hExpH->GetBinContent(i+1), hExpH->GetBinError(i+1));
	printf("%s - MC energy:\nEnergy     Events\n", newname);
	for (i=0; i<hMCE->GetNbinsX(); i++) printf("%5.3f-%5.3f    %f +- %f\n",
		hMCE->GetBinLowEdge(i+1), hMCE->GetBinLowEdge(i+1) + hMCE->GetBinWidth(i+1),
		hMCE->GetBinContent(i+1), hMCE->GetBinError(i+1));
	printf("%s - MC SiPM hits:\nN     Events\n", newname);
	for (i=0; i<hMCH->GetNbinsX(); i++) printf("%d    %f +- %f\n", i,
		hMCH->GetBinContent(i+1), hMCH->GetBinError(i+1));

	fExp->Close();
	fMC->Close();
}
