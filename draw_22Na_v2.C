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


void draw_Src(TChain *tMc, TChain *tExpA, TChain *tExpB, const char *name, const char *fname, double kSP = 0.5, double kRndm = 0.0)
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
	sprintf(str, "Monte Carlo energy deposit in %s decay with random %2.0f%%;E, MeV", name, kRndm*100);
	TH1D *hMcR = new TH1D("hMcR", str, 70, 0, 7);
	sprintf(str, "Monte Carlo number of hits from %s decay", name);
	TH1D *hMcHits = new TH1D("hMcHits", str, 20, 0, 20);
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
	
	TH1D *hTmpA = new TH1D("hTmpA", "Normalization counts A", 100, 0, 1000);
	TH1D *hTmpB = new TH1D("hTmpB", "Normalization counts B", 100, 0, 1000);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cz50("(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400");
	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
	TCut cn("SiPmCleanHits > 5");
	
	sprintf(str, "%6.4f*SiPmCleanEnergy+%6.4f*PmtCleanEnergy", kSP, 1-kSP);
	tMc->Project("hMc", str, cxyz && ccc && cVeto && cn);
	tMc->Project("hMcSiPM", "SiPmCleanEnergy", cxyz && ccc && cVeto && cn);
	tMc->Project("hMcPMT", "PmtCleanEnergy", cxyz && ccc && cVeto && cn);
	sprintf(str, "MyRandom::GausAdd(%6.4f*SiPmCleanEnergy+%6.4f*PmtCleanEnergy, %6.4f)", kSP, 1-kSP, kRndm);
	tMc->Project("hMcR", str, cxyz && ccc && cVeto && cn);
	tMc->Project("hMcHits", "SiPmCleanHits", cxyz && ccc && cVeto);
	
	sprintf(str, "%6.4f*SiPmCleanEnergy+%6.4f*PmtCleanEnergy", kSP, 1-kSP);
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cxyz && cz50 && cVeto && cn);
	tExpA->Project("hExpA", str, cxyz && ccc && cVeto && cn);
	tExpB->Project("hExpB", str, cxyz && ccc && cVeto && cn);
	tExpA->Project("hExpSiPMA", "SiPmCleanEnergy", cxyz && ccc && cVeto && cn);
	tExpB->Project("hExpSiPMB", "SiPmCleanEnergy", cxyz && ccc && cVeto && cn);
	tExpA->Project("hExpPMTA", "PmtCleanEnergy", cxyz && ccc && cVeto && cn);
	tExpB->Project("hExpPMTB", "PmtCleanEnergy", cxyz && ccc && cVeto && cn);
	tExpA->Project("hHitsA", "SiPmCleanHits", cxyz && ccc && cVeto);
	tExpB->Project("hHitsB", "SiPmCleanHits", cxyz && ccc && cVeto);
	NA = tExpA->Project("hTmpA", "SiPmCleanEnergy", "(SiPmCleanEnergy + PmtCleanEnergy) / 2 > 100");
	NB = tExpB->Project("hTmpB", "SiPmCleanEnergy", "(SiPmCleanEnergy + PmtCleanEnergy) / 2 > 100");
	
	rAB = 1.0 * NA / NB;
	printf("NA = %ld    NB = %ld    rAB = %f\n", NA, NB, rAB);
	
	hMc->Sumw2();
	hMcSiPM->Sumw2();
	hMcPMT->Sumw2();
	hMcR->Sumw2();
	hMcHits->Sumw2();
	hExpA->Sumw2();
	hExpB->Sumw2();
	hExpSiPMA->Sumw2();
	hExpSiPMB->Sumw2();
	hExpPMTA->Sumw2();
	hExpPMTB->Sumw2();
	hHitsA->Sumw2();
	hHitsB->Sumw2();
	
	hExpC->Add(hExpA, hExpB, 1.0, -rAB);
	hExpSiPMC->Add(hExpSiPMA, hExpSiPMB, 1.0, -rAB);
	hExpPMTC->Add(hExpPMTA, hExpPMTB, 1.0, -rAB);
	hHitsC->Add(hHitsA, hHitsB, 1.0, -rAB);
	hMcHits->Scale(hHitsC->Integral() / hMcHits->Integral());
	hMcR->Scale(hExpC->Integral() / hMcR->Integral());

	hMc->GetYaxis()->SetLabelSize(0.05);
	hMcSiPM->GetYaxis()->SetLabelSize(0.05);
	hMcPMT->GetYaxis()->SetLabelSize(0.05);
	hMcR->GetYaxis()->SetLabelSize(0.05);
	hMcR->SetLineColor(kRed);
	hMcHits->GetYaxis()->SetLabelSize(0.05);
	hMcHits->SetLineColor(kRed);
	hMcHits->SetMarkerColor(kRed);
	hMcHits->SetMarkerStyle(kFullCircle);
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

	hMcHits->SetStats(0);
	hHitsC->SetStats(0);
	
	TLegend *lg = new TLegend(0.65, 0.8, 0.95, 0.93);
	lg->AddEntry(hMcHits, "Monte Carlo", "L");
	lg->AddEntry(hHitsC,  "DANSS", "LP");
	lg->SetTextSize(0.035);
	
	TCanvas *cMc = new TCanvas("cMc", "Monte Carlo", 1200, 800);
	cMc->Divide(2, 2);
	cMc->cd(1);
	hMc->Fit("gaus", "", "", 1.5, 3.5);
	cMc->cd(2);
	hMcSiPM->Fit("gaus", "", "", 1.5, 3.5);
	cMc->cd(3);
	hMcPMT->Fit("gaus", "", "", 1.5, 3.5);
	cMc->cd(4);
	hMcR->Fit("gaus", "", "", 1.5, 3.5);
	sprintf(str, "%s.pdf(", fname);
	cMc->SaveAs(str);
	
	TCanvas *cExp = new TCanvas("cExp", "Data", 1200, 800);
	cExp->Divide(2, 2);
	cExp->cd(1);
	hExpC->Fit("gaus", "", "", 1.5, 3.5);
	cExp->cd(2);
	hExpSiPMC->Fit("gaus", "", "", 1.5, 3.5);
	cExp->cd(3);
	hExpPMTC->Fit("gaus", "", "", 1.5, 3.5);
	cExp->cd(4);
	hXY->Draw("colz");
	sprintf(str, "%s.pdf", fname);
	cExp->SaveAs(str);

	TCanvas *cHits = new TCanvas("cHits", "Hits", 1200, 800);
	cHits->Divide(2, 1);
	cHits->cd(1);
	hHitsC->Draw();
	hMcHits->Draw("same,hist");
	lg->Draw();
	cHits->cd(2);
	hExpC->Draw();
	hMcR->Draw("same,hist");
	cHits->Update();
	sprintf(str, "%s.pdf)", fname);
	cHits->SaveAs(str);

	sprintf(str, "%s.root", fname);
	TFile *f = new TFile(str, "RECREATE");
	if (f->IsOpen()) {
		f->cd();
		hMc->Write();
		hMcSiPM->Write();
		hMcPMT->Write();
		hMcR->Write();
		hMcHits->Write();
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
	
	TCanvas *cPRL = new TCanvas("PRL", "PRL", 800, 800);
	cPRL->SetLeftMargin(0.17);
	cPRL->SetRightMargin(0.03);
	cPRL->SetTopMargin(0.03);
	cPRL->SetBottomMargin(0.10);
	hExpC->SetLineColor(kBlack);
	hExpC->GetXaxis()->SetRange(0, 50);
	hExpC->SetTitle(";E, MeV;Events/100 keV");
	hExpC->SetStats(0);
	hExpC->GetYaxis()->SetTitleOffset(1.7);
	hExpC->Draw();
	hMcR->Draw("hits,same");
	lg->Draw();
	sprintf(str, "%s-prl.pdf", fname);
	cPRL->SaveAs(str);
}

void draw_22Na(int iType = -1, double kSP = 0.5, double kRndm = 0.0)
{
	char str[128];
	int i;
	
	TChain *tMc = new TChain("DanssEvent");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_transcodeNew.root");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_simple.root");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_newCuts_transcodeNew.root");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_dl200um_transcodeNew.root");
	tMc->AddFile("/mnt/root0/danss_root4/mc_22Na_center_transcode.root");
//	tMc->AddFile("/mnt/root0/danss_root4/MC_newGeo_G/mc_22Na_stripHeterogenity_transcode.root");
	TChain *tExpA = new TChain("DanssEvent");
	TChain *tExpB = new TChain("DanssEvent");
	switch (iType) {
	case 0:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002197.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002198.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002285.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002286.root");
		break;
	case 1:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002197_000.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002198_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002285_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002286_000.root");
		break;
	case 2:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012380.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012381.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		break;
	case 4:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020243.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020244.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		break;
	case 5:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020243_000.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020244_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020252_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020253_000.root");
		break;
	case 6:
		tExpA->AddFile("../digi.v2/danss_root4_2/danss_020243.root");
		tExpA->AddFile("../digi.v2/danss_root4_2/danss_020244.root");
		tExpB->AddFile("../digi.v2/danss_root4_2/danss_020252.root");
		tExpB->AddFile("../digi.v2/danss_root4_2/danss_020253.root");
		break;
	case 8:
		for (i=12376; i<=12407; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
		}
		break;
	case 16:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020252.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020253.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020254.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020255.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020256.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020257.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020258.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020259.root");
		break;
	case 52:
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012380.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012381.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		break;
	case 54:
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020243.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020244.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020252.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020253.root");
		break;
	case 58:
		for (i=12376; i<=12407; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
		}
		break;
	case 60:
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012364.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012365.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		break;
	default:
		printf("Types: \n");
		printf("0  - April 16, fast, old analysis\n");
		printf("1  - April 16, full, old analysis\n");
		printf("2  - February 17, fast, old analysis, small dataset\n");
		printf("4  - May 18, fast, old analysis\n");
		printf("6  - May 18, fast, new analysis prerelease\n");
		printf("8  - February 17, fast, old analysis, large dataset\n");
		printf("16 - May 18, background to background\n");
		printf("52 - February 17, fast, new analysis, small dataset\n");
		printf("54 - May 18, fast, new analysis\n");
		printf("58 - February 17, fast, new analysis, large dataset\n");
		printf("Other types are unknown.\n");
		return;
	}
	sprintf(str, "src/22Na_v2_%d-%5.3f-%5.3f", iType, kSP, kRndm);
	draw_Src(tMc, tExpA, tExpB, "^{22}Na", str, kSP, kRndm);
}

void draw_60Co(int iType = -1, double kSP = 0.5, double kRndm = 0.0)
{
	char str[128];
	int i;
	
	TChain *tMc = new TChain("DanssEvent");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_transcodeNew.root");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_newCuts_transcodeNew.root");
//	tMc->AddFile("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_dl200um_transcodeNew.root");
	tMc->AddFile("/mnt/root0/danss_root4/MC_newGeo_G/mc_60Co_center_transcode.root");
//	tMc->AddFile("/mnt/root0/danss_root4/MC_newGeo_G/mc_60Co_stripHeterogenity_transcode.root");
	TChain *tExpA = new TChain("DanssEvent");
	TChain *tExpB = new TChain("DanssEvent");
	switch (iType) {
	case 0:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002106.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002107.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002285.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002286.root");
		break;
	case 1:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002106_000.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_002107_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002285_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_002286_000.root");
		break;
	case 2:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012310.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_012311.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_012302.root");
		break;
	case 4:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020233.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020234.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020230.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020231.root");
		break;
	case 5:
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020233_000.root");
		tExpA->AddFile("/mnt/root0/danss_root4/danss_020234_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020252_000.root");
		tExpB->AddFile("/mnt/root0/danss_root4/danss_020253_000.root");
		break;
	case 6:
		tExpA->AddFile("../digi.v2/danss_root4_2/danss_020233.root");
		tExpA->AddFile("../digi.v2/danss_root4_2/danss_020234.root");
		tExpB->AddFile("../digi.v2/danss_root4_2/danss_020252.root");
		tExpB->AddFile("../digi.v2/danss_root4_2/danss_020253.root");
		break;
	case 8:
		for (i=12306; i<=12346; i++) {
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpA->AddFile(str);
//			printf("Signal: %s\n", str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root0/danss_root4/danss_%6.6d.root", i);
			tExpB->AddFile(str);
//			printf("Background: %s\n", str);
		}
		break;
	case 52:
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012310.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012311.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		break;
	case 54:
		tExpA->AddFile("/mnt/root1/danss_root4/danss_020233.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_020234.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020230.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_020231.root");
		break;
	case 58:
		for (i=12306; i<=12346; i++) {
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpA->AddFile(str);
		}
		for (i=12212; i<=12304; i++) {
			if (i == 12240) continue;
			sprintf(str, "/mnt/root1/danss_root5/danss_%6.6d.root", i);
			tExpB->AddFile(str);
		}
		break;
	case 60:
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012360.root");
		tExpA->AddFile("/mnt/root1/danss_root5/danss_012361.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012301.root");
		tExpB->AddFile("/mnt/root1/danss_root5/danss_012302.root");
		break;
	default:
		printf("Types: \n");
		printf("0  - April 16, fast, old analysis\n");
		printf("1  - April 16, full, old analysis\n");
		printf("2  - February 17, fast, old analysis, small dataset\n");
		printf("4  - May 18, fast, old analysis\n");
		printf("6  - May 18, fast, new analysis prerelease\n");
		printf("8  - February 17, fast, old analysis, large dataset\n");
		printf("52 - February 17, fast, new analysis, small dataset\n");
		printf("54 - May 18, fast, new analysis\n");
		printf("58 - February 17, fast, new analysis, large dataset\n");
		printf("Other types are unknown.\n");
		return;
	}
	sprintf(str, "src/60Co_v2_%d-%5.3f-%5.3f", iType, kSP, kRndm);
	draw_Src(tMc, tExpA, tExpB, "^{60}Co", str, kSP, kRndm);
}

void draw_all(void)
{
	int i;
	for (i=0; i<6; i++) {
		draw_22Na(i);
		draw_60Co(i);
	}
}

void draw_Amps(void)
{
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	
	TChain *tNa = new TChain("DanssEvent");
	TChain *tCo = new TChain("DanssEvent");
	TChain *tBgnd = new TChain("DanssEvent");
	tNa->AddFile("../digi.v2/danss_root4_3/danss_020243.root");
	tNa->AddFile("../digi.v2/danss_root4_3/danss_020244.root");
	tCo->AddFile("../digi.v2/danss_root4_3/danss_020233.root");
	tCo->AddFile("../digi.v2/danss_root4_3/danss_020234.root");
	tBgnd->AddFile("../digi.v2/danss_root4_3/danss_020252.root");
	tBgnd->AddFile("../digi.v2/danss_root4_3/danss_020253.root");

	TH1D *hNaSiPMA = new TH1D("hNaSiPMA", "SiPM ^{22}Na;ADC units", 50, 0, 1500);
	TH1D *hNaSiPMB = new TH1D("hNaSiPMB", "SiPM ^{22}Na;ADC units", 50, 0, 1500);
	TH1D *hCoSiPMA = new TH1D("hCoSiPMA", "SiPM ^{60}Co;ADC units", 50, 0, 1500);
	TH1D *hCoSiPMB = new TH1D("hCoSiPMB", "SiPM ^{60}Co;ADC units", 50, 0, 1500);
	TH1D *hBgndSiPM = new TH1D("hBgndSiPM", "SiPM ^{60}Co;ADC units", 50, 0, 1500);

	TH1D *hNaPMTA = new TH1D("hNaPMTA", "PMT ^{22}Na;ADC units", 50, 0, 1000);
	TH1D *hNaPMTB = new TH1D("hNaPMTB", "PMT ^{22}Na;ADC units", 50, 0, 1000);
	TH1D *hCoPMTA = new TH1D("hCoPMTA", "PMT ^{60}Co;ADC units", 50, 0, 1000);
	TH1D *hCoPMTB = new TH1D("hCoPMTB", "PMT ^{60}Co;ADC units", 50, 0, 1000);
	TH1D *hBgndPMT = new TH1D("hBgndPMT", "PMT ^{60}Co;ADC units", 50, 0, 1000);
	
	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400");
	TCut cVeto("VetoHits < 2 && VetoEnergy < 4");
	
	tNa->Project("hNaSiPMA", "SiPmCleanEnergy", cxyz && ccc && cVeto);
	tCo->Project("hCoSiPMA", "SiPmCleanEnergy", cxyz && ccc && cVeto);
	tBgnd->Project("hBgndSiPM", "SiPmCleanEnergy", cxyz && ccc && cVeto);
	tNa->Project("hNaPMTA", "PmtCleanEnergy", cxyz && ccc && cVeto);
	tCo->Project("hCoPMTA", "PmtCleanEnergy", cxyz && ccc && cVeto);
	tBgnd->Project("hBgndPMT", "PmtCleanEnergy", cxyz && ccc && cVeto);
	
	hNaSiPMA->Sumw2();
	hCoSiPMA->Sumw2();
	hBgndSiPM->Sumw2();
	hNaPMTA->Sumw2();
	hCoPMTA->Sumw2();
	hBgndPMT->Sumw2();

	hNaSiPMB->Add(hNaSiPMA, hBgndSiPM, 1.0, -0.977);
	hCoSiPMB->Add(hCoSiPMA, hBgndSiPM, 1.0, -0.977);
	hNaPMTB->Add(hNaPMTA, hBgndPMT, 1.0, -0.977);
	hCoPMTB->Add(hCoPMTA, hBgndPMT, 1.0, -0.977);

	hNaSiPMB->GetYaxis()->SetLabelSize(0.05);
	hCoSiPMB->GetYaxis()->SetLabelSize(0.05);
	hNaPMTB->GetYaxis()->SetLabelSize(0.05);
	hCoPMTB->GetYaxis()->SetLabelSize(0.05);

	TCanvas *cAmp = new TCanvas("cAmp", "Raw Amplitude", 1200, 800);
	cAmp->Divide(2, 2);
	cAmp->cd(1);
	hNaSiPMB->Fit("gaus", "", "", 200, 800);
	cAmp->cd(2);
	hCoSiPMB->Fit("gaus", "", "", 200, 200);
	cAmp->cd(3);
	hNaPMTB->Fit("gaus", "", "", 150, 600);
	cAmp->cd(4);
	hCoPMTB->Fit("gaus", "", "", 150, 600);
	cAmp->SaveAs("22Na_60Co_raw_ADC_amplitudes.pdf");
}

void draw_graphs(void)
{
	const char title[4][20] = {"Na_mean;um", "Na_sigma;um", "Co_mean;um", "Co_sigma;um"};
	const double tgt[4] = {1.787, 0.394, 2.082, 0.473};
	const double mc[4][4] = {{2.008, 1.985, 1.974, 1.933}, {0.319, 0.326, 0.330, 0.336}, 
		{2.233, 2.211, 2.185, 2.149}, {0.377, 0.385, 0.394, 0.405}};
	const double dl[4] = {0, 50, 100, 200};
	const double mm[4][2] = {{1.7, 2.1}, {0.3, 0.4}, {2.0, 2.3}, {0.35, 0.5}};
	TGraph *gr;
	TLine *ln;
	TH1D *h;
	TLegend *lg;
	int i;
	
	TCanvas *cv = new TCanvas;
	cv->Divide(2,2);
	gr = new TGraph;
	ln = new TLine;
	gr->SetMarkerStyle(kFullCircle);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(kBlue);
	ln->SetLineWidth(4);
	ln->SetLineColor(kRed);
	lg = new TLegend(0.65, 0.65, 0.95, 0.85);
	lg->AddEntry(ln, "Exp.", "L");
	lg->AddEntry(gr, "MC", "P");
	for (i=0; i<4; i++) {
		cv->cd(i+1);
		
		h = new TH1D(title[i], title[i], 100, 0, 200);
		h->SetMinimum(mm[i][0]);
		h->SetMaximum(mm[i][1]);
		h->DrawCopy();
		gr->DrawGraph(4, dl, mc[i], "P");
		ln->DrawLine(0, tgt[i], 200, tgt[i]);
		lg->Draw();
		delete h;
	}
}
