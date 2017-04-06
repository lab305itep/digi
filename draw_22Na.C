void draw_22Na(void)
{
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetLineWidth(2);
	
	TFile *fMc = new TFile("/space/danss_root3/mc_22Na_simple.root");
	TTree *tMc = (TTree *) fMc->Get("DanssEvent");
	if (!tMc) {
		printf("Can not open MC tree.\n");
		return;
	}
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/space/danss_root3/danss_data_002197_000.root");
	tExpA->AddFile("/space/danss_root3/danss_data_002198_000.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/space/danss_root3/danss_data_002193_000.root");
	tExpB->AddFile("/space/danss_root3/danss_data_002194_000.root");
	
	TH1D *hMc = new TH1D("hMc", "Monte Carlo energy deposit in ^{22}Na decay;E, MeV", 50, 0, 5);
	TH1D *hMcHits = new TH1D("hMcHits", "Monte Carlo number of hits from ^{22}Na decay", 20, 0, 20);
	TH2D *hXY = new TH2D("hXY", "XY distribution of gamma flash center;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	TH1D *hExpA = new TH1D("hExpA", "DANSS energy deposit in ^{22}Na decay;E, MeV", 50, 0, 5);
	TH1D *hExpB = new TH1D("hExpB", "DANSS energy deposit in ^{22}Na decay;E, MeV", 50, 0, 5);
	TH1D *hExpC = new TH1D("hExpC", "DANSS energy deposit in ^{22}Na decay;E, MeV", 50, 0, 5);
	TH1D *hHitsA = new TH1D("hHitsA", "Number of hits from ^{22}Na decay", 20, 0, 20);
	TH1D *hHitsB = new TH1D("hHitsB", "Number of hits from ^{22}Na decay", 20, 0, 20);
	TH1D *hHitsC = new TH1D("hHitsC", "Number of hits from ^{22}Na decay", 20, 0, 20);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cz50("(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	TCut cc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 43) * (NeutronX[1] - 43) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	
	tMc->Project("hMc", "(SiPmCleanEnergy+PmtCleanEnergy)/2");
	tMc->Project("hMcHits", "SiPmCleanHits");
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cxyz && cz50);
	tExpA->Project("hExpA", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && cc);
	tExpB->Project("hExpB", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && cc);
	tExpA->Project("hHitsA", "SiPmCleanHits", cxyz && cc);
	tExpB->Project("hHitsB", "SiPmCleanHits", cxyz && cc);
	
	hExpA->Sumw2();
	hExpB->Sumw2();
	hMcHits->Sumw2();
	hHitsA->Sumw2();
	hHitsB->Sumw2();
	
//	hExpC->Add(hExpA, hExpB, 1.0, -hExpA->Integral(3, 8) / hExpB->Integral(3, 8));
	hExpC->Add(hExpA, hExpB, 1.0, -1.0);
	hHitsC->Add(hHitsA, hHitsB, 1.0, -1.0);
	hMcHits->Scale(hHitsC->Integral() / hMcHits->Integral());

	hMc->GetYaxis()->SetLabelSize(0.05);
	hMcHits->GetYaxis()->SetLabelSize(0.05);
	hXY->GetXaxis()->SetLabelSize(0.045);
	hXY->GetYaxis()->SetLabelSize(0.045);
	hXY->GetZaxis()->SetLabelSize(0.05);
	hExpC->GetYaxis()->SetLabelSize(0.05);
	hHitsC->GetYaxis()->SetLabelSize(0.05);
	hMcHits->SetLineColor(kRed);
	hHitsC->SetLineColor(kBlue);
	hMcHits->SetMarkerColor(kRed);
	hHitsC->SetMarkerColor(kBlue);
	hMcHits->SetMarkerStyle(kFullCircle);
	hHitsC->SetMarkerStyle(kFullSquare);

	hMcHits->SetStats(0);
	hHitsC->SetStats(0);
	
	TLegend *lg = new TLegend(0.8, 0.8, 0.95, 0.95);
	lg->AddEntry(hMcHits, "Monte Carlo", "L");
	lg->AddEntry(hHitsC,  "DANSS", "LP");
	
	TCanvas *cMc = new TCanvas("cMc", "Monte Carlo", 800, 1000);
	cMc->cd();
	hMc->Fit("gaus", "", "", 1.5, 2.7);
	
	TCanvas *cXY = new TCanvas("cXY", "DANSS XY", 800, 800);
	cXY->cd();
	hXY->Draw("colz");
	
	TCanvas *cExp = new TCanvas("cExp", "Danss", 800, 1000);
	cExp->cd();
	hExpC->Fit("gaus", "", "", 1.3, 2.7);
	
	TCanvas *cHits = new TCanvas("cHits", "Hits", 800, 1000);
	cHits->cd();
	hHitsC->Draw();
	hMcHits->Draw("same,hist");
	lg->Draw();
	cHits->Update();
}
