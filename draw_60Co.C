void draw_60Co(void)
{
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetLineWidth(2);
	
	TFile *fMc = new TFile("/space/danss_root3/mc_60Co_simple.root");
	TTree *tMc = (TTree *) fMc->Get("DanssEvent");
	if (!tMc) {
		printf("Can not open MC tree.\n");
		return;
	}
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/space/danss_root3/danss_data_002106_000.root");
	tExpA->AddFile("/space/danss_root3/danss_data_002107_000.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/space/danss_root3/danss_data_002102_000.root");
	tExpB->AddFile("/space/danss_root3/danss_data_002103_000.root");
	
	TH1D *hMc = new TH1D("hMc", "Monte Carlo energy deposit in ^{60}Co decay;E, MeV", 50, 0, 5);
	TH2D *hXY = new TH2D("hXY", "XY distribution of gamma flash center;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	TH1D *hExpA = new TH1D("hExpA", "DANSS energy deposit in ^{60}Co decay;E, MeV", 50, 0, 5);
	TH1D *hExpB = new TH1D("hExpB", "DANSS energy deposit in ^{60}Co decay;E, MeV", 50, 0, 5);
	TH1D *hExpC = new TH1D("hExpC", "DANSS energy deposit in ^{60}Co decay;E, MeV", 50, 0, 5);
	
	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cz50("(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	TCut cc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	
	tMc->Project("hMc", "(SiPmCleanEnergy+PmtCleanEnergy)/2");
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cxyz && cz50);
	tExpA->Project("hExpA", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && cc);
	tExpB->Project("hExpB", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && cc);
	
	hExpA->Sumw2();
	hExpB->Sumw2();
	
	hExpC->Add(hExpA, hExpB, 1.0, -hExpA->Integral(3, 8) / hExpB->Integral(3, 8));
	hExpC->Add(hExpA, hExpB, 1.0, -1.0);
	
	hMc->GetYaxis()->SetLabelSize(0.05);
	hXY->GetXaxis()->SetLabelSize(0.045);
	hXY->GetYaxis()->SetLabelSize(0.045);
	hXY->GetZaxis()->SetLabelSize(0.05);
	hExpC->GetYaxis()->SetLabelSize(0.05);
	
	TCanvas *cMc = new TCanvas("cMc", "Monte Carlo", 800, 1000);
	cMc->cd();
	hMc->Fit("gaus", "", "", 1.7, 3.0);
	
	TCanvas *cXY = new TCanvas("cXY", "DANSS XY", 800, 800);
	cXY->cd();
	hXY->Draw("colz");
	
	TCanvas *cExp = new TCanvas("cExp", "Danss", 800, 1000);
	cExp->cd();
	hExpC->Fit("gaus", "", "", 1.6, 3.0);
}