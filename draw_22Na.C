struct HitTypeStruct {
	char 	type;
	char	z;
	char	xy;
	char	flag;
};
struct DanssEventStruct4 {
//		Common parameters
	long long	globalTime;		// time in terms of 125 MHz
	long long	number;			// event number in the file
	int		runNumber;		// the run number
	int		unixTime;		// linux time, seconds
	float		fineTime;		// fine time of the event (for hit selection)
//		Veto parameters
	int		VetoHits;		// hits in the event record
	float		VetoEnergy;		// Energy Sum of all hits
	int		VetoCleanHits;		// hits above threshold and in time window
	float		VetoCleanEnergy;	// Energy Sum of clean hits
//		PMT parameters
	int		PmtHits;		// the same as above for PMT
	float		PmtEnergy;
	int		PmtCleanHits;
	float		PmtCleanEnergy;
//		SiPM parameters
	int		SiPmHits;		// the same as above for PMT
	float		SiPmEnergy;
	int		SiPmCleanHits;
	float		SiPmCleanEnergy;
	int		SiPmEarlyHits;		// to understand random background
	float		SiPmEarlyEnergy;
//		"positron cluster" parameters
	int		PositronHits;		// hits in the cluster
	int		PositronFlags;		// Positron flags
	float		PositronMinLen;		// Minimum track length to create the cluster
	float		PositronEnergy;		// Energy sum of the cluster, corrected (SiPM+PMT)
	float		TotalEnergy;		// Event full energy correctd (SiPM+PMT)
	float		PositronX[3];		// cluster position
	int		AnnihilationGammas;	// number of possible annihilation gammas
	float		AnnihilationEnergy;	// Energy in annihilation gammas
	float		AnnihilationMax;	// Energy in the maximum annihilation hit
//		"neutron" parameters
	float		NeutronX[3];		// center of gammas position
	int		NHits;			// Number of hits
} EVT;
struct HitStruct {
	float			E[2600];
	float			T[2600];
	struct HitTypeStruct 	type[2600];
}	HitArray;

void project_hits_distrib(TChain *chain, TH1D *hist)
{
	int i, j, N;
	double r2;
	
	chain->SetBranchAddress("Data", &EVT);
	chain->SetBranchAddress("HitE", HitArray.E);
	chain->SetBranchAddress("HitType", (int *)HitArray.type);
	N = chain->GetEntries();
	
	for (i=0; i<N; i++) {
		chain->GetEntry(i);
//	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
		if (EVT.NeutronX[0] < 0 || EVT.NeutronX[1] < 0 || EVT.NeutronX[2] < 0) continue;
//	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
		r2 = 0;
		for (j=0; j<3; j++) r2 += (EVT.NeutronX[0] - 48) * (EVT.NeutronX[0] - 48);
		if (r2 >= 100) continue;
		for (j=0; j<EVT.NHits; j++) if (HitArray.type[j].type == 0) hist->Fill(HitArray.E[j]);
	}
}

void project_mc_random_energy(TChain *chain, TH1D *hist)
{
	int i, j, N;
	double r2, Esum;
	
	TRandom2 *random = new TRandom2();
	chain->SetBranchAddress("Data", &EVT);
	chain->SetBranchAddress("HitE", HitArray.E);
	chain->SetBranchAddress("HitType", (int *)HitArray.type);
	N = chain->GetEntries();
	
	for (i=0; i<N; i++) {
		chain->GetEntry(i);
//	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
		if (EVT.NeutronX[0] < 0 || EVT.NeutronX[1] < 0 || EVT.NeutronX[2] < 0) continue;
//	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
		r2 = 0;
		for (j=0; j<3; j++) r2 += (EVT.NeutronX[0] - 48) * (EVT.NeutronX[0] - 48);
		if (r2 >= 100) continue;
		Esum = 0;
		for (j=0; j<EVT.NHits; j++) if (HitArray.type[j].type == 0 || HitArray.type[j].type == 1) Esum += HitArray.E[j] * random->Gaus(1.0, 0.26);
		hist->Fill(Esum/2);
	}
}

void draw_Src(TChain *tMc, TChain *tExpA, TChain *tExpB, double rAB, const char *name, const char *fname)
{
	char str[256];
	
	gStyle->SetOptStat("i");
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetLineWidth(2);
	
	sprintf(str, "Monte Carlo energy deposit in %s decay;E, MeV", name);
	TH1D *hMc = new TH1D("hMc", str, 50, 0, 5);
	sprintf(str, "Monte Carlo energy deposit in %s decay with random;E, MeV", name);
	TH1D *hMcR = new TH1D("hMcR", str, 50, 0, 5);
	sprintf(str, "Monte Carlo number of hits from %s decay", name);
	TH1D *hMcHits = new TH1D("hMcHits", str, 20, 0, 20);
	TH1D *hMcE = new TH1D("hMcE", "Monte Carlo hit energy;E, MeV", 60, 0, 3);
	TH2D *hXY = new TH2D("hXY", "XY distribution of gamma flash center;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	sprintf(str, "DANSS energy deposit in %s decay;E, MeV", name);
	TH1D *hExpA = new TH1D("hExpA", str, 50, 0, 5);
	TH1D *hExpB = new TH1D("hExpB", str, 50, 0, 5);
	TH1D *hExpC = new TH1D("hExpC", str, 50, 0, 5);
	sprintf(str, "Number of hits from %s decay", name);
	TH1D *hHitsA = new TH1D("hHitsA", str, 20, 0, 20);
	TH1D *hHitsB = new TH1D("hHitsB", str, 20, 0, 20);
	TH1D *hHitsC = new TH1D("hHitsC", str, 20, 0, 20);
	TH1D *hEA = new TH1D("hEA", "Hit energy;E, MeV", 60, 0, 3);
	TH1D *hEB = new TH1D("hEB", "Hit energy;E, MeV", 60, 0, 3);
	TH1D *hEC = new TH1D("hEC", "Hit energy;E, MeV", 60, 0, 3);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cz50("(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	
	tMc->Project("hMc", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc);
	tMc->Project("hMcHits", "SiPmCleanHits", cxyz && ccc);
	project_hits_distrib(tMc, hMcE);
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cxyz && cz50);
	tExpA->Project("hExpA", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc);
	tExpB->Project("hExpB", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc);
	tExpA->Project("hHitsA", "SiPmCleanHits", cxyz && ccc);
	tExpB->Project("hHitsB", "SiPmCleanHits", cxyz && ccc);
	project_hits_distrib(tExpA, hEA);
	project_hits_distrib(tExpB, hEB);
	project_mc_random_energy(tMc, hMcR);
	
	hMc->Sumw2();
	hMcR->Sumw2();
	hMcHits->Sumw2();
	hMcE->Sumw2();
	hExpA->Sumw2();
	hExpB->Sumw2();
	hHitsA->Sumw2();
	hHitsB->Sumw2();
	hEA->Sumw2();
	hEB->Sumw2();
	
	hExpC->Add(hExpA, hExpB, 1.0, -rAB);
	hHitsC->Add(hHitsA, hHitsB, 1.0, -rAB);
	hEC->Add(hEA, hEB, 1.0, -rAB);
	hMcHits->Scale(hHitsC->Integral() / hMcHits->Integral());
	hMcE->Scale(hEC->Integral() / hMcE->Integral());

	hMc->GetYaxis()->SetLabelSize(0.05);
	hMcR->GetYaxis()->SetLabelSize(0.05);
	hMcHits->GetYaxis()->SetLabelSize(0.05);
	hMcE->GetYaxis()->SetLabelSize(0.05);
	hXY->GetXaxis()->SetLabelSize(0.045);
	hXY->GetYaxis()->SetLabelSize(0.045);
	hXY->GetZaxis()->SetLabelSize(0.05);
	hExpC->GetYaxis()->SetLabelSize(0.05);
	hEC->GetYaxis()->SetLabelSize(0.05);
	hHitsC->GetYaxis()->SetLabelSize(0.05);
	hMcHits->SetLineColor(kRed);
	hHitsC->SetLineColor(kBlue);
	hMcHits->SetMarkerColor(kRed);
	hHitsC->SetMarkerColor(kBlue);
	hMcHits->SetMarkerStyle(kFullCircle);
	hHitsC->SetMarkerStyle(kFullSquare);
	hMcE->SetLineColor(kRed);
	hEC->SetLineWidth(2);
	hEC->SetMarkerStyle(kFullSquare);
	hEC->SetLineColor(kBlue);
	hEC->SetMarkerColor(kBlue);

	hMcHits->SetStats(0);
	hHitsC->SetStats(0);
	hMcE->SetStats(0);
	hEC->SetStats(0);
	
	TLegend *lg = new TLegend(0.65, 0.8, 0.95, 0.93);
	lg->AddEntry(hMcHits, "Monte Carlo", "L");
	lg->AddEntry(hHitsC,  "DANSS", "LP");
	lg->SetTextSize(0.035);
	
	TCanvas *cMc = new TCanvas("cMc", "Monte Carlo", 800, 1000);
	cMc->cd();
	hMc->Fit("gaus", "", "", 1.5, 2.7);
	sprintf(str, "%s.pdf(", fname);
	cMc->SaveAs(str);

	TCanvas *cMcR = new TCanvas("cMcR", "Monte Carlo", 800, 1000);
	cMcR->cd();
	hMcR->Fit("gaus", "", "", 1.5, 2.7);
	sprintf(str, "%s.pdf", fname);
	cMcR->SaveAs(str);
	
	TCanvas *cXY = new TCanvas("cXY", "DANSS XY", 800, 800);
	cXY->cd();
	hXY->Draw("colz");
	cXY->SaveAs(str);
	
	TCanvas *cExp = new TCanvas("cExp", "Danss", 800, 1000);
	cExp->cd();
	hExpC->Fit("gaus", "", "", 1.3, 2.7);
	cExp->SaveAs(str);
	
	TCanvas *cHits = new TCanvas("cHits", "Hits", 800, 1000);
	cHits->cd();
	hHitsC->Draw();
	hMcHits->Draw("same,hist");
	lg->Draw();
	cHits->Update();
	cHits->SaveAs(str);

	TCanvas *cE = new TCanvas("cE", "Hit energy", 800, 1000);
	cE->cd();
	hEC->Draw();
	hMcE->Draw("same,hist");
	lg->Draw();
	cE->Update();
	sprintf(str, "%s.pdf)", fname);
	cE->SaveAs(str);
}

void draw_22Na(void)
{
//	TFile *fMc = new TFile("/space/danss_root3/mcold/mc_22Na_simple.root");
//	TFile *fMc = new TFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/oldTransvProfile/mc_22Na_center_simple.root");
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_simple.root");
	TChain *tExpA = new TChain("DanssEvent");
//	tExpA->AddFile("/space/danss_root3/danss_data_002197_000.root");
//	tExpA->AddFile("/space/danss_root3/danss_data_002198_000.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_020243.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_020244.root");
	TChain *tExpB = new TChain("DanssEvent");
//	tExpB->AddFile("/space/danss_root3/danss_data_002193_000.root");
//	tExpB->AddFile("/space/danss_root3/danss_data_002194_000.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_020252.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_020253.root");
	draw_Src(tMc, tExpA, tExpB, 0.977, "^{22}Na", "22Na");
}

void draw_60Co(void)
{
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_simple.root");
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_020233.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_020234.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_020252.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_020253.root");
	draw_Src(tMc, tExpA, tExpB, 1.002, "^{60}Co", "60Co");
}

void draw_old22Na(void)
{
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_simple.root");
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002197.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002198.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002193.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002194.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002306.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002307.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002118.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002119.root");
	draw_Src(tMc, tExpA, tExpB, 3124.0/3103.0, "^{22}Na", "22Na_old");
}

void draw_old60Co(void)
{
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_simple.root");
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002106.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002107.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002102.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002103.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002118.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002119.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002306.root");
//	tExpB->AddFile("/mnt/space1/danss_root4/danss_002307.root");
	draw_Src(tMc, tExpA, tExpB, 2882.0/3024.0, "^{60}Co", "60Co_old");
}

void draw_oldDummy(void)
{
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_simple.root");
	TChain *tExpA = new TChain("DanssEvent");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002620.root");
	tExpA->AddFile("/mnt/space1/danss_root4/danss_002621.root");
	TChain *tExpB = new TChain("DanssEvent");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002620_000.root");
	tExpB->AddFile("/mnt/space1/danss_root4/danss_002621_000.root");
	draw_Src(tMc, tExpA, tExpB, 1.0, "No source", "2620_cmp");
}
