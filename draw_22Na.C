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
//	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400");
		r2 = 0;
		for (j=0; j<3; j++) r2 += (EVT.NeutronX[0] - 48) * (EVT.NeutronX[0] - 48);
		if (r2 >= 400) continue;
//	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
		if (EVT.VetoCleanHits >= 2 || EVT.VetoCleanEnergy >= 4) continue;
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
//	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400");
		r2 = 0;
		for (j=0; j<3; j++) r2 += (EVT.NeutronX[0] - 48) * (EVT.NeutronX[0] - 48);
		if (r2 >= 400) continue;
//	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
		if (EVT.VetoCleanHits >= 2 || EVT.VetoCleanEnergy >= 4) continue;
		Esum = 0;
		for (j=0; j<EVT.NHits; j++) if (HitArray.type[j].type == 0 || HitArray.type[j].type == 1) Esum += HitArray.E[j] * random->Gaus(1.0, 0.26);
		hist->Fill(Esum/2);
	}
}

void draw_Src(TChain *tMc, TChain *tExpA, TChain *tExpB, double rAB, const char *name, const char *fname)
{
	char str[256];
	
//	gStyle->SetOptStat("i");
//	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
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
	sprintf(str, "SiPM energy deposit in %s decay;E, MeV", name);
	TH1D *hExpSiPMA = new TH1D("hExpSiPMA", str, 50, 0, 5);
	TH1D *hExpSiPMB = new TH1D("hExpSiPMB", str, 50, 0, 5);
	TH1D *hExpSiPMC = new TH1D("hExpSiPMC", str, 50, 0, 5);
	sprintf(str, "PMT energy deposit in %s decay;E, MeV", name);
	TH1D *hExpPMTA = new TH1D("hExpPMTA", str, 50, 0, 5);
	TH1D *hExpPMTB = new TH1D("hExpPMTB", str, 50, 0, 5);
	TH1D *hExpPMTC = new TH1D("hExpPMTC", str, 50, 0, 5);
	sprintf(str, "Number of hits from %s decay", name);
	TH1D *hHitsA = new TH1D("hHitsA", str, 20, 0, 20);
	TH1D *hHitsB = new TH1D("hHitsB", str, 20, 0, 20);
	TH1D *hHitsC = new TH1D("hHitsC", str, 20, 0, 20);
	TH1D *hEA = new TH1D("hEA", "Hit energy;E, MeV", 60, 0, 3);
	TH1D *hEB = new TH1D("hEB", "Hit energy;E, MeV", 60, 0, 3);
	TH1D *hEC = new TH1D("hEC", "Hit energy;E, MeV", 60, 0, 3);

	TCut cxyz("NeutronX[0] >= 0 && NeutronX[1] >= 0 && NeutronX[2] >= 0");
	TCut cz50("(NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 100");
	TCut ccc("(NeutronX[0] - 48) * (NeutronX[0] - 48) + (NeutronX[1] - 48) * (NeutronX[1] - 48) + (NeutronX[2] - 49.5) * (NeutronX[2] - 49.5) < 400");
	TCut cVeto("VetoCleanHits < 2 && VetoCleanEnergy < 4");
	
	tMc->Project("hMc", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc && cVeto);
	tMc->Project("hMcHits", "SiPmCleanHits", cxyz && ccc && cVeto);
	project_hits_distrib(tMc, hMcE);
	tExpA->Project("hXY", "NeutronX[1]+2:NeutronX[0]+2", cxyz && cz50 && cVeto);
	tExpA->Project("hExpA", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc && cVeto);
	tExpB->Project("hExpB", "(SiPmCleanEnergy+PmtCleanEnergy)/2", cxyz && ccc && cVeto);
	tExpA->Project("hExpSiPMA", "SiPmCleanEnergy", cxyz && ccc && cVeto);
	tExpB->Project("hExpSiPMB", "SiPmCleanEnergy", cxyz && ccc && cVeto);
	tExpA->Project("hExpPMTA", "PmtCleanEnergy", cxyz && ccc && cVeto);
	tExpB->Project("hExpPMTB", "PmtCleanEnergy", cxyz && ccc && cVeto);
	tExpA->Project("hHitsA", "SiPmCleanHits", cxyz && ccc && cVeto);
	tExpB->Project("hHitsB", "SiPmCleanHits", cxyz && ccc && cVeto);
	project_hits_distrib(tExpA, hEA);
	project_hits_distrib(tExpB, hEB);
	project_mc_random_energy(tMc, hMcR);
	
	hMc->Sumw2();
	hMcR->Sumw2();
	hMcHits->Sumw2();
	hMcE->Sumw2();
	hExpA->Sumw2();
	hExpB->Sumw2();
	hExpSiPMA->Sumw2();
	hExpSiPMB->Sumw2();
	hExpPMTA->Sumw2();
	hExpPMTB->Sumw2();
	hHitsA->Sumw2();
	hHitsB->Sumw2();
	hEA->Sumw2();
	hEB->Sumw2();
	
	hExpC->Add(hExpA, hExpB, 1.0, -rAB);
	hExpSiPMC->Add(hExpSiPMA, hExpSiPMB, 1.0, -rAB);
	hExpPMTC->Add(hExpPMTA, hExpPMTB, 1.0, -rAB);
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
	hExpSiPMC->GetYaxis()->SetLabelSize(0.05);
	hExpPMTC->GetYaxis()->SetLabelSize(0.05);
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
	hExpC->SetLineWidth(2);
	hExpC->Draw();
//	hExpC->Fit("gaus", "", "", 1.3, 2.7);
	cExp->SaveAs(str);

	TCanvas *cExpSiPM = new TCanvas("cExpSiPM", "Danss SiPM", 800, 1000);
	cExpSiPM->cd();
	hExpSiPMC->Fit("gaus", "", "", 1.3, 2.7);
	cExpSiPM->SaveAs(str);

	TCanvas *cExpPMT = new TCanvas("cExpPMT", "Danss PMT", 800, 1000);
	cExpPMT->cd();
	hExpPMTC->Fit("gaus", "", "", 1.3, 2.7);
	cExpPMT->SaveAs(str);
	
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

void draw_22Na(int iFull, int iSer)
{
	char str[128];
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_22Na_center_transcodeNew.root");
	TChain *tExpA = new TChain("DanssEvent");
	TChain *tExpB = new TChain("DanssEvent");
	switch (2*iSer + iFull) {
	case 0:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002197.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002198.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002285.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002286.root");
		break;
	case 1:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002197_000.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002198_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002285_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002286_000.root");
		break;
	case 2:
		tExpA->AddFile("/mnt/dcopy0/danss_root4/danss_012380.root");
		tExpA->AddFile("/mnt/dcopy0/danss_root4/danss_012381.root");
		tExpB->AddFile("/mnt/dcopy0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/dcopy0/danss_root4/danss_012302.root");
		break;
	case 4:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020243.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020244.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020252.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020253.root");
		break;
	case 5:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020243_000.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020244_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020252_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020253_000.root");
		break;
	default:
		printf("Not yet.\n");
		return;
	}
	sprintf(str, "22Na_%d_%s_20cm", iSer, (iFull) ? "full" : "fast");
	draw_Src(tMc, tExpA, tExpB, 0.977, "^{22}Na", str);
}

void draw_60Co(int iFull, int iSer)
{
	char str[128];
	TChain *tMc = new TChain("DanssEvent");
	tMc->AddFile("/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_60Co_center_transcodeNew.root");
	TChain *tExpA = new TChain("DanssEvent");
	TChain *tExpB = new TChain("DanssEvent");
	switch (2*iSer + iFull) {
	case 0:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002106.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002107.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002285.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002286.root");
		break;
	case 1:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002106_000.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_002107_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002285_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_002286_000.root");
		break;
	case 2:
		tExpA->AddFile("/mnt/dcopy0/danss_root4/danss_012310.root");
		tExpA->AddFile("/mnt/dcopy0/danss_root4/danss_012311.root");
		tExpB->AddFile("/mnt/dcopy0/danss_root4/danss_012301.root");
		tExpB->AddFile("/mnt/dcopy0/danss_root4/danss_012302.root");
		break;
	case 4:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020233.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020234.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020252.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020253.root");
		break;
		
	case 5:
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020233_000.root");
		tExpA->AddFile("/mnt/space1/danss_root4/danss_020234_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020252_000.root");
		tExpB->AddFile("/mnt/space1/danss_root4/danss_020253_000.root");
		break;
	default:
		printf("Not yet.\n");
		return;
	}
	sprintf(str, "60Co_%d_%s_20cm", iSer, (iFull) ? "full" : "fast");
	draw_Src(tMc, tExpA, tExpB, 0.977, "^{60}Co", str);
}

void draw_all(void)
{
	int i;
	for (i=0; i<6; i++) {
		draw_22Na(i&1, i/2);
		draw_60Co(i&1, i/2);
	}
}
