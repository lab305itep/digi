void draw_12C(void)
{
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	
	TFile f("/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mc_12C_deexitation_transcodeNew.root");
	TTree *t = (TTree *)f.Get("DanssEvent");
	if (!t) {
		printf("No DanssEvent tree\n");
		return;
	}
	TH1D *hPE = new TH1D("hPE", "Positron energy;MeV", 80, 0, 8);
	TH1D *hPEall = new TH1D("hPEall", "Positron energy, all events;MeV", 80, 0, 8);
	TH1D *hTE = new TH1D("hTE", "Total energy;MeV", 80, 0, 8);
	TH1D *hTEall = new TH1D("hTEall", "Total energy, all events;MeV", 80, 0, 8);
	TH1D *hSiPM = new TH1D("hSiPM", "SiPM energy;MeV", 66, 0, 8);
	TH1D *hSiPMall = new TH1D("hSiPMall", "SiPM energy, all events;MeV", 66, 0, 8);
	TH1D *hPMT = new TH1D("hPMT", "PMT energy;MeV", 80, 0, 8);
	TH1D *hPMTall = new TH1D("hPMTall", "PMT energy, all events;MeV", 80, 0, 8);
	
	TCut cX("(PositronX[0]>10 && PositronX[0]<90) || PositronX[0] < 0");
	TCut cY("(PositronX[1]>10 && PositronX[1]<90) || PositronX[1] < 0");
	TCut cZ("PositronX[2]>10 && PositronX[2]<90");
	TCut cXYZ = cX && cY && cZ;
	TCut cN("AnnihilationGammas == 1 || AnnihilationGammas == 2");
	
	t->Project("hPE", "(PositronEnergy-0.179)/0.929", cN && cXYZ);
	t->Project("hPEall", "(PositronEnergy-0.179)/0.929", cXYZ);
	t->Project("hTE", "(SiPmCleanEnergy + PmtCleanEnergy)/2", cN && cXYZ);
	t->Project("hTEall", "(SiPmCleanEnergy + PmtCleanEnergy)/2", cXYZ);
	t->Project("hSiPM", "SiPmCleanEnergy", cN && cXYZ);
	t->Project("hSiPMall", "SiPmCleanEnergy", cXYZ);
	t->Project("hPMT", "PmtCleanEnergy", cN && cXYZ);
	t->Project("hPMTall", "PmtCleanEnergy", cXYZ);

	TCanvas *cE = new TCanvas("cE", "Energy", 1200, 800);
	cE->SaveAs("12C.pdf[");
	cE->Divide(2, 2);
	cE->cd(1);
	hPE->Fit("gaus", "", "", 3.5, 6);
	cE->cd(2);
	hTE->Fit("gaus", "", "", 3.5, 6);
	cE->cd(3);
	hSiPM->Fit("gaus", "", "", 3.5, 6);
	cE->cd(4);
	hPMT->Fit("gaus", "", "", 3.5, 6);

	cE->SaveAs("12C.pdf");
	cE->Clear();
	cE->Divide(2, 2);
	cE->cd(1);
	hPEall->Fit("gaus", "", "", 3.5, 6);
	cE->cd(2);
	hTEall->Fit("gaus", "", "", 3.5, 6);
	cE->cd(3);
	hSiPMall->Fit("gaus", "", "", 3.5, 6);
	cE->cd(4);
	hPMTall->Fit("gaus", "", "", 3.5, 6);
	cE->SaveAs("12C.pdf");
	cE->SaveAs("12C.pdf]");
}

