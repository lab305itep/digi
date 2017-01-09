//	Draw MC monochrome positrons with and without annihilation
void draw_mc_posi(void)
{
	TFile *f[2][7];
	TH1D *h[4][7];
	TTree *t[2][7];
	TFile *fOut;
	TCanvas *cv;
	char stra[1024], strb[1024];
	int i, j;
	
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);
	
	for (i=0; i<7; i++) {
		sprintf(stra, "MC_PositronEnergyAllClusters_%1dMeV", i+1);
		sprintf(strb, "Monte Carlo reconstructed positron energy, all clusters, %1d MeV;MeV", i+1);
		h[0][i] = new TH1D(stra, strb, 100, 0, 10);
		sprintf(stra, "MC_PositronEnergyValidClusters_%1dMeV", i+1);
		sprintf(strb, "Monte Carlo reconstructed positron energy, valid clusters, %1d MeV;MeV", i+1);
		h[1][i] = new TH1D(stra, strb, 100, 0, 10);
		sprintf(stra, "MC_PositronEnergyNoAnnihilAllClusters_%1dMeV", i+1);
		sprintf(strb, "Monte Carlo reconstructed positron energy, no annihilation, all clusters, %1d MeV;MeV", i+1);
		h[2][i] = new TH1D(stra, strb, 100, 0, 10);
		sprintf(stra, "MC_PositronEnergyNoAnnihilValidClusters_%1dMeV", i+1);
		sprintf(strb, "Monte Carlo reconstructed positron energy, no annihilation, valid clusters, %1d MeV;MeV", i+1);
		h[3][i] = new TH1D(stra, strb, 100, 0, 10);
	}
	for (i=0; i<7; i++) {
		sprintf(stra, "danss_root3/mc_positron%1dMeV_simple.root", i+1);
		f[0][i] = new TFile(stra);
		t[0][i] = (TTree *) f[0][i]->Get("DanssEvent");
		if (!t[0][i]) break;
		sprintf(stra, "danss_root3/mc_positron_NoAnnihil_%1dMev_simple.root", i+1);
		f[1][i] = new TFile(stra);
		t[1][i] = (TTree *) f[1][i]->Get("DanssEvent");
		if (!t[1][i]) break;
	}
	if (i != 7) {
		printf("something wrong: i = %d (7)\n", i);
		return;
	}
	
	gROOT->cd();
	
	TCut cEdge("McX[0] > 4 && McX[0] < 96 && McX[1] > 4 && McX[1] < 96 && McX[2] > 4 && McX[2] < 96");
	TCut cValid("PositronValid > 0");
	
	for (j=0; j<2; j++) for (i=0; i < 7; i++) {
		t[j][i]->Project(h[2*j][i]->GetName(), "PositronEnergy", cEdge);
		t[j][i]->Project(h[2*j+1][i]->GetName(), "PositronEnergy", cEdge && cValid);
	}
	
	printf("Projected...\n");
	
	fOut = new TFile("MCMonochromePositrons.root", "RECREATE");
	fOut->cd();
	for (j=0; j<4; j++) for (i=0; i < 7; i++) h[j][i]->Write();
	fOut->Close();
	
	printf("File written...\n");

	for (j=0; j<2; j++) for (i=0; i < 7; i++) f[j][i]->Close();

	cv = new TCanvas("CV", "PositronEnergy", 1000, 800);
	for (i=0; i<7; i++) {
		if (i) cv->Clear();
		cv->Divide(2, 2);
		for (j=0; j<4; j++) {
			cv->cd(j+1);
			h[j][i]->Fit("gaus", "Q");
		}
		cv->SaveAs((!i) ? "mc_positron_energy.pdf(" : ((i == 6) ? "mc_positron_energy.pdf)" : "mc_positron_energy.pdf"));
	}
}
