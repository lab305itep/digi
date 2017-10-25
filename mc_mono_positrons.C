//	Draw MC monochrome positrons
TH1D *mc_mono_positrons(const char *fname, const char *hname, const char *htitle)
{
	TFile *f;
	TH1D *h;
	TTree *t;
	
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);
	
	h = new TH1D(hname, htitle, 48, 0, 12);
	f = new TFile(fname);
	t = (TTree *) f->Get("DanssEvent");
	if (!t) {
		printf("No tree DanssEvent in file %s\n", fname);
		delete h;
		return NULL;
	}
//		Set cuts
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
	TCut cSel = cX && cY && cZ && cGamma && cGammaMax;

	gROOT->cd();
	t->Project(hname, "(PositronEnergy-0.179)/0.929", cSel);
	return h;
}

void all_mono_positrons(void)
{
	const double Energy[12] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	const char *fvar[12] = {"0-5", "1", "1-5", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
	const char *dir = "/mnt/space1/danss_root4/LY_siPm18_pmt20_new/newTransvProfile";
	const char *tvar[2] = {"simple", "transcodeNew"};
	char strF[1024];
	char strH[1024];
	char strT[1024];
	int i, j;
	TH1D *h[12][2];
	double elow, ehigh;
	
	for (i=0; i<12; i++) for (j=0; j<2; j++) {
		sprintf(strF, "%s/mc_positron%sMeV_%s.root", dir, fvar[i], tvar[j]);
		sprintf(strH, "hPE%s_%s", fvar[i], tvar[j]);
		sprintf(strT, "Reconstructed positron energy for MC mono positrons at %4.1f MeV;E, MeV", Energy[i]);
		h[i][j] = mc_mono_positrons(strF, strH, strT);
	}
	TFile *fOut = new TFile("mc_mono_positrons_v2.root", "RECREATE");
	fOut->cd();
	for (i=0; i<12; i++) for (j=0; j<2; j++) h[i][j]->Write();
	
	for (i=0; i<12; i++) {
		h[i][0]->SetLineColor(kRed);
		h[i][1]->SetLineColor(kBlue);
	}
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	TCanvas *cv = new TCanvas("CV", "CV", 600, 800);
	
	cv->SaveAs("mc_mono_positrons_v2.pdf[");
	
	for (i=0; i<12; i++) {
		cv->Clear();
		elow=Energy[i] - 2 * h[i][0]->GetRMS();
		ehigh=Energy[i] + 2 * h[i][0]->GetRMS();
		h[i][0]->Fit("gaus", "", "", elow, ehigh);
		h[i][1]->Draw("same");
		cv->SaveAs("mc_mono_positrons_v2.pdf");
	}
	
	cv->SaveAs("mc_mono_positrons_v2.pdf]");
	
	delete cv;
	fOut->Close();
}

void all_mono_positrons_250keV(void)
{
	int Energy;
	char fvar[20];
	const char *dir = "/mnt/root0/danss_root4/LY_siPm18_pmt20_new/newTransvProfile/mono_positrons_250keV/";
	char strF[1024];
	char strH[1024];
	char strT[1024];
	int i;
	TH1D *h[48];
	double elow, ehigh;
	
	for (i = 0; i < 48; i++) {
		Energy = 125 + 250 * i;
		sprintf(fvar, "%d-%d", Energy/1000, Energy%1000);
		sprintf(strF, "%s/mc_positron%sMeV_transcodeNew.root", dir, fvar);
		sprintf(strH, "hPE%s", fvar);
		sprintf(strT, "Reconstructed positron energy for MC mono positrons at %6.3f MeV;E, MeV", Energy/1000.0);
		h[i] = mc_mono_positrons(strF, strH, strT);
	}
	TFile *fOut = new TFile("mc_mono_positrons_v3.root", "RECREATE");
	fOut->cd();
	for (i=0; i<48; i++) h[i]->Write();
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	TCanvas *cv = new TCanvas("CV", "CV", 600, 800);
	
	cv->SaveAs("mc_mono_positrons_v3.pdf[");
	
	for (i=0; i<48; i++) {
		cv->Clear();
		Energy = 0.125 + 0.25 * i;
		elow=Energy - 2 * h[i]->GetRMS();
		ehigh=Energy + 2 * h[i]->GetRMS();
		h[i]->Fit("gaus", "", "", elow, ehigh);
		cv->SaveAs("mc_mono_positrons_v3.pdf");
	}
	
	cv->SaveAs("mc_mono_positrons_v3.pdf]");
	
	delete cv;
	fOut->Close();
}

