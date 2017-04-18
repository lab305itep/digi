void make_fuel_hist(void)
{
	TFile *f;
	TFile *fout;
	TTree *t;
	TH1D  *h[4];
	TH1D  *hm[2];
	TH1D  *hr;
	const char fuel[4][10] = {"235U", "238U", "239Pu", "241Pu"};
	const float begin[4]   = {  0.69,   0.07,    0.21,    0.03};
	const float middle[4]  = {  0.58,   0.07,    0.30,    0.05};
	int i;
	char strs[128], strl[2048];
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut ct = cX && cY && cZ && cGamma;

	fout = new TFile("fuel_hist.root", "RECREATE");
	for (i=0; i<4; i++) {
		sprintf(strs, "PositronEnergy_%s", fuel[i]);
		sprintf(strl, "Positron spectrum for %s;MeV", fuel[i]);
		h[i] = new TH1D(strs, strl, 28, 1, 8);
		sprintf(strl, "danss_root3/mc_positron_%s_simple_newScale.root", fuel[i]);
		f = new TFile(strl);
		if (!f->IsOpen()) {
			printf("File not found %s\n", strl);
			continue;
		}
		t = (TTree *) f->Get("DanssEvent");
		if (!t) {
			printf("Tree DanssEvent not found %s\n", strl);
			continue;
		}
		fout->cd();
		t->Project(h[i]->GetName(), "PositronEnergy", ct);
		h[i]->Sumw2();
		fout->cd();
		h[i]->Write();
		f->Close();
	}
	
	hm[0] = new TH1D("PositronEnergy_Begin", "Positron spectrum at begin;MeV", 28, 1, 8);
	hm[1] = new TH1D("PositronEnergy_Middle", "Positron spectrum at middle;MeV", 28, 1, 8);
	
	for (i=0; i<4; i++) {
		hm[0]->Add(h[i], begin[i]);
		hm[1]->Add(h[i], middle[i]);
	}

	fout->cd();
	hm[0]->Write();
	hm[1]->Write();

	hr = new TH1D("Ratio", "Positron spectra ratio begin to middle;MeV", 28, 1, 8);

	hr->Divide(hm[0], hm[1]);
	fout->cd();
	hr->Write();

	fout->Close();
}
