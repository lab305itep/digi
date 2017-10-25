void make_fuel_hist(void)
{
	TFile *f;
	TFile *fout;
	TTree *t;
	TH1D  *h[2][4];
	TH1D  *hm[2][2];
	TH1D  *hr[2];
	const char fuel[4][10] = {"235U", "238U", "239Pu", "241Pu"};
	const float begin[4]   = {  0.69,   0.07,    0.21,    0.03};
	const float middle[4]  = {  0.58,   0.07,    0.30,    0.05};
	int i, j;
	char strs[128], strl[2048];
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut ct = cX && cY && cZ && cGamma;

	fout = new TFile("fuel_hist.root", "RECREATE");
	for (j=0; j<2; j++) for (i=0; i<4; i++) {
		sprintf(strs, "PositronEnergy_%s_%c", fuel[i], (j) ? 'D' : 'A');
		sprintf(strl, "Positron spectrum for %s%s;MeV", fuel[i], (j) ? ", with dead channels" : "");
		h[j][i] = new TH1D(strs, strl, 28, 1, 8);
		sprintf(strl, "danss_root3%s/mc_positron_%s_simple_newScale.root", (j) ? "/withdead-uncorr" : "", fuel[i]);
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
		t->Project(h[j][i]->GetName(), "PositronEnergy", ct);
		h[j][i]->Sumw2();
		fout->cd();
		h[j][i]->Write();
		f->Close();
	}
	
	hm[0][0] = new TH1D("PositronEnergy_Begin", "Positron spectrum at begin;MeV", 28, 1, 8);
	hm[0][1] = new TH1D("PositronEnergy_Middle", "Positron spectrum at middle;MeV", 28, 1, 8);
	hm[1][0] = new TH1D("PositronEnergy_Begin_D", "Positron spectrum at begin, with dead channels;MeV", 28, 1, 8);
	hm[1][1] = new TH1D("PositronEnergy_Middle_D", "Positron spectrum at middle, with dead channels;MeV", 28, 1, 8);
	
	for (j=0; j<2; j++) for (i=0; i<4; i++) {
		hm[j][0]->Add(h[j][i], begin[i]);
		hm[j][1]->Add(h[j][i], middle[i]);
	}

	fout->cd();
	for (j=0; j<2; j++) for (i=0; i<2; i++) hm[j][i]->Write();

	hr[0] = new TH1D("Ratio", "Positron spectra ratio begin to middle;MeV", 28, 1, 8);
	hr[1] = new TH1D("Ratio_D", "Positron spectra ratio begin to middle, with dead channels;MeV", 28, 1, 8);

	for (j=0; j<2; j++) hr[j]->Divide(hm[j][0], hm[j][1]);
	fout->cd();
	for (j=0; j<2; j++) hr[j]->Write();

	fout->Close();
}
