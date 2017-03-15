void make_mono_positrons(void)
{
	const char mfiles[10][128] = {
		"/space/danss_root3/mc_positron1-5MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron2MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron3MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron4MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron5MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron6MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron7MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron8MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron9MeV_simple_newScale.root",
		"/space/danss_root3/mc_positron10MeV_simple_newScale.root"
	};
	const double pe[10] = {1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	int i;
	char str[1024];
	char strs[128];
	TH1D *h[10];
	TFile *f[10];
	TTree *t[10];
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TFile *fOut = new TFile("mc_mono_positrons.root", "RECREATE");
	for (i=0; i<10; i++) {
		sprintf(str, "Reconstructed positron energy from %5.1f MeV;MeV;", pe[i]);
		sprintf(strs, "PositronEnergy_%5.1f", pe[i]);
		h[i] = new TH1D(strs, str, 75, 0, 15);
	}
	for (i=0; i<10; i++) {
		f[i] = new TFile(mfiles[i]);
		if (!f[i]->IsOpen()) return;
		t[i] = (TTree *) f[i]->Get("DanssEvent");
		if (!t[i]) {
			printf("Something is wrong with file %s\n", mfiles[i]);
			return;
		}
	}
	fOut->cd();
	for (i=0; i<10; i++) t[i]->Project(h[i]->GetName(), "(PositronEnergy-0.179)/0.929", cX && cY && cZ);

	for (i=0; i<10; i++) h[i]->Write();
	for (i=0; i<10; i++) f[i]->Close();
	fOut->Close();
}
