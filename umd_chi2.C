void umd_chi2(const char *fname, int last_bin = 24)
{
	const double R[3] = {10.7, 11.7, 12.7};
	const double ER[3] = {0, 0, 0};
	double Cnt[3], Err[3], chi2, tmp;
	int i;
	
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) return;
	TH1D *hUp = (TH1D*) f->Get("hUp_30");
	TH1D *hMid = (TH1D*) f->Get("hMid_30");
	TH1D *hDown = (TH1D*) f->Get("hDown_30");
	if (!(hUp && hMid && hDown)) return;
	
	Cnt[0] = hUp->IntegralAndError(1, last_bin, Err[0]);
	Cnt[1] = hMid->IntegralAndError(1, last_bin, Err[1]);
	Cnt[2] = hDown->IntegralAndError(1, last_bin, Err[2]);
	
	printf("Count rate 1-%2.0f MeV per day:\n", 0.25*last_bin + 0.75);
	printf("Up   = %7.2f +- %5.2f\n", Cnt[0], Err[0]);
	printf("Mid  = %7.2f +- %5.2f\n", Cnt[1], Err[1]);
	printf("Down = %7.2f +- %5.2f\n", Cnt[2], Err[2]);
	
	TF1 *fR2 = new TF1("fR2", "[0]/(x*x - 2.5*2.5/4)", 1, 100);
	TGraphErrors *gr = new TGraphErrors(3, R, Cnt, ER, Err);
	gr->Fit(fR2, "Q");
	chi2 = fR2->GetChisquare();
	printf("Chi^2 of 1/R^2 = %g\n", chi2);
	
	hUp->Scale(5000.0 / Cnt[0]);
	hMid->Scale(5000.0 / Cnt[1]);
	hDown->Scale(5000.0 / Cnt[2]);
	
	TH1D *hRes = (TH1D *) hUp->Clone("hMD2U_ratio");
	hRes->Add(hDown, hMid, 1, -1);
	hRes->Divide(hUp);
	TFile *ff = new TFile("umd_chi2.root", "RECREATE");
	ff->cd();
	hRes->Write();
	ff->Close();
	
	chi2 = 0;
	for (i=0; i<last_bin; i++) {
		tmp = hRes->GetBinContent(i+1) / hRes->GetBinError(i+1);
		chi2 += tmp * tmp;
	}
	printf("Chi^2 of (D - M)/U = %g\n", chi2);
	f->Close();
}
