void reactor_anomaly(void)
{
	const double reno[20] = {
		1.00, 1.02, 0.995, 1.00, 1.00, 1.015, 1.005, 0.98,  1.02, 1.055,
		1.09, 1.12, 1.13,  1.12, 1.06, 1.05,  1.04,  1.005, 0.95, 0.98
	};
	TF1 *fGaus;
	TH1D *hReno;
	TH1D *hDanss1;
	TH1D *hDanss2;
	int i;
	
	gStyle->SetOptStat(0);
	hReno = new TH1D("hReno", "RENO data on the ratio experiment to Monte-Carlo;MeV", 20, 1, 6);
	hDanss1 = new TH1D("hDanss1", "RENO ratio spoiled by DANSS resolution 30%/#sqrt{E};MeV", 20, 1, 6);
	hDanss2 = new TH1D("hDanss2", "RENO ratio spoiled by DANSS resolution 45%/#sqrt{E};MeV", 20, 1, 6);
	for (i=0; i<20; i++) {
		hDanss1->SetBinContent(i+1, 1);
		hDanss2->SetBinContent(i+1, 1);
	}

	hReno->SetMarkerColor(kBlack);
	hReno->SetMarkerStyle(20);
	hReno->SetMarkerSize(2);
	hDanss1->SetMarkerColor(kBlue);
	hDanss1->SetMarkerStyle(21);
	hDanss1->SetMarkerSize(2);
	hDanss2->SetMarkerColor(kRed);
	hDanss2->SetMarkerStyle(22);
	hDanss2->SetMarkerSize(2);
	
	for(i=0; i<20; i++) hReno->SetBinContent(i+1, reno[i]);
	
	fGaus = new TF1("fGaus", "gaus", -20.0, 20.0);
	for (i=0; i<20; i++) {
		fGaus->SetParameter(0, 1);	// amplitude
		fGaus->SetParameter(1, 1.125 + 0.25*i);	// position
		fGaus->SetParameter(2, 0.3*sqrt(1.125 + 0.25*i));	// DANSS resolution: 30%/sqrt(E)
		hDanss1->Add(fGaus, (reno[i]-1)*0.25/fGaus->Integral(-20.0, 20.0));
		fGaus->SetParameter(2, 0.45*sqrt(1.125 + 0.25*i));	// DANSS resolution: 30%/sqrt(E)
		hDanss2->Add(fGaus, (reno[i]-1)*0.25/fGaus->Integral(-20.0, 20.0));
	}
//	hReno->SetMinimum(0.8);

	hReno->DrawCopy("P0");
	hDanss1->DrawCopy("P0,same");
	hDanss2->DrawCopy("P0,same");

	TLegend *lg = new TLegend(0.6, 0.12, 0.8, 0.25);
	lg->AddEntry(hReno, "RENO data", "P");
	lg->AddEntry(hDanss1, "30%/#sqrt{E}", "P");
	lg->AddEntry(hDanss2, "45%/#sqrt{E}", "P");
	lg->Draw();
}
