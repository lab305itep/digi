{
	const double E_pos[] =  {1.104, 2.032, 2.961, 3.982, 4.802, 5.724, 6.669};
	const double EE_pos[] = {0.004, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011};
	double E0[7];
	double EE0[7];
	double E_posd[7];
	double a, b;
	TH1D *hE = new TH1D("hE", "Positron energy seen;MeV;MeV", 7, 0.5, 7.5);
	TH1D *hEd = new TH1D("hEd", "Deviation;MeV;MeV", 7, 0.5, 7.5);
	int i;
	TGraphErrors *gr;
	TGraphErrors *grd;
	TF1 *f = new TF1("F", "pol1");

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	hE->SetMinimum(0);
	hE->SetMaximum(8);
	hEd->SetMinimum(-0.5);
	hEd->SetMaximum(0.5);
	for (i=0; i<7; i++) {
		E0[i] = i + 1;
		EE0[i] = 0;
	}
	gr = new TGraphErrors(7, E0, E_pos, EE0, EE_pos);
	gr->SetMarkerStyle(21);
	gr->Fit("F");
	a = f->GetParameter(0);
	b = f->GetParameter(1);
	for (i=0; i<7; i++) E_posd[i] = E_pos[i] - a - b * E0[i];
	grd = new TGraphErrors(7, E0, E_posd, EE0, EE_pos);
	grd->SetMarkerStyle(21);

	TCanvas *cv = new TCanvas();
	cv->Divide(1, 2);
	cv->cd(1);
	hE->Draw();
	gr->Draw("p");
	cv->cd(2);
	hEd->Draw();
	grd->Draw("p");	
}

