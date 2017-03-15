void energy_resolution_4MeV(const char *fname)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetLineWidth(4);

	TH1D *h4 = new TH1D("H4", ";Cluster energy, MeV", 40, 0, 8);
	h4->SetMarkerStyle(21);
	h4->SetLineColor(kBlue);
	h4->SetFillColor(kBlue-10);
	h4->SetMarkerColor(kBlue);
	h4->GetYaxis()->SetLabelSize(0.06);

	TF1 *fG = new TF1("fG", "gaus", 2.5, 5.5);
	fG->SetLineWidth(4);
	fG->SetLineColor(kRed);

	TFile f(fname);
	if (!f.IsOpen()) return;
	TTree *t = (TTree *) f.Get("DanssEvent");
	if (!t) return;
	gROOT->cd();
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");

	t->Project("H4", "(PositronEnergy-0.179)/0.929", cX && cY && cZ);

	TCanvas *c1 = new TCanvas("CV", "Neutron Capture", 800, 800);
	h4->Fit(fG, "", "hist,e", 2.7, 5.5);
	fG->Draw("same");
	c1->Write("mc_energy_resolution_4MeV.png");

	f.Close();
}

