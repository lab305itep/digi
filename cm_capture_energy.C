void cm_capture_energy(const char *fname)
{
	int i, j;
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetLineWidth(2);
	gStyle->SetPalette(kRainBow, 0);

	TH1D *hcm = new TH1D("HDC", ";Energy of delayed event, MeV", 60, 0, 12);
	hcm->SetLineWidth(4);
	hcm->SetMarkerStyle(21);
	hcm->SetLineColor(kBlue);
	hcm->SetFillColor(kBlue-10);
	hcm->SetMarkerColor(kBlue);
	hcm->GetYaxis()->SetLabelSize(0.06);

	TH2D *hcmxy = new TH2D("HXY", ";X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	hcmxy->GetXaxis()->SetLabelSize(0.05);
	hcmxy->GetYaxis()->SetLabelSize(0.05);
	hcmxy->GetZaxis()->SetLabelSize(0.05);
//	for (i=0; i<25; i++) for (j=0; j<25; j++) hcmxy->Fill(4.0*i+2, 4.0*j+2);

	TFile f(fname);
	if (!f.IsOpen()) return;
	TTree *t = (TTree *) f.Get("DanssCm");
	if (!t) return;
	gROOT->cd();
	t->Project("HDC", "(SiPmCleanEnergy[1]+PmtCleanEnergy[1])/2", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2");
	t->Project("HXY", "NeutronX[1][1]+2:NeutronX[1][0]+2", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2 && NeutronX[1][1]>=0 && NeutronX[1][0]>=0");
	hcmxy->Fill(18., 2., -1);

	TCanvas *c1 = new TCanvas("CV", "Neutron Capture", 800, 800);
	hcm->Draw("hist,e");
	TCanvas *c2 = new TCanvas("CXY", "Neutron XY", 800, 800);
	TVirtualPad *pd2 = c2->cd(0);
	pd2->SetRightMargin(0.16);
	hcmxy->Draw("COLZ");

	f.Close();
}

