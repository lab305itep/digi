//	Crystall Ball function
double CBfunction(double *x, double *par)
{
	double F;
	
	double X = x[0];
	double N = par[0];
	double alpha = par[1];
	double n = par[2];
	double X0 = par[3];
	double sigma = par[4];
	
	double dx = (X - X0) / sigma;
	double A = exp(n*log(n/fabs(alpha)) - alpha*alpha/2);
	double B = n / fabs(alpha) - fabs(alpha);
	
	F = (dx > -alpha) ? exp(-dx*dx/2) : A * exp(-n * log(B - dx));
	return N*F;
}

//	Gd mixture:
//	155Gd	8.536 MeV	18.78%
//	157Gd	7.937MeV	81.21%
double GdFunction(double *x, double *par)
{
	const double frac = 0.1878;
	const double E155 = 8.536;
	const double E157 = 7.937;

	double parA[5];
	double parB[5];
	
	memcpy(parA, par, sizeof(parA));
	memcpy(parB, par, sizeof(parB));
	
	parA[0] = frac * par[0];
	parB[0] = (1 - frac) * par[0];
	parA[3] = par[0] + E155 - E157;
	
	return CBfunction(x, parA) + CBfunction(x, parB);
}

//	H + Gd
double CaptureFunction(double *x, double *par)
{
	double parGd[5];
	double parH[5];
	
	memcpy(parGd, par, sizeof(parGd));
	memcpy(parH, &par[5], sizeof(parH));
	
	return GdFunction(x, parGd) + CBfunction(x, parH);
}

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

	TH1D *hcm = new TH1D("HDC", ";Energy of delayed event, MeV", 120, 0, 12);
	hcm->SetLineWidth(4);
	hcm->SetMarkerStyle(21);
	hcm->SetLineColor(kBlue);
	hcm->SetFillColor(kBlue-10);
	hcm->SetMarkerColor(kBlue);
	hcm->GetYaxis()->SetLabelSize(0.06);

	TH1D *hpH = new TH1D("HPH", ";Energy of cluster in delayed event, MeV", 50, 0, 5);
	hpH->SetLineWidth(4);
	hpH->SetMarkerStyle(21);
	hpH->SetLineColor(kBlue);
	hpH->SetFillColor(kBlue-10);
	hpH->SetMarkerColor(kBlue);
	hpH->GetYaxis()->SetLabelSize(0.06);

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
	t->Project("HPH", "PositronEnergy[1]", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2 && Hits[1] < 5");
	t->Project("HXY", "NeutronX[1][1]+2:NeutronX[1][0]+2", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2 && NeutronX[1][1]>=0 && NeutronX[1][0]>=0");
	hcmxy->Fill(18., 2., -1);
	TF1 *fCapt = new TF1("fCapt", CaptureFunction, 0, 20, 10);
	fCapt->SetParNames("Gd:Const", "Gd:#alpha", "Gd:n", "Gd:mean", "Gd:#sigma", "H:Const", "H:#alpha", "H:n", "H:mean", "H:#sigma");
	fCapt->SetParameters(hcm->GetMaximum(), 1, 1, 6.7, 1, hcm->GetMaximum()/5, 1, 1, 1.7, 1);
	fCapt->SetLineColor(kRed);
	fCapt->SetLineWidth(2);
	
	TF1 *fH = new TF1("fH", CBfunction, 0, 20, 5);
	fH->SetParNames("Const", "#alpha", "n", "mean", "#sigma");
	fH->SetParameters(hpH->GetMaximum(), 1, 1, 1.7, 1);
	fH->SetLineColor(kRed);
	fH->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("CV", "Neutron Capture", 800, 800);
	hcm->Fit(fCapt, "", "", 0.7, 12.);
	TCanvas *c2 = new TCanvas("CXY", "Neutron XY", 800, 800);
	TVirtualPad *pd2 = c2->cd(0);
	pd2->SetRightMargin(0.16);
	hcmxy->Draw("COLZ");
	TCanvas *c3 = new TCanvas("CPH", "Cluster", 800, 800);
	hpH->Fit(fH, "", "", 1, 4);

	f.Close();
}
