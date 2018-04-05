TRandom2 rnd;

class MyRandom {
    public:
	inline MyRandom(void) {;};
	inline ~MyRandom(void) {;};
	static inline double Gaus(double mean = 0, double sigma = 1) 
		{
		return rnd.Gaus(mean, sigma);
	};
	static inline double GausAdd(double val, double sigma)
	{
		return rnd.Gaus(val, sqrt(val)*sigma);
	};
	static inline double GausAdd2(double val, double sigma)
	{
		return rnd.Gaus(val, val*sigma);
	};
};

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
	parA[3] = par[3] + E155 - E157;
	
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

//	H + pol2
double MCFunction(double *x, double *par)
{
	return CBfunction(x, par) + par[5] + x[0]*par[6] + x[0]*x[0]*par[7];
}

void cm_capture_energy(const char *fname, const char *mcname, double kRndm)
{
	int i, j;
	char str[2048];
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.055);
	gStyle->SetTitleYSize(0.055);
	gStyle->SetLabelSize(0.055);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetLineWidth(2);
	gStyle->SetPalette(kRainBow, 0);

	TH1D *hcm = new TH1D("HDC", "^{248}Cm source data;Energy of delayed event, MeV", 120, 0, 12);
	hcm->SetLineWidth(4);
	hcm->SetMarkerStyle(21);
	hcm->SetLineColor(kBlue);
	hcm->SetFillColor(kBlue-10);
	hcm->SetMarkerColor(kBlue);
	hcm->GetYaxis()->SetLabelSize(0.055);

	TH1D *hcmSi = new TH1D("HDCSi", "^{248}Cm source data: SiPM;Energy of delayed event, MeV", 120, 0, 12);
	hcmSi->SetLineWidth(4);
	hcmSi->SetMarkerStyle(21);
	hcmSi->SetLineColor(kBlue);
	hcmSi->SetFillColor(kBlue-10);
	hcmSi->SetMarkerColor(kBlue);
	hcmSi->GetYaxis()->SetLabelSize(0.055);

	TH1D *hcmPMT = new TH1D("HDCPMT", "^{248}Cm source data: PMT;Energy of delayed event, MeV", 120, 0, 12);
	hcmPMT->SetLineWidth(4);
	hcmPMT->SetMarkerStyle(21);
	hcmPMT->SetLineColor(kBlue);
	hcmPMT->SetFillColor(kBlue-10);
	hcmPMT->SetMarkerColor(kBlue);
	hcmPMT->GetYaxis()->SetLabelSize(0.055);

	TH1D *hMc = new TH1D("HMC", "Neutron Monte Carlo;Energy of delayed event, MeV", 120, 0, 12);
	hMc->SetLineWidth(4);
	hMc->SetMarkerStyle(21);
	hMc->SetLineColor(kBlue);
	hMc->SetFillColor(kBlue-10);
	hMc->SetMarkerColor(kBlue);
	hMc->GetYaxis()->SetLabelSize(0.06);

	TH1D *hMcSi = new TH1D("HMCSi", "Neutron Monte Carlo: SiPM;Energy of delayed event, MeV", 60, 0, 12);
	hMcSi->SetLineWidth(4);
	hMcSi->SetMarkerStyle(21);
	hMcSi->SetLineColor(kBlue);
	hMcSi->SetFillColor(kBlue-10);
	hMcSi->SetMarkerColor(kBlue);
	hMcSi->GetYaxis()->SetLabelSize(0.06);

	TH1D *hMcPMT = new TH1D("HMCPMT", "Neutron Monte Carlo: PMT;Energy of delayed event, MeV", 120, 0, 12);
	hMcPMT->SetLineWidth(4);
	hMcPMT->SetMarkerStyle(21);
	hMcPMT->SetLineColor(kBlue);
	hMcPMT->SetFillColor(kBlue-10);
	hMcPMT->SetMarkerColor(kBlue);
	hMcPMT->GetYaxis()->SetLabelSize(0.06);

	TH1D *hpH = new TH1D("HPH", "^{248}Cm source data;Energy of cluster in delayed event, MeV", 50, 0, 5);
	hpH->SetLineWidth(4);
	hpH->SetMarkerStyle(21);
	hpH->SetLineColor(kBlue);
	hpH->SetFillColor(kBlue-10);
	hpH->SetMarkerColor(kBlue);
	hpH->GetYaxis()->SetLabelSize(0.06);

	TH2D *hcmxy = new TH2D("HXY", "^{248}Cm source data;X, cm;Y, cm", 25, 0, 100, 25, 0, 100);
	hcmxy->GetXaxis()->SetLabelSize(0.05);
	hcmxy->GetYaxis()->SetLabelSize(0.05);
	hcmxy->GetZaxis()->SetLabelSize(0.05);
//	for (i=0; i<25; i++) for (j=0; j<25; j++) hcmxy->Fill(4.0*i+2, 4.0*j+2);

	TFile f(fname);
	if (!f.IsOpen()) return;
	TTree *t = (TTree *) f.Get("DanssCm");
	if (!t) return;
	TFile fMc(mcname);
	if (!fMc.IsOpen()) return;
	TTree *tMc = (TTree *) fMc.Get("DanssEvent");
	if (!tMc) return;
	
	gROOT->cd();
	t->Project("HDC", "(SiPmCleanEnergy[1]+PmtCleanEnergy[1])/2", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2");
	sprintf(str, "MyRandom::GausAdd((SiPmCleanEnergy+PmtCleanEnergy)/2, %6.4f)", kRndm);
	tMc->Project("HMC", str);
	t->Project("HDCSi", "SiPmCleanEnergy[1]", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2");
	sprintf(str, "MyRandom::GausAdd(SiPmCleanEnergy, %6.4f)", kRndm);
	tMc->Project("HMCSi", "SiPmCleanEnergy");
	t->Project("HDCPMT", "PmtCleanEnergy[1]", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2");
	sprintf(str, "MyRandom::GausAdd(PmtCleanEnergy, %6.4f)", kRndm);
	tMc->Project("HMCPMT", "PmtCleanEnergy");
	t->Project("HPH", "PositronEnergy[1]", "N>1 && gtDiff[1]/125<5 && gtDiff[1]/125>1 && Hits[1] < 5");
	t->Project("HXY", "NeutronX[1][1]+2:NeutronX[1][0]+2", "N>1 && gtDiff[1]/125<50 && gtDiff[1]/125>2 && NeutronX[1][1]>=0 && NeutronX[1][0]>=0");
	hcmxy->Fill(18., 2., -1);
	TF1 *fCapt = new TF1("fCapt", CaptureFunction, 0, 20, 10);
	fCapt->SetParNames("Gd:Const", "Gd:#alpha", "Gd:n", "Gd:mean", "Gd:#sigma", "H:Const", "H:#alpha", "H:n", "H:mean", "H:#sigma");
	fCapt->SetLineColor(kRed);
	fCapt->SetLineWidth(2);
	
	TF1 *fMC = new TF1("fMC", MCFunction, 0, 20, 8);
	fMC->SetParNames("H:Const", "H:#alpha", "H:n", "H:mean", "H:#sigma", "P0", "P1", "P2");
	fMC->SetLineColor(kRed);
	fMC->SetLineWidth(2);
	
	TF1 *fH = new TF1("fH", CBfunction, 0, 20, 5);
	fH->SetParNames("Const", "#alpha", "n", "mean", "#sigma");
	fH->SetLineColor(kRed);
	fH->SetLineWidth(2);
	
	TF1 *fGP2 = new TF1("fGP2", "gaus(0) + pol2(3)", 0, 10);
	fGP2->SetParNames("Const", "Mean", "#sigma", "P0", "P1", "P2");
	fGP2->SetLineColor(kRed);
	fGP2->SetLineWidth(2);

	TCanvas *cv = new TCanvas("CV", "Neutron Capture", 1200, 900);
	cv->Divide(2, 2);
	cv->cd(1);
//	fCapt->SetParameters(hcm->GetMaximum(), 1, 10, 6.7, 1, hcm->GetMaximum()/5, 1, 10, 1.8, 0.5);
	fGP2->SetParameters(hcm->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hcm->Fit(fGP2, "", "0", 0.9, 4);
	hcm->DrawCopy();
	TVirtualPad *pd2 = cv->cd(2);
	pd2->SetRightMargin(0.16);
	hcmxy->Draw("COLZ");
	cv->cd(3);
//	fH->SetParameters(hpH->GetMaximum(), 1, 1, 1.7, 1);
	fGP2->SetParameters(hpH->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hpH->Fit(fGP2, "", "", 0.9, 4);
	cv->cd(4);
//	fMC->SetParameters(hMc->GetMaximum()/5, 1, 1, 2.0, 0.5, 0, 0, 0);
	fGP2->SetParameters(hMc->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hMc->Fit(fGP2, "", "", 0.9, 4);
	cv->SaveAs("248Cm.pdf(");
	
	TCanvas *cvA = new TCanvas("CVA", "Neutron Capture", 1200, 900);
	cvA->Divide(2, 2);
	cvA->cd(1);
//	fCapt->SetParameters(hcmSi->GetMaximum(), 1, 5, 6.7, 1, hcmSi->GetMaximum()/5, 1, 5, 1.8, 0.5);
	fGP2->SetParameters(hcmSi->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hcmSi->Fit(fGP2, "", "", 0.9, 4);
	cvA->cd(2);
//	fCapt->SetParameters(hcmPMT->GetMaximum(), 1, 10, 6.7, 1, hcmPMT->GetMaximum()/5, 1, 10, 1.8, 0.5);
	fGP2->SetParameters(hcmPMT->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hcmPMT->Fit(fGP2, "", "", 0.9, 4);
	cvA->cd(3);
//	fMC->SetParameters(hMcSi->GetMaximum()/5, 1, 1, 2.0, 0.5, 0, 0, 0);
	fGP2->SetParameters(hMcSi->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hMcSi->Fit(fGP2, "", "", 0.9, 4);
	cvA->cd(4);
//	fMC->SetParameters(hMcPMT->GetMaximum()/5, 1, 1, 2.0, 0.5, 0, 0, 0);
	fGP2->SetParameters(hMcPMT->GetMaximum()/5, 2, 0.5, 0, 0, 0);
	hMcPMT->Fit(fGP2, "", "", 0.9, 4);
	cvA->SaveAs("248Cm.pdf)");
	
	TCanvas *prl = new TCanvas("PRL", "Neutron capture", 800, 800);
	prl->cd();
	prl->SetLeftMargin(0.20);
	prl->SetRightMargin(0.03);
	prl->SetTopMargin(0.03);
	prl->SetBottomMargin(0.12);
	hcm->GetYaxis()->SetTitleOffset(1.9);
	hcm->SetStats(0);
	hcm->SetFillStyle(kNone);
	hcm->SetLineColor(kBlack);
	hcm->SetMarkerStyle(kNone);
	hcm->SetTitle(";Delayed energy, MeV;Events/100 keV");
	hcm->Draw("e");
	hMc->SetFillStyle(kNone);
	hMc->Scale(hcm->Integral(60, 120) / hMc->Integral(60, 120));
	hMc->Draw("same,hist");
	TLegend *lg = new TLegend(0.3, 0.75, 0.5, 0.85);
	lg->AddEntry(hcm, "Experiment", "LE");
	sprintf(str, "MC+%2.0f%%/#sqrt{E}", 100*kRndm);
	lg->AddEntry(hMc, str, "L");
	lg->Draw();
	prl->Update();
	
	f.Close();
}
