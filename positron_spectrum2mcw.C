#include <TRandom2.h>
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
};

double CalculateScale(TH1 *hA, TH1 *hB, int iMin, int iMax)
{
	int i;
	double SumB2, SumAB;
	
	SumB2 = SumAB = 0;
	for (i=iMin; i<=iMax; i++) {
		SumAB += hA->GetBinContent(i) * hB->GetBinContent(i);
		SumB2 += hB->GetBinContent(i) * hB->GetBinContent(i);
	}
	return SumAB / SumB2;
}

void positron_spectrum2mcw(int RangeMin, int RangeMax, int NormMin, int NormMax, float norm = 1.0, double sigma = 0)
{
	const char fuel[4][6] = {"235U", "238U", "239Pu", "241Pu"};
	const Color_t fColor[4] = {kRed, kBlue, kGreen, kOrange};
//	const float crossection[4] = {6.39, 8.90, 4.18, 5.76};	// from Sinev
	const float crossection[4] = {6.69, 10.10, 4.36, 6.04};	// from H-M as quoted by DB 1707.07728
//	const char mcfile[4][128] = {
//		"/space/danss_root3/withdead-uncorr/mc_positron_235U_simple_newScale.root",
//		"/space/danss_root3/withdead-uncorr/mc_positron_238U_simple_newScale.root",
//		"/space/danss_root3/withdead-uncorr/mc_positron_239Pu_simple_newScale.root",
//		"/space/danss_root3/withdead-uncorr/mc_positron_241Pu_simple_newScale.root"
//	};
	const char mcfile[4][128] = {
		"mc_ibdNewGd_235U_transcode_pair.root",
		"mc_ibdNewGd_238U_transcode_pair.root",
		"mc_ibdNewGd_239Pu_transcode_pair.root",
		"mc_ibdNewGd_241Pu_transcode_pair.root"
	};
//	From Sinev
//	const double fuelmix[3][4]  = {{0.69, 0.07, 0.21, 0.03}, {0.58, 0.07, 0.30, 0.05}, {0.47, 0.07, 0.39, 0.07}};
//	From Khvatov
//	const double fuelmix[3][4]  = {{0.78, 0.07, 0.12, 0.03}, {0.69, 0.07, 0.19, 0.05}, {0.59, 0.07, 0.26, 0.07}};
//		Replace end with 2 month before the end and begin with 5th campaign begin - our calculations
//	const double fuelmix[3][4]  = {{0.795, 0.07, 0.111, 0.024}, {0.69, 0.07, 0.19, 0.05}, {0.613, 0.079, 0.244, 0.065}};
//		KNPP calculations       begin 5                      Sinev (old) middle        End 4 - 2 month
//	const double fuelmix[3][4]  = {{0.661, 0.067, 0.249, 0.023}, {0.58, 0.07, 0.30, 0.05}, {0.474, 0.074, 0.371, 0.077}};
//	KNPP fuel contribution calculations Begin 4			End 4				Begin 5
	const double fuelmix0[3][4]  = {{0.637, 0.068, 0.266, 0.028}, {0.447, 0.075, 0.389, 0.085}, {0.661, 0.067, 0.249, 0.023}};
//		Calculated   begin 5 + 2.5 months           middle (not used)    End 4 - 1.5 month Coeffs
	double fuelmix[3][4];
	const char cmppart[3][20] = {"Begin", "Middle", "End"};
	TString expname("danss_report_v4_mar18-calc");
	TFile *fMc[4];
	TTree *tMc[4];
	TH1D *hMc[4];
	TH1D *hMcMixt[3];
	TH1D *hMcRatio;
	TFile *fExp;
	TH1D *hExpw;
	TH1D *hExp;
	TH1D *hExpRatio;
	TH1D *hRatio;
	TH1D *hDiff;
	TPaveStats *pv;
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cXYZ = cX && cY && cZ && cGamma;
	int i, j;
	char strs[128];
	char strl[1024];
	double s, term, A, B;

//		Calculated   begin 5 + 2.5 months     middle (not used)    End 4 - 1.5 month
	for (i=0; i<4; i++) {
//		begin 5 + 2.5 months => <begin 5> * (1- 2.5*30/471) + <end 4> * (2.5*30/471)
		term = 2.5 * 30 / 471;
		fuelmix[0][i] = fuelmix0[2][i] * (1 - term) + fuelmix0[1][i] * term;
//		middle (not used)
		fuelmix[1][i] = (fuelmix0[0][i] + fuelmix0[1][i]) / 2.0;
//		End 4 - 1.5 month
		term = 1.5 * 30 / 471;
		fuelmix[2][i] = fuelmix0[1][i] * (1 - term) + fuelmix0[0][i] * term;
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);

	fExp = new TFile(expname + ".root");
	if (!fExp->IsOpen()) return;
	hExpw = (TH1D*) fExp->Get("hSumAfter");
	if (!hExpw) return;
	hExp = new TH1D("hExp", "Experimental positron spectrum;MeV;Events/(day*0.25 MeV)", 
		RangeMax - RangeMin + 1, hExpw->GetBinLowEdge(RangeMin), hExpw->GetBinLowEdge(RangeMax) + hExpw->GetBinWidth(RangeMax));
	for (i=0; i<RangeMax - RangeMin + 1; i++) {
		hExp->SetBinContent(i+1, hExpw->GetBinContent(i + RangeMin));
		hExp->SetBinError(i+1, hExpw->GetBinError(i + RangeMin));
	}
	hExp->SetLineColor(kBlack);
	hExp->SetMarkerColor(kBlack);
	hExp->SetMarkerStyle(20);
	hExp->SetLineWidth(4);

	hExpRatio = (TH1D*) fExp->Get("hNRatioAfter2Before");
	if (!hExpRatio) return;
	hExpRatio->Scale(1.004);	// Dead time & Power correction
	hExpRatio->SetLineWidth(3);
	hExpRatio->SetLineColor(kGreen);
	hExpRatio->GetXaxis()->SetRange(RangeMin, RangeMax);
	hExpRatio->SetMinimum(1.0);
	hExpRatio->SetMaximum(1.2);
	hExpRatio->GetYaxis()->SetLabelSize(0.05);

	gROOT->cd();
	for (i=0; i<4; i++) {
		sprintf(strs, "h%s", fuel[i]);
		sprintf(strl, "Positron spectrum of %s;MeV", fuel[i]);
		hMc[i] = new TH1D(strs, strl, RangeMax - RangeMin + 1, hExpw->GetBinLowEdge(RangeMin), hExpw->GetBinLowEdge(RangeMax) + hExpw->GetBinWidth(RangeMax));
		hMc[i]->SetLineColor(fColor[i]);
	}
	for (i=0; i<4; i++) {
		fMc[i] = new TFile(mcfile[i]);
		if (!fMc[i]->IsOpen()) break;
		tMc[i] = (TTree *) fMc[i]->Get("DanssPair");
	}
	if (i != 4) {
		printf("Mc file %s error\n", mcfile[i]);
		return;
	}
	gROOT->cd();
	sprintf(strs, "%6.3f*MyRandom::GausAdd(PositronEnergy, %8.5f)", norm, sigma);
	for (i=0; i<4; i++) tMc[i]->Project(hMc[i]->GetName(), strs, cXYZ);
	for (i=0; i<4; i++) hMc[i]->Sumw2();
	for (i=0; i<4; i++) hMc[i]->Scale(crossection[i]);
	for (j=0; j<3; j++) {
		sprintf(strs, "h%sMixt", cmppart[j]);
		hMcMixt[j] = (TH1D*) hMc[0]->Clone(strs);
		hMcMixt[j]->Reset();
		for (i=0; i<4; i++) hMcMixt[j]->Add(hMc[i], fuelmix[j][i]);
		hMcMixt[j]->SetLineColor(fColor[j]);
		hMcMixt[j]->SetLineWidth(2);
		hMcMixt[j]->SetMarkerSize(3);
		hMcMixt[j]->SetMarkerStyle(kFullStar);
		hMcMixt[j]->SetMarkerColor(fColor[j]);
		hMcMixt[j]->GetYaxis()->SetLabelSize(0.05);
		sprintf(strs, "MC spectrum %s;Positron energy, MeV;Events/(day*0.25 MeV)", cmppart[j]);
		hMcMixt[j]->SetTitle(strs);
	}
	hMcRatio = (TH1D*) hMc[0]->Clone("hMcBegin2End");
	hMcRatio->Reset();
	hMcRatio->Divide(hMcMixt[0], hMcMixt[2]);
	for (i=0; i<hMcRatio->GetNbinsX(); i++) {
		s = 0;
		A = hMcMixt[0]->GetBinContent(i + 1);
		B = hMcMixt[2]->GetBinContent(i + 1);
		for (j=0; j<4; j++) {
			term = fuelmix[0][j] * B - fuelmix[2][j] * A;
			term *= hMc[j]->GetBinError(i+1);
			s += term*term;
		}
		s = sqrt(s) / (B * B);
		hMcRatio->SetBinError(i+1, s);
	}
	hMcRatio->SetLineColor(kBlue);
	hMcRatio->SetLineWidth(2);
	hMcRatio->SetMarkerStyle(kFullStar);
	hMcRatio->SetMarkerColor(kBlue);
	hMcRatio->SetMarkerSize(3);

	hMcRatio->GetYaxis()->SetLabelSize(0.05);
	hMcRatio->SetTitle("Begin to end spectrum ratio;Positron energy, MeV;Events/(day*0.25 MeV)");
	
	TFile fSave("mc_fuel.root", "RECREATE");
	fSave.cd();
	for (j=0; j<4; j++) hMc[j]->Write();
	for (j=0; j<3; j++) hMcMixt[j]->Write();
	hMcRatio->Write();
	fSave.Close();
	
	for (j=0; j<3; j++) hMcMixt[j]->Scale(CalculateScale(hExp, hMcMixt[j], NormMin, NormMax));
	hRatio = (TH1D *) hMcMixt[0]->Clone("hRatioExpMc");
	hRatio->Divide(hExp, hMcMixt[0]);
	hRatio->SetTitle(";Positron energy, MeV;#frac{N_{EXP}}{N_{MC}}");

	hDiff = (TH1D *) hMcMixt[0]->Clone("hRatioExpMc");
	hDiff->Add(hExp, hMcMixt[0], 1, -1);
	hDiff->SetTitle("Experiment - Monte Carlo;Positron energy, MeV;#Delta Events/(day*0.25 MeV)");
	
	TCanvas *cm = new TCanvas("CM", "Fuel", 1200, 900);
	hMc[1]->SetTitle("DANSS simulated spectrum per isotope;MeV;");
	hMc[1]->Draw("hist");
	hMc[0]->Draw("hist,same");
	hMc[2]->Draw("hist,same");
	hMc[3]->Draw("hist,same");
	TLegend *lm = new TLegend(0.5, 0.7, 0.9, 0.9);
	for (i=0; i<4; i++) {
		sprintf(strs, "%s - %4.1f%% (before) / %4.1f%% (after)", fuel[i], 100*fuelmix[2][i], 100*fuelmix[0][i]);
		lm->AddEntry(hMc[i], strs, "L");
	}
	lm->Draw();
	TLatex *txt = new TLatex();
	txt->SetTextSize(0.035);
	sprintf(strs, "MC Scale=%5.3f Random+=%2.0f%%", norm, sigma*100);
	txt->DrawLatex(4.3, 95000, strs);
	sprintf(strs, "_S%5.3f_R%4.2f", norm, sigma);
	expname += strs;
	cm->SaveAs(expname + "-pos2mc.pdf[");
	cm->SaveAs(expname + "-pos2mc.pdf");
	
	TCanvas *cv = new TCanvas("CV", "Exp & MC", 1200, 900);
//	for (j=0; j<3; j++) hMcMixt[j]->Draw((j) ? "same,hist" : "hist");
	hMcMixt[0]->SetTitle("3 months period one month after campaign begin");
	hMcMixt[0]->Draw("E");
	hExp->Draw("same");
	TLegend *lg = new TLegend(0.6, 0.75, 0.9, 0.9);
	lg->AddEntry(hExp, "DANSS data", "LP");
	lg->AddEntry(hMcMixt[0], "Monte Carlo", "LP");
//	for (j=0; j<3; j++) {
//		sprintf(strs, "MC - %s", cmppart[j]);
//		lg->AddEntry(hMcMixt[j], strs, "L");
//	}
	lg->Draw();
	cv->SaveAs(expname + "-pos2mc.pdf");
	
	TCanvas *cr = new TCanvas("CR", "Ratio", 1200, 500);
	hRatio->GetYaxis()->SetLabelSize(0.08);
	hRatio->GetXaxis()->SetLabelSize(0.08);
	hRatio->GetYaxis()->SetTitleSize(0.08);
	hRatio->GetXaxis()->SetTitleSize(0.08);
	hRatio->SetLineWidth(4);
	hRatio->SetMinimum(0.92);
	hRatio->SetMaximum(1.15);
	hRatio->Draw();
	cr->SaveAs(expname + "-pos2mc.pdf");
	
	TCanvas *cd = new TCanvas("CD", "Difference", 1200, 500);
	hDiff->GetYaxis()->SetLabelSize(0.07);
	hDiff->GetXaxis()->SetLabelSize(0.07);
	hDiff->GetYaxis()->SetTitleSize(0.07);
	hDiff->GetXaxis()->SetTitleSize(0.07);
	hDiff->SetLineWidth(4);
	hDiff->SetMinimum(-10);
	hDiff->SetMaximum(10);
	hDiff->Draw();
	cd->SaveAs(expname + "-pos2mc.pdf");

	TCanvas *cnr = new TCanvas("CNR", "Spectrum ratio: Exp vs MC", 1200, 900);
	TH1D *hExpRatioR = (TH1D *)hExpRatio->Clone("hExpRatioR");
	TH1D *hOneE = (TH1D *)hExpRatio->Clone("hOne");
	TH1D *hMcRatioR = (TH1D *)hMcRatio->Clone("hMcRatioR");
	TH1D *hOneM = (TH1D *)hMcRatio->Clone("hOne");
	hExpRatioR->SetTitle("Normalized ratio before shutdown / after shutdown");
	for (i=hOneE->GetXaxis()->GetFirst(); i<=hOneE->GetXaxis()->GetLast(); i++) {
		hOneE->SetBinContent(i, 1);
		hOneE->SetBinError(i, 1.0E-10);
	}
	for (i=hOneM->GetXaxis()->GetFirst(); i<=hOneM->GetXaxis()->GetLast(); i++) {
		hOneM->SetBinContent(i, 1);
		hOneM->SetBinError(i, 1.0E-10);
	}
	hExpRatioR->Divide(hOneE, hExpRatio);
	hMcRatioR->Divide(hOneM, hMcRatio);
	hExpRatioR->SetMinimum(0.75);
	hExpRatioR->SetMaximum(1);
	hExpRatioR->Draw();
	hMcRatioR->Draw("same");
	lg = new TLegend(0.2, 0.25, 0.5, 0.4);
	lg->AddEntry(hExpRatioR, "DANSS data", "LE");
	lg->AddEntry(hMcRatioR, "Monte Carlo", "LPE");
	lg->Draw();
	cnr->SaveAs(expname + "-pos2mc.pdf");

	cnr->Clear();
	hExpRatio->Draw();
	hMcRatio->Draw("same");
	lg = new TLegend(0.2, 0.7, 0.5, 0.85);
	lg->AddEntry(hExpRatio, "DANSS data", "LE");
	lg->AddEntry(hMcRatio, "Monte Carlo", "LPE");
	lg->Draw();
	cnr->SaveAs(expname + "-pos2mc.pdf");

	cnr->Clear();
	gStyle->SetOptFit();
	hExpRatio->SetStats();
	hMcRatio->SetStats();
	hExpRatio->Fit("pol1", "", "", 1, 6);
	hExpRatio->GetFunction("pol1")->SetLineColor(kGreen);
	hMcRatio->Fit("pol1", "", "sames", 1, 6);
	hMcRatio->GetFunction("pol1")->SetLineColor(kBlue);
	lg->Draw();
	cnr->Update();
	pv = (TPaveStats *) hExpRatio->FindObject("stats");
	pv->SetLineColor(kGreen);
	pv->SetX1NDC(0.46);
	pv->SetX2NDC(0.66);
	pv->SetY1NDC(0.18);
	pv->SetY2NDC(0.32);
	
	pv = (TPaveStats *) hMcRatio->FindObject("stats");
	pv->SetLineColor(kBlue);
	pv->SetX1NDC(0.68);
	pv->SetX2NDC(0.88);
	pv->SetY1NDC(0.18);
	pv->SetY2NDC(0.32);

//		Calculate various chi2
	TF1 *fPol1 = new TF1("fPol1", "pol1", 1, 6);
	txt->SetTextSize(0.02);
	fPol1->FixParameter(0, 1);
	fPol1->FixParameter(1, 0);
	sprintf(strl, "Const=1.000  #chi^{2}/ndf = %7.2f/20", hExpRatio->Chisquare(fPol1, "r"));
	txt->DrawLatex(1.5, 1.01, strl);
	fPol1->FixParameter(0, 1.015);
	sprintf(strl, "Const=1.015  #chi^{2}/ndf = %7.2f/20", hExpRatio->Chisquare(fPol1, "r"));
	txt->DrawLatex(1.5, 1.02, strl);
	fPol1->FixParameter(0, 1.02);
	sprintf(strl, "Const=1.02  #chi^{2}/ndf = %7.2f/20", hExpRatio->Chisquare(fPol1, "r"));
	txt->DrawLatex(1.5, 1.03, strl);

	hExpRatio->Fit("pol0", "Q+", "", 1, 6);
	sprintf(strl, "Const=%5.3f  #chi^{2}/ndf = %7.2f/19", 
		((TF1*)hExpRatio->FindObject("pol0"))->GetParameter(0),
		((TF1*)hExpRatio->FindObject("pol0"))->GetChisquare());
	txt->DrawLatex(1.5, 1.12, strl);

	fPol1->FixParameter(0, ((TF1*)hMcRatio->FindObject("pol1"))->GetParameter(0));
	fPol1->FixParameter(1, ((TF1*)hMcRatio->FindObject("pol1"))->GetParameter(1));
	sprintf(strl, "MC   #chi^{2}/ndf = %7.2f/20", hExpRatio->Chisquare(fPol1, "r"));
	txt->DrawLatex(1.5, 1.13, strl);

	cnr->Update();
	cnr->SaveAs(expname + "-pos2mc.pdf");
	cnr->SaveAs(expname + "-pos2mc.pdf]");
}

void positron_spectrum2mcw3(int RangeMin, int RangeMax, int NormMin, int NormMax)
{
	const char fuel[4][6] = {"235U", "238U", "239Pu", "241Pu"};
	const char pos[3][6] = {"UP", "MID", "DOWN"};
	const Color_t fColor[4] = {kRed, kBlue, kGreen, kOrange};
	const char mcfile[4][128] = {
		"/space/danss_root3/withdead-uncorr/mc_positron_235U_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_238U_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_239Pu_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_241Pu_simple_newScale.root"
	};
	const double fuelmix[3][4]  = {{0.69, 0.07, 0.21, 0.03}, {0.58, 0.07, 0.30, 0.05}, {0.47, 0.07, 0.39, 0.07}};
	const char cmppart[3][20] = {"Begin", "Middle", "End"};
	TString expname("danss_report_v4n-calc");
	TFile *fMc[4];
	TTree *tMc[4];
	TH1D *hMc[4];
	TH1D *hMcMixt[3];
	TFile *fExp;
	TH1D *hExpw[3];
	TH1D *hExp[3];
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cXYZ = cX && cY && cZ && cGamma;
	int i, j;
	char strs[128];
	char strl[1024];
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);

	fExp = new TFile(expname + ".root");
	if (!fExp->IsOpen()) return;
	hExpw[0] = (TH1D*) fExp->Get("hUp_28");
	hExpw[1] = (TH1D*) fExp->Get("hMid_28");
	hExpw[2] = (TH1D*) fExp->Get("hDown_28");
	if (!hExpw[0] || !hExpw[1] || !hExpw[2]) return;
	for (i=0; i<3;i++) {
		sprintf(strs, "%s-r", hExpw[i]->GetName());
		hExp[i] = new TH1D(strs, "Experimental positron spectrum;MeV;Events/(day*0.25 MeV)", 
			RangeMax - RangeMin + 1, hExpw[i]->GetBinLowEdge(RangeMin), hExpw[i]->GetBinLowEdge(RangeMax) + hExpw[i]->GetBinWidth(RangeMax));
	}
	for (j=0; j<3; j++) for (i=0; i<RangeMax - RangeMin + 1; i++) {
		hExp[j]->SetBinContent(i+1, hExpw[j]->GetBinContent(i + RangeMin));
		hExp[j]->SetBinError(i+1, hExpw[j]->GetBinError(i + RangeMin));
	}
	hExp[0]->SetLineColor(kRed);
	hExp[0]->SetFillColor(kRed-10);
	hExp[1]->SetLineColor(kGreen);
	hExp[1]->SetFillColor(kGreen-10);
	hExp[2]->SetLineColor(kBlue);
	hExp[2]->SetFillColor(kBlue-10);
//	hExp->SetMarkerColor(kBlack);
//	hExp->SetMarkerStyle(20);
//	hExp->SetLineWidth(4);
	gROOT->cd();
	for (i=0; i<4; i++) {
		sprintf(strs, "h%s", fuel[i]);
		sprintf(strl, "Positron spectrum of %s;MeV", fuel[i]);
		hMc[i] = new TH1D(strs, strl, RangeMax - RangeMin + 1, hExpw[0]->GetBinLowEdge(RangeMin), hExpw[0]->GetBinLowEdge(RangeMax) + hExpw[0]->GetBinWidth(RangeMax));
		hMc[i]->SetLineColor(fColor[i]);
	}
	for (i=0; i<4; i++) {
		fMc[i] = new TFile(mcfile[i]);
		if (!fMc[i]->IsOpen()) break;
		tMc[i] = (TTree *) fMc[i]->Get("DanssEvent");
	}
	if (i != 4) {
		printf("Mc file %s error\n", mcfile[i]);
		return;
	}
	gROOT->cd();
	for (i=0; i<4; i++) tMc[i]->Project(hMc[i]->GetName(), "1.07*(PositronEnergy-0.179)/0.929", cXYZ);
	for (i=0; i<4; i++) hMc[i]->Sumw2();
	for (j=0; j<3; j++) {
		sprintf(strs, "h%sMixt", cmppart[j]);
		hMcMixt[j] = (TH1D*) hMc[0]->Clone(strs);
		hMcMixt[j]->Reset();
		for (i=0; i<4; i++) hMcMixt[j]->Add(hMc[i], fuelmix[j][i]);
		hMcMixt[j]->Scale(hExp[2]->Integral(NormMin, NormMax) / hMcMixt[j]->Integral(NormMin, NormMax));
		hMcMixt[j]->SetLineColor(fColor[j]);
		hMcMixt[j]->SetLineWidth(4);
		hMcMixt[j]->GetYaxis()->SetLabelSize(0.1);
		hMcMixt[j]->SetTitle(";Positron energy, MeV;Events/(day*0.25 MeV)");
	}

	TCanvas *cv = new TCanvas("CV", "Exp & MC", 1200, 900);
//	for (j=0; j<3; j++) hMcMixt[j]->Draw((j) ? "same,hist" : "hist");
	hExp[0]->Draw("hist,e");
	hExp[1]->Draw("hist,e,same");
	hExp[2]->Draw("hist,e,same");
	hMcMixt[1]->SetLineColor(kBlack);
	hMcMixt[1]->SetMarkerStyle(kStar);
	hMcMixt[1]->SetMarkerSize(0.02);
	hMcMixt[1]->Draw("same");
	TLegend *lg = new TLegend(0.6, 0.75, 0.9, 0.9);
	for (i=0; i<3; i++) {
		sprintf(strs, "DANSS %s", pos[i]);
		lg->AddEntry(hExp[i], strs, "LF");
	}
	lg->AddEntry(hMcMixt[1], "Monte Carlo", "LP");
//	for (j=0; j<3; j++) {
//		sprintf(strs, "MC - %s", cmppart[j]);
//		lg->AddEntry(hMcMixt[j], strs, "L");
//	}
	lg->Draw();
}
