void positron_spectrum2mcw(int RangeMin, int RangeMax, int NormMin, int NormMax)
{
	const char fuel[4][6] = {"235U", "238U", "239Pu", "241Pu"};
	const Color_t fColor[4] = {kRed, kBlue, kGreen, kOrange};
	const char mcfile[4][128] = {
		"/space/danss_root3/withdead-uncorr/mc_positron_235U_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_238U_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_239Pu_simple_newScale.root",
		"/space/danss_root3/withdead-uncorr/mc_positron_241Pu_simple_newScale.root"
	};
	const double fuelmix[3][4]  = {{0.69, 0.07, 0.21, 0.03}, {0.58, 0.07, 0.30, 0.05}, {0.47, 0.07, 0.39, 0.07}};
	const char cmppart[3][20] = {"Begin", "Middle", "End"};
	const char expname[] = "danss_report_v4n-calc.root";
	TFile *fMc[4];
	TTree *tMc[4];
	TH1D *hMc[4];
	TH1D *hMcMixt[3];
	TFile *fExp;
	TH1D *hExpw;
	TH1D *hExp;
	TH1D *hRatio;
	TH1D *hDiff;
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

	fExp = new TFile(expname);
	if (!fExp->IsOpen()) return;
	hExpw = (TH1D*) fExp->Get("hSum4");
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
		tMc[i] = (TTree *) fMc[i]->Get("DanssEvent");
	}
	if (i != 4) {
		printf("Mc file %s error\n", mcfile[i]);
		return;
	}
	gROOT->cd();
	for (i=0; i<4; i++) tMc[i]->Project(hMc[i]->GetName(), "1.04*(PositronEnergy-0.179)/0.929", cXYZ);
	for (i=0; i<4; i++) hMc[i]->Sumw2();
	for (j=0; j<3; j++) {
		sprintf(strs, "h%sMixt", cmppart[j]);
		hMcMixt[j] = (TH1D*) hMc[0]->Clone(strs);
		hMcMixt[j]->Reset();
		for (i=0; i<4; i++) hMcMixt[j]->Add(hMc[i], fuelmix[j][i]);
		hMcMixt[j]->Scale(hExp->Integral(NormMin, NormMax) / hMcMixt[j]->Integral(NormMin, NormMax));
		hMcMixt[j]->SetLineColor(fColor[j]);
		hMcMixt[j]->SetLineWidth(2);
		hMcMixt[j]->GetYaxis()->SetLabelSize(0.05);
		hMcMixt[j]->SetTitle(";Positron energy, MeV;Events/(day*0.25 MeV)");
	}
	TFile fSave("mc_fuel.root", "RECREATE");
	fSave.cd();
	for (j=0; j<3; j++) hMcMixt[j]->Write();
	fSave.Close();
	
	hRatio = (TH1D *) hMcMixt[1]->Clone("hRatioExpMc");
	hRatio->Divide(hExp, hMcMixt[2]);
	hRatio->SetTitle(";Positron energy, MeV;#frac{N_{EXP}}{N_{MC}}");

	hDiff = (TH1D *) hMcMixt[2]->Clone("hRatioExpMc");
	hDiff->Add(hExp, hMcMixt[2], 1, -1);
	
	TCanvas *cm = new TCanvas("CM", "Fuel", 1200, 900);
	hMc[2]->Draw("hist");
	hMc[0]->Draw("hist,same");
	hMc[1]->Draw("hist,same");
	hMc[3]->Draw("hist,same");
	TLegend *lm = new TLegend(0.7, 0.7, 0.9, 0.9);
	for (i=0; i<4; i++) lm->AddEntry(hMc[i], fuel[i], "L");
	lm->Draw();
	
	TCanvas *cv = new TCanvas("CV", "Exp & MC", 1200, 900);
//	for (j=0; j<3; j++) hMcMixt[j]->Draw((j) ? "same,hist" : "hist");
	hMcMixt[1]->Draw("hist");
	hExp->Draw("same");
	TLegend *lg = new TLegend(0.6, 0.75, 0.9, 0.9);
	lg->AddEntry(hExp, "DANSS data", "LP");
	lg->AddEntry(hMcMixt[2], "Monte Carlo", "L");
//	for (j=0; j<3; j++) {
//		sprintf(strs, "MC - %s", cmppart[j]);
//		lg->AddEntry(hMcMixt[j], strs, "L");
//	}
	lg->Draw();
	
	TCanvas *cr = new TCanvas("CR", "Ratio", 1200, 500);
	hRatio->GetYaxis()->SetLabelSize(0.08);
	hRatio->GetXaxis()->SetLabelSize(0.08);
	hRatio->GetYaxis()->SetTitleSize(0.08);
	hRatio->GetXaxis()->SetTitleSize(0.08);
	hRatio->SetLineWidth(4);
	hRatio->SetMinimum(0.92);
	hRatio->SetMaximum(1.15);
	hRatio->Draw();
	
	TCanvas *cd = new TCanvas("CD", "Difference", 1200, 500);
	hDiff->GetYaxis()->SetLabelSize(0.07);
	hDiff->GetXaxis()->SetLabelSize(0.07);
	hDiff->GetYaxis()->SetTitleSize(0.07);
	hDiff->GetXaxis()->SetTitleSize(0.07);
	hDiff->SetLineWidth(4);
	hDiff->SetMinimum(-10);
	hDiff->SetMaximum(10);
	hDiff->Draw();
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
	const char expname[] = "danss_report_v4n-calc.root";
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

	fExp = new TFile(expname);
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
	for (i=0; i<4; i++) tMc[i]->Project(hMc[i]->GetName(), "1.04*(PositronEnergy-0.179)/0.929", cXYZ);
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
