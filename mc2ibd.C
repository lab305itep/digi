TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
TCut cR1("Distance < 45");
TCut cR2("Distance < 55");
TCut cR = cR2 && (cRXY || cR1);
TCut c20("gtDiff > 2");
TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
TCut cGammaMax("AnnihilationMax < 0.8");
TCut cPe("PositronEnergy > 1");
TCut cN("NeutronEnergy > 3.5");

void mc2ibd_gtDiff(TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	TLatex *txt = new TLatex();

	TH1D *hExp = (TH1D *) fBgnd->Get("hgtDiffA-diff");
	if (!hExp) {
		printf("Histogram hgtDiffA-diff not found in %s\n", fBgnd->GetName());
		return;
	}
	hExp->SetTitle("Time between prompt and delayed events (IBD);us;Events/us");
	hExp->SetLineColor(kBlack);
	hExp->SetMinimum(0);

	gROOT->cd();
	TH1D *hMC = new TH1D("hgtDiffMC", "Time between prompt and delayed events (MC);us;Events/us", 50, 0, 50);
	hMC->SetLineColor(kBlack);
	tMC->Project(hMC->GetName(), "gtDiff", cX && cY && cZ && cR && cGamma && cGammaMax && cPe && cN);
	hMC->Sumw2();
	
	TF1 *f2Exp = new TF1("f2Exp", "[0]*(exp(-x/[1]) - exp(-x/[2]))", 0, 100);
	f2Exp->SetParNames("Const.", "#tau_{capt}", "#tau_{th}");
	f2Exp->SetLineColor(kBlue);
	
	TF1 *fExp = new TF1("fExp", "expo", 0, 100);
	fExp->SetLineColor(kRed);
	
	cv->Clear();
	cv->Divide(1, 2);
	cv->cd(1);
	f2Exp->SetParameters(hExp->GetMaximum(), 13, 5);
	hExp->Fit(f2Exp, "", "", 1, 50);
	hExp->Fit(fExp, "+", "", 15, 50);
	TLegend *lg = new TLegend(0.15, 0.2, 0.4, 0.35);
	lg->AddEntry(hExp, "Data", "LE");
	lg->AddEntry(f2Exp, "Fit with two exponents", "L");
	lg->AddEntry(fExp, "Capture exponent fit", "L");
	lg->Draw();
	sprintf(str, "#tau_{capt}=%5.2f us", -1/fExp->GetParameter(1));
	txt->DrawLatex(20, hExp->GetMaximum()/5, str);
	sprintf(str, "#chi^{2}/NDF=%6.1f/34", fExp->GetChisquare());
	txt->DrawLatex(20, hExp->GetMaximum()/10, str);
	
	cv->cd(2);
	f2Exp->SetParameters(hMC->GetMaximum(), 13, 5);
	hMC->Fit(f2Exp, "", "", 1, 50);
	hMC->Fit(fExp, "+", "", 15, 50);
	lg->Draw();
	sprintf(str, "#tau_{capt}=%5.2f us", -1/fExp->GetParameter(1));
	txt->DrawLatex(20, hMC->GetMaximum()/5, str);
	sprintf(str, "#chi^{2}/NDF=%6.1f/34", fExp->GetChisquare());
	txt->DrawLatex(20, hMC->GetMaximum()/10, str);
}

void mc2ibd_R1(TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	TLatex *txt = new TLatex();

	TH1D *hExp = (TH1D *) fBgnd->Get("hR1A-diff");
	if (!hExp) {
		printf("Histogram hR1A-diff not found in %s\n", fBgnd->GetName());
		return;
	}
	hExp->SetTitle("Distance between positron and neutron, 2D case;cm;Events/4cm");
	hExp->SetLineColor(kBlack);
	hExp->SetLineWidth(3);

	gROOT->cd();
	TH1D *hMC = new TH1D("hR1MC", "Distance between positron and neutron, 2D case (MC);cm;Events/4cm", 40, 0, 160);
	hMC->SetLineColor(kBlue);
	tMC->Project(hMC->GetName(), "Distance", cX && cY && cZ && !cRXY && c20 && cGamma && cGammaMax && cPe && cN);
	hMC->Sumw2();
	hMC->Scale(hExp->Integral() / hMC->Integral());
	
	cv->Clear();
	hExp->DrawCopy();
	hMC->Draw("hist,same");
	TLegend *lg = new TLegend(0.65, 0.8, 0.89, 0.89);
	lg->AddEntry(hExp, "IBD", "LE");
	lg->AddEntry(hMC, "MC", "L");
	lg->Draw();
}

void mc2ibd_R2(TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	TLatex *txt = new TLatex();

	TH1D *hExp = (TH1D *) fBgnd->Get("hR2A-diff");
	if (!hExp) {
		printf("Histogram hR2A-diff not found in %s\n", fBgnd->GetName());
		return;
	}
	hExp->SetTitle("Distance between positron and neutron, 3D case;cm;Events/4cm");
	hExp->SetLineColor(kBlack);
	hExp->SetLineWidth(3);

	gROOT->cd();
	TH1D *hMC = new TH1D("hR2MC", "Distance between positron and neutron, 3D case (MC);cm;Events/4cm", 40, 0, 160);
	hMC->SetLineColor(kBlue);
	tMC->Project(hMC->GetName(), "Distance", cX && cY && cZ && cRXY && c20 && cGamma && cGammaMax && cPe && cN);
	hMC->Sumw2();
	hMC->Scale(hExp->Integral() / hMC->Integral());
	
	cv->Clear();
	hExp->DrawCopy();
	hMC->Draw("hist,same");
	TLegend *lg = new TLegend(0.65, 0.8, 0.89, 0.89);
	lg->AddEntry(hExp, "IBD", "LE");
	lg->AddEntry(hMC, "MC", "L");
	lg->Draw();
}

void mc2ibd_NE(TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	TLatex *txt = new TLatex();

	TH1D *hExp = (TH1D *) fBgnd->Get("hNEA-diff");
	if (!hExp) {
		printf("Histogram hNEA-diff not found in %s\n", fBgnd->GetName());
		return;
	}
	hExp->SetTitle("Delayed event energy;MeV;Events/200 keV");
	hExp->SetLineColor(kBlack);
	hExp->SetLineWidth(3);

	gROOT->cd();
	TH1D *hMC = new TH1D("hNEMC", "Delayed event energy (MC);MeV;Events/200 keV", 45, 3, 12);
	hMC->SetLineColor(kBlue);
	tMC->Project(hMC->GetName(), "NeutronEnergy", cX && cY && cZ && cR && c20 && cGamma && cGammaMax && cPe);
	hMC->Sumw2();
	hMC->Scale(hExp->Integral(15, 45) / hMC->Integral(15,45));
	
	cv->Clear();
	hMC->Draw("hist");
	hExp->DrawCopy("same");
	TLegend *lg = new TLegend(0.65, 0.8, 0.89, 0.89);
	lg->AddEntry(hExp, "IBD", "LE");
	lg->AddEntry(hMC, "MC", "L");
	lg->Draw();
}

void mc2ibd_PX(TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	int i;
	TH1D *hExp[3];
	TH1D *hMC[3];
	TH1D *hR[3];
	TCut cut;
	const int Color[] = {kRed, kBlue, kGreen};

	for (i=0; i<3; i++) {
		sprintf(str, "hP%cA-diff", 'X'+i);
		hExp[i] = (TH1D *) fBgnd->Get(str);
		if (!hExp[i]) {
			printf("Histogram hP%cA-diff not found in %s\n", 'X'+i, fBgnd->GetName());
			return;
		}
		gROOT->cd();
		sprintf(str, "hMCP%c", 'X'+i);
		hMC[i] = (TH1D*) hExp[i]->Clone(str);
		hMC[i]->Reset();
		sprintf(str, "PositronX[%d]+%4.1f", i, (i==2) ? 0.5 : 2.0);
		switch (i) {
		case 0:
			cut = cY && "PositronX[0]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
			break;
		case 1:
			cut = cX && "PositronX[1]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
			break;
		default:
			cut = cX && cY && cR && c20 && cGamma && cGammaMax && cPe && cN;
		}
		gROOT->cd();
		tMC->Project(hMC[i]->GetName(), str, cut);
		hMC[i]->Sumw2();
		if (i == 2) {
			hMC[i]->Rebin(4);
			hExp[i]->Rebin(4);
		}
		sprintf(str, "hRatioP%c", 'X'+i);
		hR[i] = (TH1D*) hExp[i]->Clone(str);
		hR[i]->Divide(hExp[i], hMC[i]);
		hR[i]->SetTitle("IBD to MC Ratio;cm;Ratio");
		hR[i]->SetLineColor(Color[i]);
	}
	cv->Clear();
	hR[0]->SetMinimum(0.9*hR[2]->GetMinimum());
	for (i=0; i<3; i++) hR[i]->DrawCopy((i) ? "same" : "");
	TLegend *lg = new TLegend(0.2, 0.7, 0.4, 0.85);
	for (i=0; i<3; i++) {
		sprintf(str, "%c", 'X'+i);
		lg->AddEntry(hR[i], str, "LE");
	}
	lg->Draw();
}

void mc2ibd_PPX(char X, TChain *tMC, TFile *fBgnd, TCanvas *cv)
{
	char str[1024];
	TH1D *hExp;
	TH1D *hMC;
	TCut cut;

	sprintf(str, "hP%cA-diff", X);
	hExp = (TH1D *) fBgnd->Get(str);
	if (!hExp) {
		printf("Histogram hP%cA-diff not found in %s\n", X, fBgnd->GetName());
		return;
	}
	gROOT->cd();
	sprintf(str, "hMCPP%c", X);
	hMC = new TH1D(str, "MC XYZ", hExp->GetNbinsX(), 0, 100);
	sprintf(str, "PositronX[%d]+%4.1f", X - 'X', (X=='Z') ? 0.5 : 2.0);
	switch (X) {
	case 'X':
		cut = cY && "PositronX[0]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
		break;
	case 'Y':
		cut = cX && "PositronX[1]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
		break;
	default:
		cut = cX && cY && cR && c20 && cGamma && cGammaMax && cPe && cN;
	}
	tMC->Project(hMC->GetName(), str, cut);
	hMC->Sumw2();
	cv->Clear();
	hExp->SetMarkerStyle(kFullStar);
	hExp->SetMarkerColor(kBlue);
	hExp->SetLineColor(kBlue);
	hExp->SetMarkerSize(3);
	hExp->SetMarkerStyle(kDot);
	hMC->SetLineColor(kBlack);
	hMC->SetLineWidth(3);
	hMC->Scale(hExp->Integral() / hMC->Integral());
	hExp->Draw();
	hMC->Draw("same");
	TLegend *lg = new TLegend(0.5, 0.2, 0.65, 0.35);
	lg->AddEntry(hExp, "IBD", "P");
	lg->AddEntry(hMC,  "MC",  "LE");
	lg->Draw();
	cv->Update();
}

TChain *make_mc_tree(const char *fname = NULL)
{
	int i;
	const char files[4][256] = {
		"mc_ibd_235U_transcode_pair.root",
		"mc_ibd_238U_transcode_pair.root",
		"mc_ibd_239Pu_transcode_pair.root",
		"mc_ibd_241Pu_transcode_pair.root"
	};

	TChain *tMC = new TChain("DanssPair", "DanssPair");
	if (fname != NULL && strcasecmp(fname, "all")) {
		tMC->AddFile(fname);
	} else {
		for (i=0; i<4; i++) tMC->AddFile(files[i]);
	}
	if (!tMC->GetEntries()) return NULL;
	return tMC;
}

void mc2ibd(const char *mcfile, const char *bgndfile)
{
	int i;

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	TFile *fBgnd = new TFile(bgndfile);
	if (!fBgnd->IsOpen()) return;
	TChain *tMC = make_mc_tree(mcfile);
	if (!tMC) {
		printf("DanssPair tree not found in %s\n", mcfile);
		return;
	}
	TCanvas *cv = new TCanvas("CV", "CV", 1200, 900);
	cv->SaveAs("mc2ibd.pdf[");
	
	mc2ibd_gtDiff(tMC, fBgnd, cv);
	cv->SaveAs("mc2ibd.pdf");

	mc2ibd_R1(tMC, fBgnd, cv);
	cv->SaveAs("mc2ibd.pdf");

	mc2ibd_R2(tMC, fBgnd, cv);
	cv->SaveAs("mc2ibd.pdf");

	mc2ibd_NE(tMC, fBgnd, cv);
	cv->SaveAs("mc2ibd.pdf");

//	mc2ibd_PX(tMC, fBgnd, cv);
//	cv->SaveAs("mc2ibd.pdf");
	
	for(i=0; i<3; i++) {
		mc2ibd_PPX('X'+i, tMC, fBgnd, cv);
		cv->SaveAs("mc2ibd.pdf");
	}
	
	cv->SaveAs("mc2ibd.pdf]");

	fBgnd->Close();
}

void ibd_xyzratio(const char *fileA, const char *fileB, const char *title)
{
	TH1 *hA[3];
	TH1 *hB[3];
	TH1 *hR[3];
	int i;
	char str[256];
	const int Color[] = {kRed, kBlue, kGreen};
	
	TFile *fA = new TFile(fileA);
	TFile *fB = new TFile(fileB);
	
	for (i=0; i<3; i++) {
		sprintf(str, "hP%cA-diff", 'X'+i);
		hA[i] = (TH1 *) fA->Get(str);
		hB[i] = (TH1 *) fB->Get(str);
		if (!hA[i] || !hB[i]) return;
		sprintf(str, "hRatio%c", 'X'+i);
		if (i == 2) {
			hA[i]->Rebin(4);
			hB[i]->Rebin(4);
		}
		hR[i] = (TH1 *) hA[i]->Clone(str);
		hR[i]->Divide(hA[i], hB[i]);
		sprintf(str, "Ratio %s;cm;Ratio", title);
		hR[i]->SetTitle(str);
		hR[i]->SetStats(0);
		hR[i]->SetLineColor(Color[i]);
	}
	TCanvas *cv = new TCanvas("CV", "CV", 1200, 900);
	hR[0]->SetMinimum(0.9*hR[2]->GetMinimum());
	for (i=0; i<3; i++) hR[i]->Draw((i) ? "same" : "");
	TLegend *lg = new TLegend(0.12, 0.73, 0.2, 0.87);
	for (i=0; i<3; i++) {
		sprintf(str, "%c", 'X'+i);
		lg->AddEntry(hR[i], str, "LE");
	}
	lg->Draw();
	cv->SaveAs("XYZ-ratio.pdf");
	
	fA->Close();
	fB->Close();
}

void mc_xyzratio(const char *fname = NULL)
{
	TH1D *hA[3];
	TH1D *hB[3];
	TH1D *hR[3];
	int i;
	char str[256];
	const int Color[] = {kRed, kBlue, kGreen};
	TCut cut;

	TChain *tMC = make_mc_tree(fname);
	if (!tMC) {
		printf("DanssPair tree not found\n");
		return;
	}

	gROOT->cd();
	for (i=0; i<3; i++) {
		sprintf(str, "hMCgt3P%c", 'X'+i);
		hA[i] = new TH1D(str, "MC gt 3 MeV", 25, 0, 100);
		sprintf(str, "hMClt3P%c", 'X'+i);
		hB[i] = new TH1D(str, "MC lt 3 MeV", 25, 0, 100);
		sprintf(str, "PositronX[%d]+%4.1f", i, (i==2) ? 0.5 : 2.0);
		switch (i) {
		case 0:
			cut = cY && "PositronX[0]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
			break;
		case 1:
			cut = cX && "PositronX[1]>=0" && cZ && cR && c20 && cGamma && cGammaMax && cPe && cN;
			break;
		default:
			cut = cX && cY && cR && c20 && cGamma && cGammaMax && cPe && cN;
		}
		gROOT->cd();
		tMC->Project(hA[i]->GetName(), str, cut && "PositronEnergy > 3");
		tMC->Project(hB[i]->GetName(), str, cut && "PositronEnergy < 3");
		hA[i]->Sumw2();
		hB[i]->Sumw2();
		sprintf(str, "hRatioMCP%c", 'X'+i);
		hR[i] = (TH1D*) hA[i]->Clone(str);
		hR[i]->Divide(hA[i], hB[i]);
		hR[i]->SetTitle("MC Ratio (PE gt 3 MeV) / (PE lt 3 MeV);cm;Ratio");
		hR[i]->SetLineColor(Color[i]);
		hR[i]->SetStats(0);
	}
	TCanvas *cv = new TCanvas("CV", "CV", 1200, 900);
	for (i=0; i<3; i++) hR[i]->Draw((i) ? "same" : "");
	TLegend *lg = new TLegend(0.12, 0.73, 0.2, 0.87);
	for (i=0; i<3; i++) {
		sprintf(str, "%c", 'X'+i);
		lg->AddEntry(hR[i], str, "LE");
	}
	lg->Draw();
	cv->SaveAs("XYZ-ratioMC.pdf");
}
