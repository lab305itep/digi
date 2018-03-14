TFile *fData;

void change_file_suffix(char *to, int len, const char *from, const char *where, const char *what)
{
	char *ptr;
	
	strncpy(to, from, len - strlen(what));
	ptr = strstr(to, where);
	if (ptr) *ptr = '\0';
	strcat(to, what);
}

int sum_of_spectra(TH1D *hSum, const char *posmask, int permask, double bgScale = 1.0, double *days = NULL)
{
#include "positions.h"
	int N;
	TH1D *hSig;
	TH1D *hConst;
	TH1D *hBgnd;
	int i;
	char str[1024];
	double tSum;
	double dt;
	char *ptr;
	int Cnt;
//		Background tail correction (mHz)
//	TF1 fBgndN("fBgndN", "0.0163-0.0007*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.0722-0.0034*x", 0, 100);
//		based on oct16 - jan18 statistics
	TF1 fBgndN("fBgndN", "0.01426-0.000613*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07142-0.003486*x", 0, 100);

	N = sizeof(positions) / sizeof(positions[0]);
	hSum->Reset();
	tSum = 0;
	for (i=0; i<N; i++) {
		ptr = strchr(posmask, positions[i].name[0]);
		if (!(ptr && ((1 << positions[i].period) & permask))) continue;
		sprintf(str, "%s_hSig-diff", positions[i].name);
		hSig = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hConst", positions[i].name);
		hConst = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		hBgnd = (TH1D*) fData->Get(str);
		if (!(hSig && hConst && hBgnd)) continue;
		dt = hConst->GetBinContent(1) / 1000.0;	// seconds * 10^3
		tSum += dt;
		hSum->Add(hSig);
		hSum->Add(&fBgndN, -dt);
		hSum->Add(hBgnd, -positions[i].bgnd * bgScale);
		hSum->Add(&fBgndC, dt * positions[i].bgnd * bgScale);
	}
	
	if (days) *days = tSum / 86.4;
	if (tSum == 0) return 0;
	Cnt = hSum->Integral();
	hSum->Scale(86.4 / tSum);
	return Cnt;
}

int sum_of_spectral(TH1D *hSum, const char *posmask, int pfrom, int pto, double bgScale = 1.0, double *days = NULL)
{
#include "positions.h"
	int N;
	TH1D *hSig;
	TH1D *hConst;
	TH1D *hBgnd;
	int i;
	char str[1024];
	double tSum;
	double dt;
	char *ptr;
	int Cnt;
//		Background tail correction (mHz)
//	TF1 fBgndN("fBgndN", "0.0163-0.0007*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.0722-0.0034*x", 0, 100);
//		based on oct16 - jan18 statistics
	TF1 fBgndN("fBgndN", "0.01426-0.000613*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07142-0.003486*x", 0, 100);

	N = sizeof(positions) / sizeof(positions[0]);
	hSum->Reset();
	tSum = 0;
	for (i=pfrom; i<pto && i<N; i++) {
		ptr = strchr(posmask, positions[i].name[0]);
		if (!(ptr && positions[i].period)) continue;
		sprintf(str, "%s_hSig-diff", positions[i].name);
		hSig = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hConst", positions[i].name);
		hConst = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		hBgnd = (TH1D*) fData->Get(str);
		if (!(hSig && hConst && hBgnd)) continue;
		dt = hConst->GetBinContent(1) / 1000.0;	// seconds * 10^3
		tSum += dt;
		hSum->Add(hSig);
		hSum->Add(&fBgndN, -dt);
		hSum->Add(hBgnd, -positions[i].bgnd * bgScale);
		hSum->Add(&fBgndC, dt * positions[i].bgnd * bgScale);
	}
	
	if (days) *days = tSum / 86.4;
	if (tSum == 0) return 0;
	Cnt = hSum->Integral();
	hSum->Scale(86.4 / tSum);
	return Cnt;
}

void sum_of_raw(TH1D *hSumSig, TH1D *hSumBgnd, const char *posmask, int permask)
{
#include "positions.h"
	int N;
	TH1D *hSig;
	TH1D *hConst;
	TH1D *hBgnd;
	int i;
	char str[1024];
	double tSum;
	double dt;
	char *ptr;

	N = sizeof(positions) / sizeof(positions[0]);
	hSumSig->Reset();
	hSumBgnd->Reset();
	tSum = 0;
	for (i=0; i<N; i++) {
		ptr = strchr(posmask, positions[i].name[0]);
		if (!(ptr && ((1 << positions[i].period) & permask))) continue;
		sprintf(str, "%s_hSig-diff", positions[i].name);
		hSig = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hConst", positions[i].name);
		hConst = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		hBgnd = (TH1D*) fData->Get(str);
		if (!(hSig && hConst && hBgnd)) continue;
		dt = hConst->GetBinContent(1) / 1000.0;	// seconds * 10^3
		tSum += dt;
		hSumSig->Add(hSig);
		hSumBgnd->Add(hBgnd);
	}
	
	if (tSum == 0) return;
	hSumSig->Scale(86.4 / tSum);
	hSumBgnd->Scale(86.4 / tSum);
}

void draw_spectra_page(TCanvas *cv, const char *title, int periodmask, double bgScale)
{
	TLegend *lg;
	char strs[128];
	char strl[1024];
	double val, err;
	int Cnt, n;
	TLatex txt;
	const double OtherBlockFraction = 0.0060;	// Distances to other reactors: 160, 336 and 478 m

	cv->Clear();
	lg = new TLegend(0.35, 0.65, 0.9, 0.9);
	sprintf(strs, "hUp_%d", periodmask);
	sprintf(strl, "Positron spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUp = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hUp, "u", periodmask, bgScale);
	Cnt = n;
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -OtherBlockFraction);
	hUp->Write();
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.08);
	hUp->SetTitle(title);
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, hUp->FindBin(7.999), err);
	sprintf(strs, "  Up: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hUp, strs, "l");

	sprintf(strs, "hMid_%d", periodmask);
	sprintf(strl, "Positron spectrum %s, MIDDLE;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hMid = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hMid, "m", periodmask, bgScale);
	Cnt += n;
	hMid->Add(hTmp, -OtherBlockFraction);
	hMid->Write();
	hMid->SetLineColor(kGreen);
	hMid->SetFillColor(kGreen-10);
	hMid->Draw("same,hist,e");
	val = hMid->IntegralAndError(1, hMid->FindBin(7.999), err);
	sprintf(strs, " Mid: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hMid, strs, "l");

	sprintf(strs, "hDown_%d", periodmask);
	sprintf(strl, "Positron spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDown = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hDown, "d", periodmask, bgScale);
	Cnt += n;
	hDown->Add(hTmp, -OtherBlockFraction);
	hDown->Write();
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, hDown->FindBin(7.999), err);
	sprintf(strs, "Down: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hDown, strs, "l");

	hUp->Draw("axis,same");
	lg->Draw();
	sprintf(strs, "%d events", Cnt);
	txt.DrawLatexNDC(0.6, 0.5, strs);
	cv->Update();
	delete hTmp;
}

void draw_spectra_pagel(TCanvas *cv, const char *title, int pfrom, int pto, double bgScale)
{
	TLegend *lg;
	char strs[128];
	char strl[1024];
	double val, err;
	int Cnt, n;
	TLatex txt;
	const double OtherBlockFraction = 0.0060;	// Distances to other reactors: 160, 336 and 478 m

	cv->Clear();
	lg = new TLegend(0.35, 0.65, 0.9, 0.9);
	sprintf(strs, "hUp_%d_%d", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUp = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hUp, "u", pfrom, pto, bgScale);
	Cnt = n;
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -OtherBlockFraction);
	hUp->Write();
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.08);
	hUp->SetTitle(title);
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, hUp->FindBin(7.999), err);
	sprintf(strs, "  Up: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hUp, strs, "l");

	sprintf(strs, "hMid_%d_%d", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, MIDDLE;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hMid = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hMid, "m", pfrom, pto, bgScale);
	Cnt += n;
	hMid->Add(hTmp, -OtherBlockFraction);
	hMid->Write();
	hMid->SetLineColor(kGreen);
	hMid->SetFillColor(kGreen-10);
	hMid->Draw("same,hist,e");
	val = hMid->IntegralAndError(1, hMid->FindBin(7.999), err);
	sprintf(strs, " Mid: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hMid, strs, "l");

	sprintf(strs, "hDown_%d_%d", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDown = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hDown, "d", pfrom, pto, bgScale);
	Cnt += n;
	hDown->Add(hTmp, -OtherBlockFraction);
	hDown->Write();
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, hDown->FindBin(7.999), err);
	sprintf(strs, "Down: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hDown, strs, "l");

	hUp->Draw("axis,same");
	lg->Draw();
	sprintf(strs, "%d events", Cnt);
	txt.DrawLatexNDC(0.6, 0.5, strs);
	cv->Update();
	delete hTmp;
}

void draw_tail_hist(const char *title, const char *posmask)
{
	TLegend *lg;
	char strs[128];
	char strl[1024];
	TPaveStats *pv;
	double y1, y2;

	lg = new TLegend(0.35, 0.65, 0.60, 0.9);
	sprintf(strs, "hTailN_%s", posmask);
	sprintf(strl, "Positron raw spectrum %s;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hN = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hTailC_%s", posmask);
	sprintf(strl, "Muon raw spectrum %s;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hC = new TH1D(strs, strl, 60, 1, 16);
	sum_of_raw(hN, hC, posmask, 30);
	hN->Write();
	hC->Write();
	hC->SetLineColor(kRed);
	hC->GetYaxis()->SetLabelSize(0.08);
	hC->SetTitle(title);
	hC->SetTitleSize(0.08);
	hC->GetXaxis()->SetRange(hC->FindBin(8.001), hC->FindBin(15.999));
	hC->SetMinimum(0);
	hC->Fit("pol1", "", "", 10, 16);
	hN->SetLineColor(kBlue);
	hN->Fit("pol1", "", "sames", 10, 16);
	gPad->Update();
	pv = (TPaveStats *)hC->FindObject("stats");
	pv->SetLineColor(kRed);
	pv->SetTextColor(kRed);
	y1 = pv->GetY1NDC();
	y2 = pv->GetY2NDC();
	pv->SetY1NDC(2 * y1 - y2);
	pv->SetY2NDC(y1);
	pv->Draw();
	pv = (TPaveStats *)hN->FindObject("stats");
	pv->SetLineColor(kBlue);
	pv->SetTextColor(kBlue);
	pv->Draw();
	lg->AddEntry(hN, "Neutrino", "l");
	lg->AddEntry(hC, "Cosmic", "l");
	lg->Draw();
}

void draw_single_ratio(const char *nameA, const char *nameB, const char *name, const char *title, double min=0.6, double max=1.2, int last=60)
{
	TH1D *hA = (TH1D *) gROOT->FindObject(nameA);
	TH1D *hB = (TH1D *) gROOT->FindObject(nameB);
	if (!hA || !hB) {
		printf("Can not find hist: %s or/and %s. Step %s\n", nameA, nameB, title);
		return;
	}
	TH1D *hAB = (TH1D*) hA->Clone(name);
	hAB->SetName(name);
	hAB->SetTitle(title);
	hAB->Divide(hA, hB);
	hAB->SetLineColor(kBlue);
	hAB->Write();	
	hAB->SetMinimum(min);
	hAB->SetMaximum(max);
	hAB->GetYaxis()->SetLabelSize(0.08);
	hAB->GetXaxis()->SetRange(1, last);
	hAB->Fit("pol0", "", "", 1, 8);
}

void draw_normalized_ratio(const int maskA, const int maskB, const char *name, const char *title, double min=0.6, double max=1.2, int last=60, double bgScale = 1.0)
{
	char strs[128];
	double daysAU, daysAM, daysAD, daysBU, daysBM, daysBD;
//	int i;
//	const char posmask[3][3] = {"u", "m", "d"};
	const double Nfactor[3] = {1.0, (11.7*11.7)/(10.7*10.7), (12.7*12.7)/(10.7*10.7)};
	const double OtherBlockFraction = 0.0060;	// Distances to other reactors: 160, 336 and 478 m
//		Book
	TH1D *hAU = new TH1D("hNRAU", "", 60, 1, 16);
	TH1D *hAM = new TH1D("hNRAM", "", 60, 1, 16);
	TH1D *hAD = new TH1D("hNRAD", "", 60, 1, 16);
	TH1D *hBU = new TH1D("hNRBU", "", 60, 1, 16);
	TH1D *hBM = new TH1D("hNRBM", "", 60, 1, 16);
	TH1D *hBD = new TH1D("hNRBD", "", 60, 1, 16);
	TH1D *hAB = new TH1D(name, title, 60, 1, 16);
//		Make initial sum
	sum_of_spectra(hAU, "u", maskA, bgScale, &daysAU);
	sum_of_spectra(hAM, "m", maskA, bgScale, &daysAM);
	sum_of_spectra(hAD, "d", maskA, bgScale, &daysAD);
	sum_of_spectra(hBU, "u", maskB, bgScale, &daysBU);
	sum_of_spectra(hBM, "m", maskB, bgScale, &daysBM);
	sum_of_spectra(hBD, "d", maskB, bgScale, &daysBD);
//		Subtract other blocks
	TH1D *hTmpA = (TH1D *) hAU->Clone("hTmpA");
	TH1D *hTmpB = (TH1D *) hBU->Clone("hTmpB");
	hAU->Add(hTmpA, -OtherBlockFraction);
	hAM->Add(hTmpA, -OtherBlockFraction);
	hAD->Add(hTmpA, -OtherBlockFraction);
	hBU->Add(hTmpB, -OtherBlockFraction);
	hBM->Add(hTmpB, -OtherBlockFraction);
	hBD->Add(hTmpB, -OtherBlockFraction);
//		Add spectra with R^2 weights
	hTmpA->Reset();
	hTmpA->Add(hAU, Nfactor[0]*daysAU);
	hTmpA->Add(hAM, Nfactor[1]*daysAM);
	hTmpA->Add(hAD, Nfactor[2]*daysAD);
	hTmpA->Scale(1.0 / (daysAU + daysAM + daysAD));
	hTmpB->Reset();
	hTmpB->Add(hBU, Nfactor[0]*daysBU);
	hTmpB->Add(hBM, Nfactor[1]*daysBM);
	hTmpB->Add(hBD, Nfactor[2]*daysBD);
	hTmpB->Scale(1.0 / (daysBU + daysBM + daysBD));
//		Draw
	hAB->Divide(hTmpA, hTmpB);
	hAB->SetLineColor(kBlue);
	hAB->Write();	
	hAB->SetMinimum(min);
	hAB->SetMaximum(max);
	hAB->GetYaxis()->SetLabelSize(0.08);
	hAB->GetXaxis()->SetRange(1, last);
	hAB->Fit("pol1", "", "", 1, 8);
//		Clean
	delete hAU;
	delete hAM;
	delete hAD;
	delete hBU;
	delete hBD;
	delete hBM;
	delete hTmpA;
	delete hTmpB;
}

void draw_period_ratios(int m1, int m2, const char *title)
{
	char nameA[64];
	char nameB[64];
	char nameR[64];
	TH1D *hRUp;
	TH1D *hRMid;
	TH1D *hRDown;
	TH1D *hA;
	TH1D *hB;
	
	sprintf(nameA, "hUp_%d", 1 << m1);
	sprintf(nameB, "hUp_%d", 1 << m2);
	hA = (TH1D *) gROOT->FindObject(nameA);
	hB = (TH1D *) gROOT->FindObject(nameB);
	if (!hA || !hB) {
		printf("Can not find hist: %s or/and %s. Step %s\n", nameA, nameB, title);
		return;
	}
	sprintf(nameR, "hRUp_%d_%d", 1 << m1, 1 << m2);
	hRUp = (TH1D*) hA->Clone(nameR);
	hRUp->SetTitle(title);
	hRUp->Divide(hA, hB);
	hRUp->SetLineColor(kRed);
	hRUp->SetMarkerColor(kRed);
	hRUp->SetMarkerSize(2);
	hRUp->SetMarkerStyle(kFullTriangleUp);
	
	sprintf(nameA, "hMid_%d", 1 << m1);
	sprintf(nameB, "hMid_%d", 1 << m2);
	hA = (TH1D *) gROOT->FindObject(nameA);
	hB = (TH1D *) gROOT->FindObject(nameB);
	if (!hA || !hB) {
		printf("Can not find hist: %s or/and %s. Step %s\n", nameA, nameB, title);
		return;
	}
	sprintf(nameR, "hRMid_%d_%d", 1 << m1, 1 << m2);
	hRMid = (TH1D*) hA->Clone(nameR);
	hRMid->SetTitle(title);
	hRMid->Divide(hA, hB);
	hRMid->SetLineColor(kGreen);
	hRMid->SetMarkerColor(kGreen);
	hRMid->SetMarkerSize(2);
	hRMid->SetMarkerStyle(kFullCircle);
	
	sprintf(nameA, "hDown_%d", 1 << m1);
	sprintf(nameB, "hDown_%d", 1 << m2);
	hA = (TH1D *) gROOT->FindObject(nameA);
	hB = (TH1D *) gROOT->FindObject(nameB);
	if (!hA || !hB) {
		printf("Can not find hist: %s or/and %s. Step %s\n", nameA, nameB, title);
		return;
	}
	sprintf(nameR, "hRDown_%d_%d", 1 << m1, 1 << m2);
	hRDown = (TH1D*) hA->Clone(nameR);
	hRDown->SetTitle(title);
	hRDown->Divide(hA, hB);
	hRDown->SetLineColor(kBlue);
	hRDown->SetMarkerColor(kBlue);
	hRDown->SetMarkerSize(2);
	hRDown->SetMarkerStyle(kFullTriangleDown);
	
	hRUp->SetMinimum(0.9);
	hRUp->SetMaximum(1.4);
	hRUp->GetXaxis()->SetRange(1, 28);
	hRUp->DrawCopy();
	hRMid->DrawCopy("same");
	hRDown->DrawCopy("same");
	TLegend *lg = new TLegend(0.3, 0.7, 0.5, 0.85);
	lg->AddEntry(hRUp, "Up", "LP");
	lg->AddEntry(hRMid, "Middle", "LP");
	lg->AddEntry(hRDown, "down", "LP");
	
	delete hRUp;
	delete hRMid;
	delete hRDown;
}

void danss_calc_ratio2w(const char *fname, double bgScale = 5.6/2.5)
{
	TCanvas *cv;
	TFile *f;
	TFile *fOut;
	char str[1024];
	TLatex *txt;
	double val, err;
	char pname[1024];
	char rname[2014];
	TLegend *lg;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.08);
	gStyle->SetTitleYSize(0.08);
	gStyle->SetLabelSize(0.08);
	gStyle->SetPadLeftMargin(0.2);
	gStyle->SetPadBottomMargin(0.2);
//	gStyle->SetLineWidth(2);
	
	cv = new TCanvas("CV", "Results", 2400, 1200);
	fData = new TFile(fname);
	if (!fData->IsOpen()) return;
	change_file_suffix(pname, sizeof(pname), fname, ".root", "-calc.pdf");
	change_file_suffix(rname, sizeof(pname), fname, ".root", "-calc.root");
	sprintf(str, "%s[", pname);
	cv->Print(str);
	fOut = new TFile(rname, "RECREATE");
	fOut->cd();
	txt = new TLatex();
	
//	Page 1:	Spectra All
	draw_spectra_page(cv, "Oct 16-Jan 18", 0x1E, bgScale);
//	cv->Print(pname);
	
//	Page 1a:	Spectra All but April-June 16.
	draw_spectra_page(cv, "Oct 16-Jan 18", 0x1C, bgScale);
	cv->Print(pname);

//	Page 1b:	before Sep17
	draw_spectra_pagel(cv, "Oct 16-Sep 17", 0, 143, bgScale);
	cv->Print(pname);

//	Page 2:	Spectra April-June
//	draw_spectra_page(cv, "April-June 16", 2, bgScale);
//	cv->Print(pname);

//	Page 3:	Spectra October-November
	draw_spectra_page(cv, "Oct 16 - Feb 17", 4, bgScale);
	cv->Print(pname);

//	Page 4:	Spectra December-January
	draw_spectra_page(cv, "March-July 17", 8, bgScale);
	cv->Print(pname);

//	Page 5:	Spectra February-March
	draw_spectra_page(cv, "Aug. 17-Jan. 18", 0x10, bgScale);
	cv->Print(pname);

//	Some more spectra pages
	draw_spectra_page(cv, "Oct 16 - Jul 17", 12, bgScale);
	cv->Print(pname);

	draw_spectra_page(cv, "Mar 17 - Jan 18", 24, bgScale);
	cv->Print(pname);

	draw_spectra_pagel(cv, "Aug 17 - Jan 18", 144, 193, bgScale);
	cv->Print(pname);

//	Page 6: Total Spectrum (all positions)
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	TH1D *hSum = new TH1D("hSum", "Positron spectrum Oct 16 - Jan 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum, "umdrs", 0x1E, bgScale);
	hSum->Write();
	hSum->Draw("e");

	cv->cd(2);
	TH1D *hSum1 = new TH1D("hSum1", "Positron spectrum Oct 16 - Sep 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectral(hSum1, "umdrs", 0, 143, bgScale);
	hSum1->Write();
	hSum1->Draw("e");

	cv->cd(3);
	TH1D *hSum2 = new TH1D("hSum2", "Positron spectrum Oct 16-Feb 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum2, "umdrs", 4, bgScale);
	hSum2->Write();
	hSum2->Draw("e");

	cv->cd(4);
	TH1D *hSum3 = new TH1D("hSum3", "Positron spectrum March-July 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum3, "umdrs", 8, bgScale);
	hSum3->Write();
	hSum3->Draw("e");

	cv->cd(5);
	TH1D *hSum4 = new TH1D("hSum4", "Positron spectrum Aug 17 - Jan 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum4, "umdrs", 0x10, bgScale);
	hSum4->Write();
	hSum4->Draw("e");

	cv->cd(6);
	TH1D *hSum23 = new TH1D("hSum23", "Positron spectrum Oct 16 - July 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	hSum23->Add(hSum2, hSum3);
	hSum23->Scale(0.5);
	hSum23->Write();
//	TH1D *hSum34 = new TH1D("hSum34", "Positron spectrum March-Sep 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
//	hSum34->Add(hSum3, hSum4);
//	hSum34->Write();
	hSum23->SetLineColor(kRed);
	hSum4->SetLineColor(kBlue);
	lg = new TLegend(0.6, 0.7, 0.9, 0.8);
	lg->AddEntry(hSum23, "Oct 16 - July 17", "LE");
	lg->AddEntry(hSum4, "Aug 17 - Jan 18", "LE");
	hSum23->Draw("e");
	hSum4->Draw("e,same");
	lg->Draw();

	cv->Update();
	cv->Print(pname);

//	Page 7:	Down/Up
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hDown_4", "hUp_4", "hDownUp_4", "Ratio Down/Up Oct 16-Feb 17;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(2);
	draw_single_ratio("hDown_12", "hUp_12", "hDownUp_12", "Ratio Down/Up Oct 16-Jul 17;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(3);
	draw_single_ratio("hDown_0_143", "hUp_0_143", "hDownUp_0_143", "Ratio Down/Up Oct 16-Sep 17;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(4);
	draw_single_ratio("hDown_28", "hUp_28", "hDownUp_28", "Ratio Down/Up Oct 16-Jan 18;Positron energy, MeV", 0.6, 0.9, 28);

	cv->Update();
	cv->Print(pname);
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hDown_8", "hUp_8", "hDownUp_8", "Ratio Down/Up Feb-July 17;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(2);
	draw_single_ratio("hDown_16", "hUp_16", "hDownUp_16", "Ratio Down/Up Aug 17-Jan 18;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(3);
	draw_single_ratio("hDown_24", "hUp_24", "hDownUp_24", "Ratio Down/Up Mar 17-Jan 18;Positron energy, MeV", 0.6, 0.9, 28);

	cv->cd(4);
	draw_single_ratio("hDown_144_193", "hUp_144_193", "hDownUp_144_193", "Ratio Down/Up Oct 17-Jan 18;Positron energy, MeV", 0.6, 0.9, 28);

	cv->Update();
	cv->Print(pname);

//	Page 7a: Main plot
	cv->Clear();
	draw_single_ratio("hDown_30", "hUp_30", "hDownUpAll", "Ratio Down/Up Oct 16-Jan 18 (ALL);Positron energy, MeV;#frac{N_{DOWN}}{N_{UP}}", 0.6, 0.9, 28);
	cv->Update();
	cv->Print(pname);

//	Page 8:	Mid/Up
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hMid_28", "hUp_28", "hMidUp0", "Ratio Middle/Up All but Apr-June 16;Positron energy, MeV");
	draw_single_ratio("hMid_30", "hUp_30", "hMidUpAll", "Ratio Middle/Up All;Positron energy, MeV", 0.7, 1, 28);

	cv->cd(2);
	draw_single_ratio("hMid_0_143", "hUp_0_143", "hMidUp1", "Ratio Middle/Up Oct 16-Sep 17;Positron energy, MeV", 0.7, 1, 28);

	cv->cd(3);
	draw_single_ratio("hMid_4", "hUp_4", "hMidUp2", "Ratio Middle/Up Oct 16 - Feb 17;Positron energy, MeV", 0.7, 1, 28);

	cv->cd(4);
	draw_single_ratio("hMid_8", "hUp_8", "hMidUp3", "Ratio Middle/Up March-July 17;Positron energy, MeV", 0.7, 1, 28);

	cv->cd(5);
	draw_single_ratio("hMid_16", "hUp_16", "hMidUp4", "Ratio Middle/Up Aug 17-Jan 18;Positron energy, MeV", 0.7, 1, 28);

	cv->Update();
	cv->Print(pname);

//	Page 8.1: Double ratios Down/Up between periods.
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hDownUp_4", "hDownUp_28", "hDownUpdr1a", "Double ratio Down/Up Oct 16-Feb 17 / All;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(2);
	draw_single_ratio("hDownUp_12", "hDownUp_28", "hDownUpdr2a", "Double ratio Down/Up Oct 16-Jul 17 / All;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(3);
	draw_single_ratio("hDownUp_0_143", "hDownUp_28", "hDownUpdr3a", "Double ratio Down/Up Oct 16-Sep 17 / All;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(4);
	draw_single_ratio("hDownUp_4", "hDownUp_24", "hDownUpdr1r", "Double ratio Down/Up Oct 16-Feb 17 / Mar 17-Jan 18;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(5);
	draw_single_ratio("hDownUp_12", "hDownUp_16", "hDownUpdr2r", "Double ratio Down/Up Oct 16-Jul 17 / Aug 17-Jan 18;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(6);
	draw_single_ratio("hDownUp_0_143", "hDownUp_144_193", "hDownUpdr3r", "Double ratio Down/Up Oct 16-Sep 17 / Oct 17-Jan 18;Positron energy, MeV", 0.8, 1.2, 28);

	cv->Update();
	cv->Print(pname);
//	Page 9: Period ratios
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hSum4", "hSum23", "hRatio4_23", "Unnormalized Ratio Aug 17-Jan 18/Oct 16-Feb 17;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hSum3", "hSum2", "hRatio3_2", "Unnormalized Ratio March-July 17/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hSum4", "hSum2", "hRatio4_2", "Unnormalized Ratio Aug 17-Jan 18/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(4);
//	draw_single_ratio("hSum2", "hSum1", "hRatio2_1", "Unnormalized Ratio Oct 16 - Feb 17/April-June 16;Positron energy, MeV");

//	cv->cd(5);
//	draw_single_ratio("hSum3", "hSum1", "hRatio3_1", "Unnormalized Ratio March-July 17/April-June 16;Positron energy, MeV");

//	cv->cd(6);
//	draw_single_ratio("hSum4", "hSum1", "hRatio4_1", "Unnormalized Ratio Aug-Sep 17/April-June 16;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 9a: after shutdown / before shutdown
	cv->Clear();
	draw_normalized_ratio(16, 8, "hNRatio43", "Normalized ratio after shutdown / before shutdown;Positron energy, MeV", 0.9, 1.4, 28, bgScale);
	cv->Update();
	cv->Print(pname);

//	Page 9b: after shutdown / before shutdown
	cv->Clear();
	draw_period_ratios(4, 3, "Ratio after shutdown / before shutdown;Positron energy, MeV");
	cv->Update();
	cv->Print(pname);

//	Page 10: Without background subtraction
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_tail_hist("all", "umd");
	
	cv->cd(2);
	draw_tail_hist("up", "u");

	cv->cd(3);
	draw_tail_hist("middle", "m");

	cv->cd(4);
	draw_tail_hist("down", "d");
	
	cv->Update();
	cv->Print(pname);

	sprintf(str, "%s]", pname);
	cv->Print(str);
	fData->Close();
	delete cv;
	fOut->Close();
}

void danss_draw_sp2(const char *fname, double bgScale = 5.6/2.5)
{
	TCanvas *cv;
	TFile *f;
	TFile *fOut;
	char str[1024];
	TLatex *txt;
	double val, err;
	char pname[1024];
	TLegend *lg;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
//	gStyle->SetLineWidth(2);
	
	cv = new TCanvas("CV", "Results", 1200, 900);
	fData = new TFile(fname);
	if (!fData->IsOpen()) return;
	change_file_suffix(pname, sizeof(pname), fname, ".root", "-sp2.png");
	txt = new TLatex();
	
	lg = new TLegend(0.45, 0.75, 0.9, 0.9);
	TH1D *hUp = new TH1D("hUp", "Positron spectrum All, UP;Positron energy, MeV;Events/(day*0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hUp, "u", 30, bgScale);
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -0.0045);
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.06);
	hUp->SetTitle("");
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, hUp->FindBin(7.999), err);
	sprintf(str, "Up:     %5.0f #pm %2.0f / day", val, err);
	lg->AddEntry(hUp, str, "l");

	TH1D *hDown = new TH1D("hDown", "Positron spectrum All, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hDown, "d", 30, bgScale);
	hDown->Add(hTmp, -0.0045);
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, hDown->FindBin(7.999), err);
	sprintf(str, "Down: %5.0f #pm %2.0f / day", val, err);
	lg->AddEntry(hDown, str, "l");
	
	TFile *fMc = new TFile("mc_fuel_middle.root");
	TH1D *hMc = (TH1D*) fMc->Get("hMcMiddleMixt");
	if (hMc) {
		hMc->SetMarkerStyle(kFullStar);
		hMc->SetMarkerColor(kGreen+2);
		hMc->SetMarkerSize(2);
		hMc->Scale(hUp->Integral(1, hUp->FindBin(7.999)) / hMc->Integral(1, hMc->FindBin(7.999)));
		hMc->Draw("same");
		lg->AddEntry(hMc, "Monte Carlo", "p");
	}

	hUp->Draw("axis,same");
	lg->Draw();
	cv->Print(pname);
}
