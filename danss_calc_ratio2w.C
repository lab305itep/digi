TFile *fData;

void change_file_suffix(char *to, int len, const char *from, const char *where, const char *what)
{
	char *ptr;
	
	strncpy(to, from, len - strlen(what));
	ptr = strstr(to, where);
	if (ptr) *ptr = '\0';
	strcat(to, what);
}

//	Makes an integer out of date.
//	A string dd.mm.yy expected.
//	result = dd + 31*mm + 366*yy
//	We don't care if numbers are not sequential
int date2int(const char *str)
{
	int d, m, y;
	d = 10 * (str[0] - '0') + str[1] - '0';
	m = 10 * (str[3] - '0') + str[4] - '0';
	y = 10 * (str[6] - '0') + str[7] - '0';	
	return d + 31*m + 366*y;
}

// Check if name like "xxxx_dd.mm.yy" is in the range from-to
int is_in_date_range(const char *name, const char *from, const char *to)
{
	int iname, ifrom, ito;
	const char *ptr;
	
	if (!(from && to)) return true;
	ifrom = date2int(from);
	ito = date2int(to);
	ptr = strrchr(name, '_');
	if (!ptr) {
		ptr = name;
	} else {
		ptr++;
	}
	iname = date2int(ptr);
	return (iname >= ifrom && iname < ito);
}

int sum_of_spectra(TH1D *hSum, TH1D *hSub, const char *posmask, int permask, double bgScale = 1.0, double *days = NULL)
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
//	TF1 fBgndN("fBgndN", "0.01426-0.000613*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.07142-0.003486*x", 0, 100);
//		based on oct16 - mar18 statistics
	TF1 fBgndN("fBgndN", "0.01446-0.000621*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07252-0.003560*x", 0, 100);

	N = sizeof(positions) / sizeof(positions[0]);
	hSum->Reset();
	hSub->Reset();
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
		hSub->Add(hBgnd, positions[i].bgnd * bgScale);
		hSub->Add(&fBgndN, dt);
		hSub->Add(&fBgndC, -dt * positions[i].bgnd * bgScale);
	}
	hSum->Add(hSub, -1);
	if (days) *days = tSum / 86.4;
	if (tSum == 0) return 0;
	Cnt = hSum->Integral();
	hSum->Scale(86.4 / tSum);
	hSub->Scale(86.4 / tSum);
	return Cnt;
}

int sum_of_spectral(TH1D *hSum, TH1D *hSub, const char *posmask, const char *pfrom, const char *pto, double bgScale = 1.0, double *days = NULL)
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
//	TF1 fBgndN("fBgndN", "0.01426-0.000613*x", 0, 100);
//	TF1 fBgndC("fBgndC", "0.07142-0.003486*x", 0, 100);
//		based on oct16 - mar18 statistics
	TF1 fBgndN("fBgndN", "0.01446-0.000621*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07252-0.003560*x", 0, 100);

	N = sizeof(positions) / sizeof(positions[0]);
	hSum->Reset();
	hSub->Reset();
	tSum = 0;
	for (i=0; i<N; i++) if (is_in_date_range(positions[i].name, pfrom, pto) && positions[i].period) {	// ignore reactor not @100%
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
		hSub->Add(hBgnd, positions[i].bgnd * bgScale);
		hSub->Add(&fBgndN, dt);
		hSub->Add(&fBgndC, -dt * positions[i].bgnd * bgScale);
	}
	hSum->Add(hSub, -1);
	if (days) *days = tSum / 86.4;
	if (tSum == 0) return 0;
	Cnt = hSum->Integral();
	hSum->Scale(86.4 / tSum);
	hSub->Scale(86.4 / tSum);
	return Cnt;
}

void sum_of_raw(TH1D *hSumSig, TH1D *hSumBgnd, const char *posmask, int permask, const char *pfrom = NULL, const char *pto = NULL)
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
		if (!(is_in_date_range(positions[i].name, pfrom, pto) && ptr && ((1 << positions[i].period) & permask))) continue;
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
	sprintf(strs, "hUpSub_%d", periodmask);
	sprintf(strl, "Background spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUpSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hUp, hUpSub, "u", periodmask, bgScale);
	Cnt = n;
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -OtherBlockFraction);
	hUp->Write();
	hUpSub->Write();
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
	sprintf(strs, "hMidSub_%d", periodmask);
	sprintf(strl, "Background spectrum %s, MIDDLE;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hMidSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hMid, hMidSub, "m", periodmask, bgScale);
	Cnt += n;
	hMid->Add(hTmp, -OtherBlockFraction);
	hMid->Write();
	hMidSub->Write();
	hMid->SetLineColor(kGreen);
	hMid->SetFillColor(kGreen-10);
	hMid->Draw("same,hist,e");
	val = hMid->IntegralAndError(1, hMid->FindBin(7.999), err);
	sprintf(strs, " Mid: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hMid, strs, "l");

	sprintf(strs, "hDown_%d", periodmask);
	sprintf(strl, "Positron spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDown = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hDownSub_%d", periodmask);
	sprintf(strl, "Background spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDownSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hDown, hDownSub, "d", periodmask, bgScale);
	Cnt += n;
	hDown->Add(hTmp, -OtherBlockFraction);
	hDown->Write();
	hDownSub->Write();
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, hDown->FindBin(7.999), err);
	sprintf(strs, "Down: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hDown, strs, "l");

	hUpSub->SetLineColor(kBlack);
	hUpSub->SetFillColor(kGray);
	hUpSub->Draw("same,hist,e");
	lg->AddEntry(hUpSub, "Background for Up position", "l");

	hUp->Draw("axis,same");
	lg->Draw();
	sprintf(strs, "%d events", Cnt);
	txt.DrawLatexNDC(0.6, 0.5, strs);
	cv->Update();
	delete hTmp;
}

void draw_spectra_pagel(TCanvas *cv, const char *title, const char *pfrom, const char *pto, double bgScale)
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
	sprintf(strs, "hUp_%s_%s", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUp = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hUpSub_%s_%s", pfrom, pto);
	sprintf(strl, "Background spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUpSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hUp, hUpSub, "u", pfrom, pto, bgScale);
	Cnt = n;
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -OtherBlockFraction);
	hUp->Write();
	hUpSub->Write();
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.08);
	hUp->SetTitle(title);
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, hUp->FindBin(7.999), err);
	sprintf(strs, "  Up: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hUp, strs, "l");

	sprintf(strs, "hMid_%s_%s", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, MIDDLE;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hMid = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hMidSub_%s_%s", pfrom, pto);
	sprintf(strl, "Background spectrum %s, MIDDLE;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hMidSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hMid, hMidSub, "m", pfrom, pto, bgScale);
	Cnt += n;
	hMid->Add(hTmp, -OtherBlockFraction);
	hMid->Write();
	hMidSub->Write();
	hMid->SetLineColor(kGreen);
	hMid->SetFillColor(kGreen-10);
	hMid->Draw("same,hist,e");
	val = hMid->IntegralAndError(1, hMid->FindBin(7.999), err);
	sprintf(strs, " Mid: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hMid, strs, "l");

	sprintf(strs, "hDown_%s_%s", pfrom, pto);
	sprintf(strl, "Positron spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDown = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hDownSub_%s_%s", pfrom, pto);
	sprintf(strl, "Background spectrum %s, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hDownSub = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectral(hDown, hDownSub, "d", pfrom, pto, bgScale);
	Cnt += n;
	hDown->Add(hTmp, -OtherBlockFraction);
	hDown->Write();
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, hDown->FindBin(7.999), err);
	sprintf(strs, "Down: %d events %5.0f #pm%4.0f / day", n, val, err);
	lg->AddEntry(hDown, strs, "l");

	hUpSub->SetLineColor(kBlack);
	hUpSub->SetFillColor(kGray);
	hUpSub->Draw("same,hist,e");
	lg->AddEntry(hUpSub, "Background for Up position", "l");

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
	TLatex txt;

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
	sprintf(strs, "Integral 1-7 MeV = %5.1f events/day", hN->GetFunction("pol1")->Integral(hN->FindBin(1.001), hN->FindBin(6.999)));
	txt.SetTextColor(kBlue);
	txt.DrawLatex(11, 5, strs);
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
	hAB->GetYaxis()->SetTitle("");
	hAB->Fit("pol0", "", "", 1, 8);
}

void draw_normalized_ratio(const char* beginA, const char *endA, const char* beginB, const char *endB, 
	const char *name, const char *title, double min=0.6, double max=1.2, int last=60, double bgScale = 1.0)
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
	TH1D *hAUb = new TH1D("hNRAUb", "", 60, 1, 16);
	TH1D *hAMb = new TH1D("hNRAMb", "", 60, 1, 16);
	TH1D *hADb = new TH1D("hNRADb", "", 60, 1, 16);
	TH1D *hBUb = new TH1D("hNRBUb", "", 60, 1, 16);
	TH1D *hBMb = new TH1D("hNRBMb", "", 60, 1, 16);
	TH1D *hBDb = new TH1D("hNRBDb", "", 60, 1, 16);
	TH1D *hAB = new TH1D(name, title, 60, 1, 16);
//		Make initial sum
	sum_of_spectral(hAU, hAUb, "u", beginA, endA, bgScale, &daysAU);
	sum_of_spectral(hAM, hAMb, "m", beginA, endA, bgScale, &daysAM);
	sum_of_spectral(hAD, hADb, "d", beginA, endA, bgScale, &daysAD);
	sum_of_spectral(hBU, hBUb, "u", beginB, endB, bgScale, &daysBU);
	sum_of_spectral(hBM, hBMb, "m", beginB, endB, bgScale, &daysBM);
	sum_of_spectral(hBD, hBDb, "d", beginB, endB, bgScale, &daysBD);
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

void draw_signal_spectra(const char *from, const char *to, double bgScale = 5.6/2.5)
{
	char strs[128];
	char strl[1024];
	TLatex txt;

	sprintf(strs, "hOff_%s_%s", from, to);
	sprintf(strl, "Positron raw spectrum %s-%s;Positron energy, MeV;Events / (day * 0.25 MeV)", from, to);
	TH1D *hN = new TH1D(strs, strl, 60, 1, 16);
	sprintf(strs, "hOffCosm_%s_%s", from, to);
	sprintf(strl, "Positron raw spectrum (cosmic) %s-%s;Positron energy, MeV;Events / (day * 0.25 MeV)", from, to);
	TH1D *hC = new TH1D(strs, strl, 60, 1, 16);
	sum_of_raw(hN, hC, "umdrs", 31, from, to);
	hN->Add((TH1D*)gROOT->FindObject("hUp_28"), -0.006);		// subtract blocks 1-3
	TF1 fBgndN("fBgndN", "0.01446-0.000621*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.07252-0.003560*x", 0, 100);
	hC->Add(&fBgndC, -86.4);
	hC->Scale(bgScale * 0.025);
	hC->Add(&fBgndN, 86.4);
	hN->Write();
	hC->Write();
	hC->SetLineColor(kRed);
	hC->GetYaxis()->SetLabelSize(0.08);
	hC->SetTitleSize(0.08);
	hC->SetMinimum(0);
	hN->SetLineColor(kBlue);
	hN->Draw();
	hC->Draw("same");
	TLegend *lg = new TLegend(0.65, 0.75, 0.90, 0.9);
	lg->AddEntry(hN, "Neutrino", "le");
	sprintf(strs, "Cosmic * %4.1f%%", bgScale * 2.5);
	lg->AddEntry(hC, strs, "le");
	lg->Draw();
	sprintf(strs, "Integral 1-7 MeV = %5.1f events/day", hN->Integral(hN->FindBin(1.001), hN->FindBin(6.999)));
	txt.SetTextColor(kBlue);
	txt.DrawLatex(6, 15, strs);
	sprintf(strs, "Integral 1-7 MeV = %5.1f events/day", hC->Integral(hN->FindBin(1.001), hN->FindBin(6.999)));
	txt.SetTextColor(kRed);
	txt.DrawLatex(6, 10, strs);
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
	
	printf("Spectra pages\n");
//	Page 1:	Spectra All
	draw_spectra_page(cv, "Apr 16-Mar 18", 0x1E, bgScale);
	cv->Print(pname);
	
//	Page 1a:	Spectra All but April-June 16.
	draw_spectra_page(cv, "Oct 16-Mar 18", 0x1C, bgScale);
	cv->Print(pname);

//	Page 1b:	before Sep17
	draw_spectra_pagel(cv, "Oct 16-Sep 17", "01.10.16", "30.09.17", bgScale);
	cv->Print(pname);
	
//	3 months Before reactor shut down
	draw_spectra_pagel(cv, "Apr-July 17", "05.04.17", "07.07.17", bgScale);
	cv->Print(pname);

//	3 months After reactor on + 1 month
	draw_spectra_pagel(cv, "Oct-Dec 17", "20.09.17", "20.12.17", bgScale);
	cv->Print(pname);

//	Page 2:	Spectra April-June
	draw_spectra_page(cv, "April-June 16", 2, bgScale);
	cv->Print(pname);

//	Page 3:	Spectra October-November
	draw_spectra_page(cv, "Oct 16 - Feb 17", 4, bgScale);
	cv->Print(pname);

//	Page 4:	Spectra December-January
	draw_spectra_page(cv, "March-July 17", 8, bgScale);
	cv->Print(pname);

//	Page 5:	Spectra February-March
	draw_spectra_page(cv, "Aug. 17-Mar. 18", 0x10, bgScale);
	cv->Print(pname);

//	Some more spectra pages
	draw_spectra_page(cv, "Oct 16 - Jul 17", 12, bgScale);
	cv->Print(pname);

	draw_spectra_page(cv, "Mar 17 - Mar 18", 24, bgScale);
	cv->Print(pname);

	draw_spectra_pagel(cv, "Aug 17 - Jan 18", "01.08.17", "31.01.18", bgScale);
	cv->Print(pname);

	printf("Sum spectra page\n");
//	Page 6: Total Spectrum (all positions)
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	TH1D *hSum = new TH1D("hSum", "Positron spectrum Oct 16 - Mar 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	TH1D *hSumSub = new TH1D("hSumSub", "Background spectrum Oct 16 - Mar 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum, hSumSub, "umdrs", 0x1E, bgScale);
	hSum->Write();
	hSumSub->Write();
	hSum->Draw("e");

	cv->cd(2);
	TH1D *hSum1 = new TH1D("hSumAfter", "Positron spectrum Oct - Dec 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	TH1D *hSum1Sub = new TH1D("hSumAfterSub", "Background spectrum Oct - Dec 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectral(hSum1, hSum1Sub, "umdrs", "20.09.17", "20.12.17", bgScale);
	hSum1->Write();
	hSum1Sub->Write();
	hSum1->Draw("e");

	cv->cd(3);
	TH1D *hSum2 = new TH1D("hSum2", "Positron spectrum Oct 16-Feb 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	TH1D *hSum2Sub = new TH1D("hSum2Sub", "Background spectrum Oct 16-Feb 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum2, hSum2Sub, "umdrs", 4, bgScale);
	hSum2->Write();
	hSum2Sub->Write();
	hSum2->Draw("e");

	cv->cd(4);
	TH1D *hSum3 = new TH1D("hSum3", "Positron spectrum March-July 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	TH1D *hSum3Sub = new TH1D("hSum3Sub", "Background spectrum March-July 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum3, hSum3Sub, "umdrs", 8, bgScale);
	hSum3->Write();
	hSum3Sub->Write();
	hSum3->Draw("e");

	cv->cd(5);
	TH1D *hSum4 = new TH1D("hSum4", "Positron spectrum Aug 17 - Mar 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	TH1D *hSum4Sub = new TH1D("hSum4Sub", "Background spectrum Aug 17 - Mar 18;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum4, hSum4Sub, "umdrs", 0x10, bgScale);
	hSum4->Write();
	hSum4Sub->Write();
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
	lg->AddEntry(hSum4, "Aug 17 - Mar 18", "LE");
	hSum23->Draw("e");
	hSum4->Draw("e,same");
	lg->Draw();

	cv->Update();
	cv->Print(pname);
	
	printf("Ratio pages\n");
//	Page 7:	Down/Up
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hDown_4", "hUp_4", "hDownUp_4", "Ratio Down/Up Oct 16-Feb 17;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(2);
	draw_single_ratio("hDown_12", "hUp_12", "hDownUp_12", "Ratio Down/Up Oct 16-Jul 17;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(3);
	draw_single_ratio("hDown_01.10.16_30.09.17", "hUp_01.10.16_30.09.17", "hDownUp_01.10.16_30.09.17", "Ratio Down/Up Oct 16-Sep 17;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(4);
	draw_single_ratio("hDown_28", "hUp_28", "hDownUp_28", "Ratio Down/Up Oct 16-Mar 18;Positron energy, MeV", 0.6, 0.9, 24);

	cv->Update();
	cv->Print(pname);
	
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hDown_8", "hUp_8", "hDownUp_8", "Ratio Down/Up Feb-July 17;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(2);
	draw_single_ratio("hDown_16", "hUp_16", "hDownUp_16", "Ratio Down/Up Aug 17-Mar 18;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(3);
	draw_single_ratio("hDown_24", "hUp_24", "hDownUp_24", "Ratio Down/Up Mar 17-Mar 18;Positron energy, MeV", 0.6, 0.9, 24);

	cv->cd(4);
	draw_single_ratio("hDown_01.08.17_31.01.18", "hUp_01.08.17_31.01.18", "hDownUp_01.08.17_31.01.18", "Ratio Down/Up Aug 17-Jan 18;Positron energy, MeV", 0.6, 0.9, 24);

	cv->Update();
	cv->Print(pname);

//	Page 7a: Main plot
	cv->Clear();
	draw_single_ratio("hDown_28", "hUp_28", "hDownUpAll", "Ratio Down/Up Oct 16-Mar 18 (ALL);Positron energy, MeV;#frac{N_{DOWN}}{N_{UP}}", 0.65, 0.8, 24);
	cv->Update();
	cv->Print(pname);

//	Page 8:	Mid/Up
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hMid_30", "hUp_30", "hMidUpAll", "Ratio Middle/Up All;Positron energy, MeV", 0.7, 1, 24);
	draw_single_ratio("hMid_28", "hUp_28", "hMidUp0", "Ratio Middle/Up All but Apr-June 16;Positron energy, MeV", 0.7, 1, 24);

	cv->cd(2);
	draw_single_ratio("hMid_01.10.16_30.09.17", "hUp_01.10.16_30.09.17", "hMidUp1", "Ratio Middle/Up Oct 16-Sep 17;Positron energy, MeV", 0.7, 1, 24);

	cv->cd(3);
	draw_single_ratio("hMid_4", "hUp_4", "hMidUp2", "Ratio Middle/Up Oct 16 - Feb 17;Positron energy, MeV", 0.7, 1, 24);

	cv->cd(4);
	draw_single_ratio("hMid_8", "hUp_8", "hMidUp3", "Ratio Middle/Up March-July 17;Positron energy, MeV", 0.7, 1, 24);

	cv->cd(5);
	draw_single_ratio("hMid_16", "hUp_16", "hMidUp4", "Ratio Middle/Up Aug 17-Mar 18;Positron energy, MeV", 0.7, 1, 24);

	cv->Update();
	cv->Print(pname);

	printf("Dauble ratios\n");
//	Page 8.1: Double ratios Down/Up between periods.
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hDownUp_4", "hDownUp_28", "hDownUpdr1a", "Double ratio Down/Up Oct 16-Feb 17 / All;Positron energy, MeV", 0.8, 1.2, 24);

	cv->cd(2);
	draw_single_ratio("hDownUp_12", "hDownUp_28", "hDownUpdr2a", "Double ratio Down/Up Oct 16-Jul 17 / All;Positron energy, MeV", 0.8, 1.2, 24);

	cv->cd(3);
	draw_single_ratio("hDownUp_01.10.16_30.09.27", "hDownUp_28", "hDownUpdr3a", "Double ratio Down/Up Oct 16-Sep 17 / All;Positron energy, MeV", 0.8, 1.2, 24);

	cv->cd(4);
	draw_single_ratio("hDownUp_4", "hDownUp_24", "hDownUpdr1r", "Double ratio Down/Up Oct 16-Feb 17 / Mar 17-Mar 18;Positron energy, MeV", 0.8, 1.2, 24);

	cv->cd(5);
	draw_single_ratio("hDownUp_12", "hDownUp_16", "hDownUpdr2r", "Double ratio Down/Up Oct 16-Jul 17 / Aug 17-Mar 18;Positron energy, MeV", 0.8, 1.2, 24);

	cv->cd(6);
	draw_single_ratio("hDownUp_01.10.16_30.09.27", "hDownUp_01.08.17_31.01.18", "hDownUpdr3r", 
		"Double ratio Down/Up Oct 16-Sep 17 / Aug 17-Jan 18;Positron energy, MeV", 0.8, 1.2, 24);

	cv->Update();
	cv->Print(pname);
//	Page 9: Period ratios
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	draw_single_ratio("hSum4", "hSum23", "hRatio4_23", "Unnormalized Ratio Aug 17-Mar 18/Oct 16-Feb 17;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hSum3", "hSum2", "hRatio3_2", "Unnormalized Ratio March-July 17/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hSum4", "hSum2", "hRatio4_2", "Unnormalized Ratio Aug 17-Mar 18/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(4);
//	draw_single_ratio("hSum2", "hSum1", "hRatio2_1", "Unnormalized Ratio Oct 16 - Feb 17/April-June 16;Positron energy, MeV");

//	cv->cd(5);
//	draw_single_ratio("hSum3", "hSum1", "hRatio3_1", "Unnormalized Ratio March-July 17/April-June 16;Positron energy, MeV");

//	cv->cd(6);
//	draw_single_ratio("hSum4", "hSum1", "hRatio4_1", "Unnormalized Ratio Aug-Sep 17/April-June 16;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 9a: after shutdown / before shutdown - use proper periods.
	cv->Clear();
	draw_normalized_ratio("20.09.17", "20.12.17", "05.04.17", "07.07.17", "hNRatioAfter2Before", 
		"Normalized ratio after shutdown / before shutdown;Positron energy, MeV", 0.9, 1.4, 24, bgScale);
	cv->Update();
	cv->Print(pname);

//	Page 9b: after shutdown / before shutdown
	cv->Clear();
	draw_period_ratios(4, 3, "Ratio after shutdown / before shutdown;Positron energy, MeV");
	cv->Update();
	cv->Print(pname);

	printf("Raw pages\n");
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

//		Reactor Off
	cv->Clear();
	draw_signal_spectra("10.07.17", "18.08.17", bgScale);
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
	TH1D *hUpSub = new TH1D("hUpSub", "Background spectrum All, UP;Positron energy, MeV;Events/(day*0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hUp, hUpSub, "u", 30, bgScale);
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
	TH1D *hDownSub = new TH1D("hDownSub", "Background spectrum All, DOWN;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hDown, hDownSub, "d", 30, bgScale);
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

void CompareSpectra(const char *nameA, const char *nameB)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	TFile *fA = new TFile(nameA);
	TFile *fB = new TFile(nameB);
	
	TH1 *hUpA   = (TH1 *) fA->Get("hUp_30");
	TH1 *hUpB   = (TH1 *) fB->Get("hUp_30");
	TH1 *hMidA  = (TH1 *) fA->Get("hMid_30");
	TH1 *hMidB  = (TH1 *) fB->Get("hMid_30");
	TH1 *hDownA = (TH1 *) fA->Get("hDown_30");
	TH1 *hDownB = (TH1 *) fB->Get("hDown_30");
	
	if (!(hUpA && hUpB && hMidA && hMidB && hDownA && hDownB)) return;
	
	TH1 *hDiffUp    = (TH1 *)hUpA->Clone("hDiffUp");
	hDiffUp->SetTitle("Position Up, normalized difference;MeV;Difference");
	TH1 *hRatioUp   = (TH1 *)hUpA->Clone("hRatioUp");
	hRatioUp->SetTitle("Position Up, ratio;MeV;Ratio");
	TH1 *hDiffMid   = (TH1 *)hUpA->Clone("hDiffMid");
	hDiffMid->SetTitle("Position Mid, normalized difference;MeV;Difference");
	TH1 *hRatioMid  = (TH1 *)hUpA->Clone("hRatioMid");
	hRatioMid->SetTitle("Position Mid, ratio;MeV;Ratio");
	TH1 *hDiffDown  = (TH1 *)hUpA->Clone("hDiffDown");
	hDiffDown->SetTitle("Position Down, normalized difference;MeV;Difference");
	TH1 *hRatioDown = (TH1 *)hUpA->Clone("hRatioDown");
	hRatioDown->SetTitle("Position Down, ratio;MeV;Ratio");

	hRatioUp->Divide(hUpA, hUpB);
	hRatioMid->Divide(hMidA, hMidB);
	hRatioDown->Divide(hDownA, hDownB);

	hDiffUp->Add(hUpA, hUpB, 1, -hUpA->Integral()/hUpB->Integral());
	hDiffMid->Add(hMidA, hMidB, 1, -hMidA->Integral()/hMidB->Integral());
	hDiffDown->Add(hDownA, hDownB, 1, -hDownA->Integral()/hDownB->Integral());

	hRatioUp->GetXaxis()->SetRange(1, 24);
	hRatioMid->GetXaxis()->SetRange(1, 24);
	hRatioDown->GetXaxis()->SetRange(1, 24);
	hDiffUp->GetXaxis()->SetRange(1, 24);
	hDiffMid->GetXaxis()->SetRange(1, 24);
	hDiffDown->GetXaxis()->SetRange(1, 24);
	
	TCanvas *cv = new TCanvas("CV", "CV", 1200, 800);
	cv->Divide(3, 2);
	cv->cd(1);
	hRatioUp->Fit("pol0");
	cv->cd(2);
	hRatioMid->Fit("pol0");
	cv->cd(3);
	hRatioDown->Fit("pol0");
	cv->cd(4);
	hDiffUp->Fit("pol0");
	cv->cd(5);
	hDiffMid->Fit("pol0");
	cv->cd(6);
	hDiffDown->Fit("pol0");
	
}
