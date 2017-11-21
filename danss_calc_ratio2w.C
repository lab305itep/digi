TFile *fData;

void change_file_suffix(char *to, int len, const char *from, const char *where, const char *what)
{
	char *ptr;
	
	strncpy(to, from, len - strlen(what));
	ptr = strstr(to, where);
	if (ptr) *ptr = '\0';
	strcat(to, what);
}

int sum_of_spectra(TH1D *hSum, const char *posmask, int permask, double bgScale = 1.0)
{
#include "positions.h"
	int N;
	TH1D *hSig;
	TH1D *hSigF;
	TH1D *hBgnd;
	int i;
	char str[1024];
	double tSum;
	double dt;
	char *ptr;
	int Cnt;
//		Background tail correction
	TF1 fBgndN("fBgndN", "0.0163-0.0007*x", 0, 100);
	TF1 fBgndC("fBgndC", "0.0722-0.0034*x", 0, 100);

	N = sizeof(positions) / sizeof(positions[0]);
	hSum->Reset();
	tSum = 0;
	for (i=0; i<N; i++) {
		ptr = strchr(posmask, positions[i].name[0]);
		if (!(ptr && ((1 << positions[i].period) & permask))) continue;
		sprintf(str, "%s_hSig-diff", positions[i].name);
		hSig = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hSig", positions[i].name);
		hSigF = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		hBgnd = (TH1D*) fData->Get(str);
		if (!(hSig && hSigF && hBgnd)) continue;
		dt = hSig->Integral() / hSigF->Integral();
		tSum += dt;
		hSum->Add(hSig);
		hSum->Add(&fBgndN, -dt);
		hSum->Add(hBgnd, -positions[i].bgnd * bgScale);
		hSum->Add(&fBgndC, dt * positions[i].bgnd * bgScale);
	}
	
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
	TH1D *hSigF;
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
		sprintf(str, "%s_hSig", positions[i].name);
		hSigF = (TH1D*) fData->Get(str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		hBgnd = (TH1D*) fData->Get(str);
		if (!(hSig && hSigF && hBgnd)) continue;
		dt = hSig->Integral() / hSigF->Integral();
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
	draw_spectra_page(cv, "April 16-Sep 17", 0x1E, bgScale);
	cv->Print(pname);
	
//	Page 1a:	Spectra All but April-June 16.
	draw_spectra_page(cv, "Oct 16-Sep 17", 0x1C, bgScale);
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
	draw_spectra_page(cv, "August-September 17", 0x10, bgScale);
	cv->Print(pname);

//	Page 6: Total Spectrum (all positions)
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	TH1D *hSum = new TH1D("hSum", "Positron spectrum April 16 - Sep 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum, "umdrs", 0x1E, bgScale);
	hSum->Write();
	hSum->Draw("e");

	cv->cd(2);
	TH1D *hSum1 = new TH1D("hSum1", "Positron spectrum April-June 16;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum1, "umdrs", 2, bgScale);
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
	TH1D *hSum4 = new TH1D("hSum4", "Positron spectrum Aug-Sep 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum4, "umdrs", 0x10, bgScale);
	hSum4->Write();
	hSum4->Draw("e");

	cv->cd(6);
	TH1D *hSum12 = new TH1D("hSum12", "Positron spectrum April 16- Feb 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	hSum12->Add(hSum1, hSum2);
	hSum12->Write();
	TH1D *hSum34 = new TH1D("hSum34", "Positron spectrum March-Sep 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	hSum34->Add(hSum3, hSum4);
	hSum34->Write();
	hSum12->SetLineColor(kRed);
	hSum34->SetLineColor(kBlue);
	lg = new TLegend(0.6, 0.7, 0.9, 0.8);
	lg->AddEntry(hSum12, "April 16 - Feb 17", "LE");
	lg->AddEntry(hSum34, "March-Sep 17", "LE");
	hSum12->Draw("e");
	hSum34->Draw("e,same");
	lg->Draw();

	cv->Update();
	cv->Print(pname);

//	Page 7:	Down/Up
	cv->Clear();
	cv->Divide(3, 2);
	
	cv->cd(1);
	draw_single_ratio("hDown_30", "hUp_30", "hDownUpAll", "Ratio Down/Up All;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hDown_28", "hUp_28", "hDownUp0", "Ratio Down/Up All but Apr-June 16;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hDown_2", "hUp_2", "hDownUp1", "Ratio Down/Up April-June 16;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hDown_4", "hUp_4", "hDownUp2", "Ratio Down/Up Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hDown_8", "hUp_8", "hDownUp3", "Ratio Down/Up Feb-July 17;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hDown_16", "hUp_16", "hDownUp4", "Ratio Down/Up Aug-Sep 17;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 7a: Main plot
	cv->Clear();
	draw_single_ratio("hDown_28", "hUp_28", "hDownUp0", "Ratio Down/Up All but Apr-June 16;Positron energy, MeV;#frac{N_{DOWN}}{N_{UP}}", 0.6, 0.9, 28);
	cv->Update();
	cv->Print(pname);

//	Page 8:	Mid/Up
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hMid_30", "hUp_30", "hMidUpAll", "Ratio Middle/Up All;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hMid_28", "hUp_28", "hMidUp0", "Ratio Middle/Up All but Apr-June 16;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hMid_2", "hUp_2", "hMidUp1", "Ratio Middle/Up April-June 16;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hMid_4", "hUp_4", "hMidUp2", "Ratio Middle/Up Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hMid_8", "hUp_8", "hMidUp3", "Ratio Middle/Up March-July 17;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hMid_16", "hUp_16", "hMidUp4", "Ratio Middle/Up Aug-Sep 17;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 8.1: Double ratios Down/Up and Mid/Up between periods. Use March-July 17 as a base
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hDownUp1", "hDownUp3", "hDownUpdr1", "Ratio Down/Up Apr-June 16 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(2);
	draw_single_ratio("hDownUp2", "hDownUp3", "hDownUpdr2", "Ratio Down/Up Oct 16 - Feb 17 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(3);
	draw_single_ratio("hDownUp4", "hDownUp3", "hDownUpdr4", "Ratio Down/Up Aug-Sep 17 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(4);
	draw_single_ratio("hMidUp1", "hMidUp3", "hMidUpdr1", "Ratio Middle/Up Apr-June 16 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(5);
	draw_single_ratio("hMidUp2", "hMidUp3", "hMidUpdr2", "Ratio Middle/Up Oct 16 - Feb 17 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->cd(6);
	draw_single_ratio("hMidUp4", "hMidUp3", "hMidUpdr4", "Ratio Middle/Up Aug-Sep 17 to March-July 17;Positron energy, MeV", 0.8, 1.2, 28);

	cv->Update();
	cv->Print(pname);
//	Page 9: Period ratios
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hSum34", "hSum12", "hRatio34_12", "Ratio Mar-Sep 17/Apr 16-Feb 17;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hSum3", "hSum2", "hRatio3_2", "Ratio March-July 17/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hSum4", "hSum2", "hRatio4_2", "Ratio Aug-Sep 17/Oct 16 - Feb 17;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hSum2", "hSum1", "hRatio2_1", "Ratio Oct 16 - Feb 17/April-June 16;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hSum3", "hSum1", "hRatio3_1", "Ratio March-July 17/April-June 16;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hSum4", "hSum1", "hRatio4_1", "Ratio Aug-Sep 17/April-June 16;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 9a: after shutdown / before shutdown
	cv->Clear();
	draw_single_ratio("hSum4", "hSum3", "hRatio4_3", "Ratio after shutdown / before shutdown;Positron energy, MeV", 0.9, 1.4, 28);
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
