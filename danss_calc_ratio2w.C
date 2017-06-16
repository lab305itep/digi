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

	cv->Clear();
	lg = new TLegend(0.35, 0.65, 0.9, 0.9);
	sprintf(strs, "hUp_%d", periodmask);
	sprintf(strl, "Positron spectrum %s, UP;Positron energy, MeV;Events / (day * 0.25 MeV)", title);
	TH1D *hUp = new TH1D(strs, strl, 60, 1, 16);
	n = sum_of_spectra(hUp, "u", periodmask, bgScale);
	Cnt = n;
	TH1D *hTmp = (TH1D *) hUp->Clone("hTmp");	// keep to subtract block #3
	hUp->Add(hTmp, -0.0045);
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
	hMid->Add(hTmp, -0.0045);
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
	hDown->Add(hTmp, -0.0045);
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


void draw_single_ratio(const char *nameA, const char *nameB, const char *name, const char *title)
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
	hAB->SetMinimum(0.6);
	hAB->SetMaximum(1.2);
	hAB->GetYaxis()->SetLabelSize(0.08);
	hAB->Fit("pol0", "", "", 1, 8);
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
	draw_spectra_page(cv, "April 16-May 17", 0x1E, bgScale);
	cv->Print(pname);
	
//	Page 1a:	Spectra All but April-June 16.
	draw_spectra_page(cv, "Oct 16-May 17", 0x1C, bgScale);
	cv->Print(pname);

//	Page 2:	Spectra April-June
	draw_spectra_page(cv, "April-June 16", 2, bgScale);
	cv->Print(pname);

//	Page 3:	Spectra October-November
	draw_spectra_page(cv, "October-December 16", 4, bgScale);
	cv->Print(pname);

//	Page 4:	Spectra December-January
	draw_spectra_page(cv, "January-March 17", 8, bgScale);
	cv->Print(pname);

//	Page 5:	Spectra February-March
	draw_spectra_page(cv, "April-May17", 0x10, bgScale);
	cv->Print(pname);

//	Page 6: Total Spectrum (all positions)
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	TH1D *hSum = new TH1D("hSum", "Positron spectrum April 16 - May 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum, "umdrs", 0x1E, bgScale);
	hSum->Write();
	hSum->Draw("e");

	cv->cd(2);
	TH1D *hSum1 = new TH1D("hSum1", "Positron spectrum April-June 16;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum1, "umdrs", 2, bgScale);
	hSum1->Write();
	hSum1->Draw("e");

	cv->cd(3);
	TH1D *hSum2 = new TH1D("hSum2", "Positron spectrum October-December 16;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum2, "umdrs", 4, bgScale);
	hSum2->Write();
	hSum2->Draw("e");

	cv->cd(4);
	TH1D *hSum3 = new TH1D("hSum3", "Positron spectrum January-March 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum3, "umdrs", 8, bgScale);
	hSum3->Write();
	hSum3->Draw("e");

	cv->cd(5);
	TH1D *hSum4 = new TH1D("hSum4", "Positron spectrum April-May 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	sum_of_spectra(hSum4, "umdrs", 0x10, bgScale);
	hSum4->Write();
	hSum4->Draw("e");

	cv->cd(6);
	TH1D *hSum12 = new TH1D("hSum12", "Positron spectrum April-December 16;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	hSum12->Add(hSum1, hSum2);
	hSum12->Write();
	TH1D *hSum34 = new TH1D("hSum34", "Positron spectrum January-May 17;Positron energy, MeV;Events / (day * 0.25 MeV)", 60, 1, 16);
	hSum34->Add(hSum3, hSum4);
	hSum34->Write();
	hSum12->SetLineColor(kRed);
	hSum34->SetLineColor(kBlue);
	lg = new TLegend(0.6, 0.7, 0.9, 0.8);
	lg->AddEntry(hSum12, "April-December 16", "LE");
	lg->AddEntry(hSum34, "January-May 17", "LE");
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
	draw_single_ratio("hDown_28", "hUp_28", "hDownUp0", "Ratio Down/Up October 16 - May 17;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hDown_2", "hUp_2", "hDownUp1", "Ratio Down/Up April-June 16;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hDown_4", "hUp_4", "hDownUp2", "Ratio Down/Up October-December 16;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hDown_8", "hUp_8", "hDownUp3", "Ratio Down/Up January-March 17;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hDown_16", "hUp_16", "hDownUp4", "Ratio Down/Up April-May 17;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 8:	Mid/Up
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hMid_30", "hUp_30", "hMidUpAll", "Ratio Middle/Up All;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hMid_28", "hUp_28", "hMidUp0", "Ratio Middle/Up October 16 - May 17;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hMid_2", "hUp_2", "hMidUp1", "Ratio Middle/Up April-June 16;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hMid_4", "hUp_4", "hMidUp2", "Ratio Middle/Up October-December 17;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hMid_8", "hUp_8", "hMidUp3", "Ratio Middle/Up January-March 17;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hMid_16", "hUp_16", "hMidUp4", "Ratio Middle/Up April-May 17;Positron energy, MeV");

	cv->Update();
	cv->Print(pname);

//	Page 9: Period ratios
	cv->Clear();
	cv->Divide(3, 2);

	cv->cd(1);
	draw_single_ratio("hSum34", "hSum12", "hRatio34_12", "Ratio 17/16;Positron energy, MeV");

	cv->cd(2);
	draw_single_ratio("hSum3", "hSum2", "hRatio3_2", "Ratio January-March 17/October-December 16;Positron energy, MeV");

	cv->cd(3);
	draw_single_ratio("hSum4", "hSum2", "hRatio4_2", "Ratio April-May 17/October-December 16;Positron energy, MeV");

	cv->cd(4);
	draw_single_ratio("hSum2", "hSum1", "hRatio2_1", "Ratio October-December 16/April-June 16;Positron energy, MeV");

	cv->cd(5);
	draw_single_ratio("hSum3", "hSum1", "hRatio3_1", "Ratio January-March 17/April-June 16;Positron energy, MeV");

	cv->cd(6);
	draw_single_ratio("hSum4", "hSum1", "hRatio4_1", "Ratio April-May 17/April-June 16;Positron energy, MeV");

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
