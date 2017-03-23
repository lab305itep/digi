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
	char *ptr;
	
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
		tSum += hSig->Integral() / hSigF->Integral();
		hSum->Add(hSig);
		hSum->Add(hBgnd, -positions[i].bgnd * bgScale);
	}
	if (tSum == 0) return 0;
	hSum->Scale(86.4 / tSum);
	return 1;
}

void danss_calc_ratio2w(const char *fname, double bgScale = 1.0)
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
	gStyle->SetLineWidth(2);
	
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
	
//	Page 1:	Spectra April-January
	cv->Clear();
	lg = new TLegend(0.55, 0.65, 0.9, 0.9);
	TH1D *hUp = new TH1D("hUp", "Positron spectrum April-January, UP;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hUp, "u", 14, bgScale);
	hUp->Write();
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.08);
	hUp->SetTitle("April-January");
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, 35, err);
	sprintf(str, "Up: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hUp, str, "l");

	TH1D *hMid = new TH1D("hMid", "Positron spectrum April-January, MIDDLE;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hMid, "m", 14, bgScale);
	hMid->Write();
	hMid->SetLineColor(kGreen);
	hMid->SetFillColor(kGreen-10);
	hMid->Draw("same,hist,e");
	val = hMid->IntegralAndError(1, 35, err);
	sprintf(str, "Middle: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hMid, str, "l");

	TH1D *hDown = new TH1D("hDown", "Positron spectrum April-January, DOWN;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hDown, "d", 14, bgScale);
	hDown->Write();
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, 35, err);
	sprintf(str, "Down: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hDown, str, "l");

	hUp->Draw("axis,same");
	lg->Draw();
	cv->Update();
	cv->Print(pname);

//	Page 2:	Spectra April-June
	cv->Clear();
	lg = new TLegend(0.55, 0.65, 0.9, 0.9);
	TH1D *hUp1 = new TH1D("hUp1", "Positron spectrum April-June, UP;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hUp1, "u", 2, bgScale);
	hUp1->Write();
	hUp1->SetLineColor(kRed);
	hUp1->SetFillColor(kRed-10);
	hUp1->GetYaxis()->SetLabelSize(0.08);
	hUp1->SetTitle("April-June");
	hUp1->Draw("hist,e");
	val = hUp1->IntegralAndError(1, 35, err);
	sprintf(str, "Up: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hUp1, str, "l");

	TH1D *hMid1 = new TH1D("hMid1", "Positron spectrum April-June, MIDDLE;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hMid1, "m", 2, bgScale);
	hMid1->Write();
	hMid1->SetLineColor(kGreen);
	hMid1->SetFillColor(kGreen-10);
	hMid1->Draw("same,hist,e");
	val = hMid1->IntegralAndError(1, 35, err);
	sprintf(str, "Middle: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hMid1, str, "l");

	TH1D *hDown1 = new TH1D("hDown1", "Positron spectrum April-June, DOWN;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hDown1, "d", 2, bgScale);
	hDown1->Write();
	hDown1->SetLineColor(kBlue);
	hDown1->SetFillColor(kBlue-10);
	hDown1->Draw("same,hist,e");
	val = hDown1->IntegralAndError(1, 35, err);
	sprintf(str, "Down: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hDown1, str, "l");

	hUp1->Draw("axis,same");
	lg->Draw();
	cv->Update();
	cv->Print(pname);

//	Page 3:	Spectra October-November
	cv->Clear();
	lg = new TLegend(0.55, 0.65, 0.9, 0.9);
	TH1D *hUp2 = new TH1D("hUp2", "Positron spectrum October-November, UP;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hUp2, "u", 4, bgScale);
	hUp2->Write();
	hUp2->SetLineColor(kRed);
	hUp2->SetFillColor(kRed-10);
	hUp2->GetYaxis()->SetLabelSize(0.08);
	hUp2->SetTitle("October-November");
	hUp2->Draw("hist,e");
	val = hUp2->IntegralAndError(1, 35, err);
	sprintf(str, "Up: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hUp2, str, "l");

	TH1D *hMid2 = new TH1D("hMid2", "Positron spectrum October-November, MIDDLE;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hMid2, "m", 4, bgScale);
	hMid2->Write();
	hMid2->SetLineColor(kGreen);
	hMid2->SetFillColor(kGreen-10);
	hMid2->Draw("same,hist,e");
	val = hMid2->IntegralAndError(1, 35, err);
	sprintf(str, "Middle: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hMid2, str, "l");

	TH1D *hDown2 = new TH1D("hDown2", "Positron spectrum October-November, DOWN;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hDown2, "d", 4, bgScale);
	hDown2->Write();
	hDown2->SetLineColor(kBlue);
	hDown2->SetFillColor(kBlue-10);
	hDown2->Draw("same,hist,e");
	val = hDown2->IntegralAndError(1, 35, err);
	sprintf(str, "Down: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hDown2, str, "l");

	hUp2->Draw("axis,same");
	lg->Draw();
	cv->Update();
	cv->Print(pname);

//	Page 4:	Spectra December-January
	cv->Clear();
	lg = new TLegend(0.55, 0.65, 0.9, 0.9);
	TH1D *hUp3 = new TH1D("hUp3", "Positron spectrum December-January, UP;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hUp3, "u", 8, bgScale);
	hUp3->Write();
	hUp3->SetLineColor(kRed);
	hUp3->SetFillColor(kRed-10);
	hUp3->GetYaxis()->SetLabelSize(0.08);
	hUp3->SetTitle(" December-January");
	hUp3->Draw("hist,e");
	val = hUp3->IntegralAndError(1, 35, err);
	sprintf(str, "Up: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hUp3, str, "l");

	TH1D *hMid3 = new TH1D("hMid3", "Positron spectrum December-January, MIDDLE;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hMid3, "m", 8, bgScale);
	hMid3->Write();
	hMid3->SetLineColor(kGreen);
	hMid3->SetFillColor(kGreen-10);
	hMid3->Draw("same,hist,e");
	val = hMid3->IntegralAndError(1, 35, err);
	sprintf(str, "Middle: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hMid3, str, "l");

	TH1D *hDown3 = new TH1D("hDown3", "Positron spectrum December-January, DOWN;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hDown3, "d", 8, bgScale);
	hDown3->Write();
	hDown3->SetLineColor(kBlue);
	hDown3->SetFillColor(kBlue-10);
	hDown3->Draw("same,hist,e");
	val = hDown3->IntegralAndError(1, 35, err);
	sprintf(str, "Down: %5.0f #pm%4.0f", val, err);
	lg->AddEntry(hDown3, str, "l");

	hUp3->Draw("axis,same");
	lg->Draw();
	cv->Update();
	cv->Print(pname);

//	Page 5: Total Spectrum (all positions)
	cv->Clear();
	cv->Divide(2, 2);

	cv->cd(1);
	TH1D *hSum = new TH1D("hSum", "Positron spectrum April-January;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hSum, "umdrs", 14, bgScale);
	hSum->Write();
	hSum->Draw("e");

	cv->cd(2);
	TH1D *hSum1 = new TH1D("hSum1", "Positron spectrum April-June;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hSum1, "umdrs", 2, bgScale);
	hSum1->Write();
	hSum1->Draw("e");

	cv->cd(3);
	TH1D *hSum2 = new TH1D("hSum2", "Positron spectrum October-November;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hSum2, "umdrs", 4, bgScale);
	hSum2->Write();
	hSum2->Draw("e");

	cv->cd(4);
	TH1D *hSum3 = new TH1D("hSum3", "Positron spectrum December-January;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hSum3, "umdrs", 8, bgScale);
	hSum3->Write();
	hSum3->Draw("e");

	cv->Update();
	cv->Print(pname);

//	Page 6:	Down/Up
	cv->Clear();
	cv->Divide(2, 2);
	
	cv->cd(1);
	TH1D *hDownUp = new TH1D("hDownUp", "Ratio Down/Up April-January;Positron energy, MeV", 55, 1, 12);
	hDownUp->Divide(hDown, hUp);
	hDownUp->Write();
	hDownUp->SetMinimum(0.5);
	hDownUp->SetMaximum(1.0);
	hDownUp->GetYaxis()->SetLabelSize(0.08);
	hDownUp->Fit("pol0", "", "", 1, 8);

	cv->cd(2);
	TH1D *hDownUp1 = new TH1D("hDownUp1", "Ratio Down/Up April-June;Positron energy, MeV", 55, 1, 12);
	hDownUp1->Divide(hDown1, hUp1);
	hDownUp1->Write();
	hDownUp1->SetMinimum(0.5);
	hDownUp1->SetMaximum(1.0);
	hDownUp1->GetYaxis()->SetLabelSize(0.08);
	hDownUp1->Fit("pol0", "", "", 1, 8);

	cv->cd(3);
	TH1D *hDownUp2 = new TH1D("hDownUp2", "Ratio Down/Up October-November;Positron energy, MeV", 55, 1, 12);
	hDownUp2->Divide(hDown2, hUp2);
	hDownUp2->Write();
	hDownUp2->SetMinimum(0.5);
	hDownUp2->SetMaximum(1.0);
	hDownUp2->GetYaxis()->SetLabelSize(0.08);
	hDownUp2->Fit("pol0", "", "", 1, 8);

	cv->cd(4);
	TH1D *hDownUp3 = new TH1D("hDownUp3", "Ratio Down/Up December-January;Positron energy, MeV", 55, 1, 12);
	hDownUp3->Divide(hDown3, hUp3);
	hDownUp3->Write();
	hDownUp3->SetMinimum(0.5);
	hDownUp3->SetMaximum(1.0);
	hDownUp3->GetYaxis()->SetLabelSize(0.08);
	hDownUp3->Fit("pol0", "", "", 1, 8);

	cv->Update();
	cv->Print(pname);

//	Page 7:	1:2:3 Mid/Up
	cv->Clear();
	cv->Divide(2, 2);
	
	cv->cd(1);
	TH1D *hMidUp = new TH1D("hMidUp", "Ratio Mid/Up April-January;Positron energy, MeV", 55, 1, 12);
	hMidUp->Divide(hMid, hUp);
	hMidUp->Write();
	hMidUp->SetMinimum(0.6);
	hMidUp->SetMaximum(1.2);
	hMidUp->GetYaxis()->SetLabelSize(0.08);
	hMidUp->Fit("pol0", "", "", 1, 8);

	cv->cd(2);
	TH1D *hSum12 = new TH1D("hRat12", "Ratio April-June / October-November;Positron energy, MeV", 55, 1, 12);
	hSum12->Divide(hSum1, hSum2);
	hSum12->Write();
	hSum12->SetMinimum(0.7);
	hSum12->SetMaximum(1.3);
	hSum12->GetYaxis()->SetLabelSize(0.08);
	hSum12->Fit("pol0", "", "", 1, 8);

	cv->cd(3);
	TH1D *hSum13 = new TH1D("hRat13", "Ratio April-June / December-January;Positron energy, MeV", 55, 1, 12);
	hSum13->Divide(hSum1, hSum3);
	hSum13->Write();
	hSum13->SetMinimum(0.7);
	hSum13->SetMaximum(1.3);
	hSum13->GetYaxis()->SetLabelSize(0.08);
	hSum13->Fit("pol0", "", "", 1, 8);

	cv->cd(4);
	TH1D *hSum23 = new TH1D("hRat23", "Ratio October-November / December-January;Positron energy, MeV", 55, 1, 12);
	hSum23->Divide(hSum2, hSum3);
	hSum23->Write();
	hSum23->SetMinimum(0.7);
	hSum23->SetMaximum(1.3);
	hSum23->GetYaxis()->SetLabelSize(0.08);
	hSum23->Fit("pol0", "", "", 1, 8);

	cv->Update();
	cv->Print(pname);

	sprintf(str, "%s]", pname);
	cv->Print(str);
	fOut->Close();
	fData->Close();
	delete cv;
}

void danss_draw_sp2(const char *fname, double bgScale = 1)
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
	gStyle->SetLineWidth(2);
	
	cv = new TCanvas("CV", "Results", 1200, 900);
	fData = new TFile(fname);
	if (!fData->IsOpen()) return;
	change_file_suffix(pname, sizeof(pname), fname, ".root", "-sp2.png");
	txt = new TLatex();
	
	lg = new TLegend(0.45, 0.75, 0.9, 0.9);
	TH1D *hUp = new TH1D("hUp", "Positron spectrum April-January, UP;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hUp, "u", 14, bgScale);
	hUp->SetLineColor(kRed);
	hUp->SetFillColor(kRed-10);
	hUp->GetYaxis()->SetLabelSize(0.06);
	hUp->SetTitle("");
	hUp->Draw("hist,e");
	val = hUp->IntegralAndError(1, 35, err);
	sprintf(str, "Up:     %5.0f #pm %2.0f / day", val, err);
	lg->AddEntry(hUp, str, "l");

	TH1D *hDown = new TH1D("hDown", "Positron spectrum April-January, DOWN;Positron energy, MeV;Events per day per 0.2 MeV", 55, 1, 12);
	sum_of_spectra(hDown, "d", 14, bgScale);
	hDown->SetLineColor(kBlue);
	hDown->SetFillColor(kBlue-10);
	hDown->Draw("same,hist,e");
	val = hDown->IntegralAndError(1, 35, err);
	sprintf(str, "Down: %5.0f #pm %2.0f / day", val, err);
	lg->AddEntry(hDown, str, "l");
	
	TFile *fMc = new TFile("mc_fuel_middle.root");
	TH1D *hMc = (TH1D*) fMc->Get("hMcMiddleMixt");
	if (hMc) {
		hMc->SetMarkerStyle(kFullStar);
		hMc->SetMarkerColor(kGreen+2);
		hMc->SetMarkerSize(2);
		hMc->Scale(hUp->Integral(1, 35) / hMc->Integral(1, 35));
		hMc->Draw("same");
		lg->AddEntry(hMc, "Monte Carlo", "p");
	}

	hUp->Draw("axis,same");
	lg->Draw();
	cv->Print(pname);
}
