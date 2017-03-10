TFile *fData;

void change_file_suffix(char *to, int len, const char *from, const char *where, const char *what)
{
	char *ptr;
	
	strncpy(to, from, len - strlen(what));
	ptr = strstr(to, where);
	if (ptr) *ptr = '\0';
	strcat(to, what);
}

int sum_of_spectra(TH1D *hSum, const char *posmask, int permask)
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
		hSum->Add(hBgnd, -positions[i].bgnd);
	}
	if (tSum == 0) return 0;
	hSum->Scale(86.4 / tSum);
	return 1;
}

void danss_calc_ratio2(const char *fname)
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
	lg = new TLegend(0.6, 0.7, 0.9, 0.9);
	
	TH1D *hUp = new TH1D("hUp", "Positron spectrum, April-January;MeV;Events per day per 0.2 MeV", 35, 1, 8);
	sum_of_spectra(hUp, "u", 14);
	hUp->SetLineWidth(2);
	hUp->SetLineColor(kRed);
	hUp->Draw();
	val = hUp->IntegralAndError(1, 35, err);
	sprintf(str, "UP: %5.0f#pm%4.0f", val, err);
	lg->AddEntry(hUp, str, "l");

	TH1D *hMid = new TH1D("hMid", "Positron spectrum April-January;MeV;Events per day per 0.2 MeV", 35, 1, 8);
	sum_of_spectra(hMid, "m", 14);
	hMid->SetLineWidth(2);
	hMid->SetLineColor(kGreen);
	hMid->Draw("same");
	val = hMid->IntegralAndError(1, 35, err);
	sprintf(str, "MID: %5.0f#pm%4.0f", val, err);
	lg->AddEntry(hMid, str, "l");

	TH1D *hDown = new TH1D("hDown", "Positron spectrum April-January;MeV;Events per day per 0.2 MeV", 35, 1, 8);
	sum_of_spectra(hDown, "d", 14);
	hDown->SetLineWidth(2);
	hDown->SetLineColor(kBlue);
	hDown->Draw("same");
	val = hDown->IntegralAndError(1, 35, err);
	sprintf(str, "DOWN: %5.0f#pm%4.0f", val, err);
	lg->AddEntry(hDown, str, "l");

	lg->Draw();	
	cv->Update();
	cv->Print(pname);

	sprintf(str, "%s]", pname);
	cv->Print(str);
	fOut->Close();
	fData->Close();
	delete cv;
}

