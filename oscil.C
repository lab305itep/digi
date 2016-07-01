#define BINMIN	6
#define BINMAX  35

void oscil(char *spfile, char *mcfile)
{
	double r;
	int i;

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);

//	Read files
	TFile *fsp = new TFile(spfile);
	if (!fsp->IsOpen()) {
		printf("File %s no opened.\n", spfile);
		return;
	}
	TFile *fmc = new TFile(mcfile);
	if (!fmc->IsOpen()) {
		printf("File %s no opened.\n", mcfile);
		return;
	}
	TH1D *hmcs = (TH1D*) fmc->Get("specEres");
	TH1D *hmcdu = (TH1D*) fmc->Get("ratioFinal_DOWNtoUP");
	TH1D *hmcmu = (TH1D*) fmc->Get("ratioFinal_MIDtoUP");
	TH1D *hmcdm = (TH1D*) fmc->Get("ratioFinal_DOWNtoMID");
	if (!hmcs || !hmcdu || !hmcmu || !hmcdm) {
		printf("Bad mc file %s\n", mcfile);
		return;
	}
	TH1D *hUp = (TH1D*) fsp->Get("hsupa");
	TH1D *hMid = (TH1D*) fsp->Get("hsstucka");
	TH1D *hDown = (TH1D*) fsp->Get("hsdowna");
	TH1D *hrud = (TH1D*) fsp->Get("hrud");
	TH1D *hrum = (TH1D*) fsp->Get("hrus");
	TH1D *hrdm = (TH1D*) fsp->Get("hrds");
	if (!hUp || !hDown || !hMid || !hrud || !hrum || !hrdm) {
		printf("Bad spectra file %s\n", spfile);
		return;		
	}

	TH1D *hone = new TH1D("HONE", "One", 40, 0, 8);
	for (i=1; i<=40; i++) {
		hone->SetBinContent(i, 1);
		hone->SetBinError(i, 0);
	}

	TH1D *hrdu = new TH1D("hrdu", "ration DOWN to UP", 40, 0, 8);
	TH1D *hrmu = new TH1D("hrmu", "ration MID to UP", 40, 0, 8);
	hrdm->GetListOfFunctions()->Clear();

	hrdu->Divide(hone, hrud);
	hrmu->Divide(hone, hrum);

	TCanvas *cvs = new TCanvas("CVS", "Spectra", 1200, 800);
	cvs->Divide(3,2);

//	Draw spectra
	hmcs->SetLineColor(kBlack);
	
	cvs->cd(1);
	r = hUp->Integral(BINMIN, BINMAX) / hmcs->Integral(BINMIN, BINMAX);
	hmcs->Scale(r);
	hUp->SetLineWidth(4);
	hUp->DrawCopy();
	hmcs->DrawCopy("same,hist");

	cvs->cd(2);
	r = hDown->Integral(BINMIN, BINMAX) / hmcs->Integral(BINMIN, BINMAX);
	hmcs->Scale(r);
	hDown->SetLineWidth(4);
	hDown->DrawCopy();
	hmcs->DrawCopy("same,hist");

	cvs->cd(3);
	r = hMid->Integral(BINMIN, BINMAX) / hmcs->Integral(BINMIN, BINMAX);
	hmcs->Scale(r);
	hMid->SetLineWidth(4);
	hMid->DrawCopy();
	hmcs->DrawCopy("same,hist");

//	Draw ratios
	
	cvs->cd(4);
	r = hrdu->Integral(BINMIN, BINMAX) / hmcdu->Integral(BINMIN, BINMAX);
	hmcdu->Scale(r);
	hrdu->SetStats(0);
	hrdu->SetMinimum(0.5);
	hrdu->SetMaximum(1.0);
	hrdu->SetLineWidth(4);
	hrdu->SetLineColor(kGreen);
	hrdu->DrawCopy();
	hmcdu->DrawCopy("hist,same");
	
	cvs->cd(5);
	r = hrmu->Integral(BINMIN, BINMAX) / hmcmu->Integral(BINMIN, BINMAX);
	hmcmu->Scale(r);
	hrmu->SetStats(0);
	hrmu->SetMinimum(0.5);
	hrmu->SetMaximum(1.0);
	hrmu->SetLineWidth(4);
	hrmu->SetLineColor(kGreen);
	hrmu->DrawCopy();
	hmcmu->DrawCopy("hist,same");

	cvs->cd(6);
	r = hrdm->Integral(BINMIN, BINMAX) / hmcdm->Integral(BINMIN, BINMAX);
	hmcdm->Scale(r);
	hrdm->SetStats(0);
	hrdm->SetMinimum(0.5);
	hrdm->SetMaximum(1.0);
	hrdm->SetLineWidth(4);
	hrdm->SetLineColor(kGreen);
	hrdm->DrawCopy();
	hmcdm->DrawCopy("hist,same");
}

