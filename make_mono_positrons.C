#define NUMBER	12
#define NH	6
void make_mono_positrons(void)
{
	TCanvas *cv;
	const char mfiles[NUMBER][128] = {
		"/space/danss_root3/withdead/mc_positron0-5MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron1MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron1-5MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron2MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron3MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron4MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron5MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron6MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron7MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron8MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron9MeV_simple.root",
		"/space/danss_root3/withdead/mc_positron10MeV_simple.root"
	};
	const double pe[NUMBER] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	const char nhits[NH][64] = {"PositronHits == 1", "PositronHits == 2", "PositronHits == 3", "PositronHits == 4", "PositronHits == 5", "PositronHits > 5"};
	int i, j;
	char str[1024];
	char strs[128];
	TH1D *h[NUMBER][NH];
	TFile *f[NUMBER];
	TTree *t[NUMBER];
	TGraphErrors *gr[NH];
	TH1 *hc;
	TPaveStats *pv;
	TVirtualPad *vc;
	double mn, mx;
	double E[NH][NUMBER];
	double EE[NH][NUMBER];
	const double epe[NUMBER] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
//	TCut cValid("(PositronFlags & 4) == 0");
	TFile *fOut = new TFile("mc_mono_positrons.root", "RECREATE");
	for (i=0; i<NUMBER; i++) for (j=0; j<NH; j++) {
		sprintf(str, "Reconstructed positron energy from %5.1f MeV, %s;MeV;", pe[i], nhits[j]);
		sprintf(strs, "PositronEnergy_%5.1f_%d", pe[i], j);
		h[i][j] = new TH1D(strs, str, 64, 0, 16);
	}
	for (i=0; i<NUMBER; i++) {
		f[i] = new TFile(mfiles[i]);
		if (!f[i]->IsOpen()) return;
		t[i] = (TTree *) f[i]->Get("DanssEvent");
		if (!t[i]) {
			printf("Something is wrong with file %s\n", mfiles[i]);
			return;
		}
	}
	fOut->cd();
	for (i=0; i<NUMBER; i++) for (j=0; j<NH; j++)  t[i]->Project(h[i][j]->GetName(), "PositronEnergy", cX && cY && cZ && nhits[j]);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TF1 *fGaus = new TF1("fGaus", "gaus", 0, 16);
	
	for (j=0; j<NH; j++) {
		sprintf(strs, "CV%d", j);
		cv  = new TCanvas(strs, nhits[j], 1200, 900);
		cv->Divide(4, 3);
		
		for (i=0; i<NUMBER; i++) {
			cv->cd(i + 1);
			fGaus->SetParameters(10000, pe[i], 0.6);
			mn = pe[i] - sqrt(pe[i]);
			if (mn < 0) mn = 0;
			mx = pe[i] + 1.5*sqrt(pe[i]);
			if (mx > 16) mx  = 16;
			h[i][j]->Fit(fGaus, "Q", "0", mn, mx);
			E[j][i] = fGaus->GetParameter(1) - pe[i];
			EE[j][i] = fGaus->GetParError(1);
			h[i][j]->DrawCopy();
		}
	}

	cv = new TCanvas("CV", "Graphs", 1200, 800);
	cv->Divide(3, 2);
	TH1D *hFrame = new TH1D("hFrame", ";E_{MC},MeV;#Delta E, MeV", 22, 0, 11);
	hFrame->SetMinimum(-1.5);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	for (j=0; j<NH; j++) {
		gr[j] = new TGraphErrors(NUMBER, pe, E[j], epe, EE[j]);
		gr[j]->SetMarkerColor(kBlue);
		gr[j]->SetLineColor(kBlue);
		gr[j]->SetMarkerStyle(kFullCircle);
		vc = cv->cd(j + 1);
		hFrame->SetTitle(nhits[j]);
		hc = hFrame->DrawCopy();
		gr[j]->Draw("p");
		gr[j]->Fit("pol1");
		vc->Update();
		pv = (TPaveStats *)gr[j]->FindObject("stats");
		pv->SetTextSize(0.05);
		pv->SetX1NDC(0.35);
		pv->SetY1NDC(0.7);
		pv->Draw();
		vc->Update();
	}

	for (i=0; i<NUMBER; i++) for (j=0; j<NH; j++) h[i][j]->Write();
	for (i=0; i<NUMBER; i++) f[i]->Close();
	fOut->Close();
}

TGraphErrors *make_mono_positrons_all(const char *subdir = "")
{
	TCanvas *cv;
	const char mfiles[NUMBER][128] = {
		"/space/danss_root3/%s/mc_positron0-5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron1MeV_simple.root",
		"/space/danss_root3/%s/mc_positron1-5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron2MeV_simple.root",
		"/space/danss_root3/%s/mc_positron3MeV_simple.root",
		"/space/danss_root3/%s/mc_positron4MeV_simple.root",
		"/space/danss_root3/%s/mc_positron5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron6MeV_simple.root",
		"/space/danss_root3/%s/mc_positron7MeV_simple.root",
		"/space/danss_root3/%s/mc_positron8MeV_simple.root",
		"/space/danss_root3/%s/mc_positron9MeV_simple.root",
		"/space/danss_root3/%s/mc_positron10MeV_simple.root"
	};
	const double pe[NUMBER] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	double spe[NUMBER];
	int i;
	char str[1024];
	char strs[128];
	TH1D *h[NUMBER];
	TFile *f[NUMBER];
	TTree *t[NUMBER];
	TGraphErrors *gr;
	TH1 *hc;
	TPaveStats *pv;
	TVirtualPad *vc;
	double mn, mx;
	double E[NUMBER];
	double EE[NUMBER];
	double SE[NUMBER];
	double SEE[NUMBER];
	const double epe[NUMBER] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
//	TCut cValid("(PositronFlags & 4) == 0");
	TFile *fOut = new TFile("mc_mono_positrons.root", "UPDATE");
	for (i=0; i<NUMBER; i++) {
		sprintf(str, "Reconstructed positron energy from %5.1f MeV;MeV;", pe[i]);
		sprintf(strs, "PositronEnergy_%5.1f", pe[i]);
		h[i] = new TH1D(strs, str, 64, 0, 16);
	}
	for (i=0; i<NUMBER; i++) {
		sprintf(str, mfiles[i], subdir);
		f[i] = new TFile(str);
		if (!f[i]->IsOpen()) return NULL;
		t[i] = (TTree *) f[i]->Get("DanssEvent");
		if (!t[i]) {
			printf("Something is wrong with file %s\n", mfiles[i]);
			return NULL;
		}
	}
	fOut->cd();
	for (i=0; i<NUMBER; i++) t[i]->Project(h[i]->GetName(), "PositronEnergy", cX && cY && cZ);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetStatW(0.5);
	gStyle->SetStatH(0.25);
	TF1 *fGaus = new TF1("fGaus", "gaus", 0, 16);
	
	cv  = new TCanvas("CVA", "All clusters", 1200, 900);
	cv->Divide(4, 3);
		
	for (i=0; i<NUMBER; i++) {
		cv->cd(i + 1);
		fGaus->SetParameters(10000, pe[i], 0.6);
		mn = pe[i] - sqrt(pe[i]);
		if (mn < 0) mn = 0;
		mx = pe[i] + 1.5*sqrt(pe[i]);
		if (mx > 16) mx  = 16;
		h[i]->Fit(fGaus, "Q", "0", mn, mx);
		E[i] = fGaus->GetParameter(1) - pe[i];
		EE[i] = fGaus->GetParError(1);
		SE[i] = fGaus->GetParameter(2);
		SEE[i] = fGaus->GetParError(2);
		h[i]->DrawCopy();
	}

	cv = new TCanvas("CVG", "Graph", 1200, 800);
	TH1D *hFrame = new TH1D("hFrame", ";E_{MC},MeV;#Delta E, MeV", 22, 0, 11);
	hFrame->SetMinimum(-1.5);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	gr = new TGraphErrors(NUMBER, pe, E, epe, EE);
	gr->SetMarkerColor(kBlue);
	gr->SetLineColor(kBlue);
	gr->SetMarkerStyle(kFullCircle);
	hc = hFrame->DrawCopy();
	gr->Draw("p");
	gr->Fit("pol1");
	cv->Update();
	pv = (TPaveStats *)gr->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.35);
	pv->SetY1NDC(0.7);
	pv->Draw();
	cv->Update();

	for (i=0; i<NUMBER;i++) spe[i] = sqrt(pe[i]);
	cv = new TCanvas("CVGE", "Graph Error", 1200, 800);
	hFrame = new TH1D("hFrameE", ";#sqrt{E_{MC}}, #sqrt{MeV};#sigma, MeV", 10, 0, 3.5);
	hFrame->SetMinimum(0);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	gr = new TGraphErrors(NUMBER, spe, SE, epe, SEE);
	gr->SetMarkerColor(kBlue);
	gr->SetLineColor(kBlue);
	gr->SetMarkerStyle(kFullCircle);
	hc = hFrame->DrawCopy();
	gr->Draw("p");
	gr->Fit("pol1");
	cv->Update();
	pv = (TPaveStats *)gr->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.35);
	pv->SetY1NDC(0.7);
	pv->Draw();
	cv->Update();

	for (i=0; i<NUMBER; i++) h[i]->Write();
	for (i=0; i<NUMBER; i++) f[i]->Close();
	fOut->Close();
	return gr;
}

TGraphErrors *make_mono_positrons_sum(const char *subdir = "")
{
	TCanvas *cv;
	const char mfiles[NUMBER][128] = {
		"/space/danss_root3/%s/mc_positron0-5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron1MeV_simple.root",
		"/space/danss_root3/%s/mc_positron1-5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron2MeV_simple.root",
		"/space/danss_root3/%s/mc_positron3MeV_simple.root",
		"/space/danss_root3/%s/mc_positron4MeV_simple.root",
		"/space/danss_root3/%s/mc_positron5MeV_simple.root",
		"/space/danss_root3/%s/mc_positron6MeV_simple.root",
		"/space/danss_root3/%s/mc_positron7MeV_simple.root",
		"/space/danss_root3/%s/mc_positron8MeV_simple.root",
		"/space/danss_root3/%s/mc_positron9MeV_simple.root",
		"/space/danss_root3/%s/mc_positron10MeV_simple.root"
	};
	const double pe[NUMBER] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	double spe[NUMBER];
	int i;
	char str[1024];
	char strs[128];
	TH1D *h[NUMBER];
	TFile *f[NUMBER];
	TTree *t[NUMBER];
	TGraphErrors *gr;
	TH1 *hc;
	TPaveStats *pv;
	TVirtualPad *vc;
	double mn, mx;
	double E[NUMBER];
	double EE[NUMBER];
	double SE[NUMBER];
	double SEE[NUMBER];
	const double epe[NUMBER] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
//	TCut cX("PositronX[0] < 0 || (PositronX[0] > 12 && PositronX[0] < 84)");
//	TCut cY("PositronX[1] < 0 || (PositronX[1] > 12 && PositronX[1] < 84)");
//	TCut cZ("PositronX[2] > 13.5 && PositronX[2] < 85.5");
//	TCut cValid("(PositronFlags & 4) == 0");
	TFile *fOut = new TFile("mc_mono_positrons.root", "UPDATE");
	for (i=0; i<NUMBER; i++) {
		sprintf(str, "Reconstructed positron energy from %5.1f MeV;MeV;", pe[i]);
		sprintf(strs, "SumEnergy_%5.1f", pe[i]);
		h[i] = new TH1D(strs, str, 64, 0, 16);
	}
	for (i=0; i<NUMBER; i++) {
		sprintf(str, mfiles[i], subdir);
		f[i] = new TFile(str);
		if (!f[i]->IsOpen()) return NULL;
		t[i] = (TTree *) f[i]->Get("DanssEvent");
		if (!t[i]) {
			printf("Something is wrong with file %s\n", mfiles[i]);
			return NULL;
		}
	}
	fOut->cd();
	for (i=0; i<NUMBER; i++) t[i]->Project(h[i]->GetName(), "TotalEnergy", cX && cY && cZ);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetStatW(0.5);
	gStyle->SetStatH(0.25);
	TF1 *fGaus = new TF1("fGaus", "gaus", 0, 16);
	
	cv  = new TCanvas("CVA", "All clusters", 1200, 900);
	cv->Divide(4, 3);
		
	for (i=0; i<NUMBER; i++) {
		cv->cd(i + 1);
		fGaus->SetParameters(10000, pe[i], 0.6);
		mn = pe[i] - sqrt(pe[i]);
		if (mn < 0) mn = 0;
		mx = pe[i] + 1.5*sqrt(pe[i]);
		if (mx > 16) mx  = 16;
		h[i]->Fit(fGaus, "Q", "0", mn, mx);
		E[i] = fGaus->GetParameter(1) - pe[i];
		EE[i] = fGaus->GetParError(1);
		SE[i] = fGaus->GetParameter(2);
		SEE[i] = fGaus->GetParError(2);
		h[i]->DrawCopy();
	}

	cv = new TCanvas("CVG", "Graph", 1200, 800);
	TH1D *hFrame = new TH1D("hFrame", ";E_{MC},MeV;#Delta E, MeV", 22, 0, 11);
	hFrame->SetMinimum(-1.5);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	gr = new TGraphErrors(NUMBER, pe, E, epe, EE);
	gr->SetMarkerColor(kBlue);
	gr->SetLineColor(kBlue);
	gr->SetMarkerStyle(kFullCircle);
	hc = hFrame->DrawCopy();
	gr->Draw("p");
	gr->Fit("pol1");
	cv->Update();
	pv = (TPaveStats *)gr->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.35);
	pv->SetY1NDC(0.7);
	pv->Draw();
	cv->Update();

	for (i=0; i<NUMBER;i++) spe[i] = sqrt(pe[i]);
	cv = new TCanvas("CVGE", "Graph Error", 1200, 800);
	cv->SetBottomMargin(0.12);
	hFrame = new TH1D("hFrameE", ";#sqrt{E_{MC}}, #sqrt{MeV};#sigma, MeV", 10, 0, 3.5);
	hFrame->SetMinimum(0);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	gr = new TGraphErrors(NUMBER, spe, SE, epe, SEE);
	gr->SetMarkerColor(kBlue);
	gr->SetLineColor(kBlue);
	gr->SetMarkerStyle(kFullCircle);
	hc = hFrame->DrawCopy();
	gr->Draw("p");
	gr->Fit("pol1");
	cv->Update();
	pv = (TPaveStats *)gr->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.45);
	pv->SetY1NDC(0.73);
	pv->Draw();
	cv->Update();

	for (i=0; i<NUMBER; i++) h[i]->Write();
	for (i=0; i<NUMBER; i++) f[i]->Close();
	fOut->Close();
	return gr;
}

void draw_error_graphs(void)
{
	TPaveStats *pv;
	TGraphErrors *grA = make_mono_positrons_all("wd-corrA");
	TGraphErrors *grB = make_mono_positrons_all("wd-corrB");
	TCanvas *cv = new TCanvas("SIGMA", "Sigmas", 900, 800);
	cv->SetBottomMargin(0.12);
	hFrame = new TH1D("hFrameE", ";#sqrt{E_{MC}}, #sqrt{MeV};#sigma, MeV", 10, 0, 3.5);
	TLegend *lg = new TLegend(0.4, 0.12, 0.95, 0.26);
	lg->SetTextSize(0.05);
	hFrame->SetMinimum(0);
	hFrame->SetMaximum(1.5);
	hFrame->GetXaxis()->SetLabelSize(0.05);
	hFrame->GetYaxis()->SetLabelSize(0.05);
	hFrame->GetXaxis()->SetTitleSize(0.05);
	hFrame->GetYaxis()->SetTitleSize(0.05);
	grA->SetMarkerColor(kBlue);
	grA->SetLineColor(kBlue);
	grA->SetMarkerStyle(kFullCircle);
	grA->GetFunction("pol1")->SetLineColor(kBlue);
	grB->SetMarkerColor(kGreen);
	grB->SetLineColor(kGreen);
	grB->SetMarkerStyle(kFullSquare);
	grB->GetFunction("pol1")->SetLineColor(kGreen);
	hFrame->DrawCopy();
	grA->Draw("p");
	grB->Draw("p");
	cv->Update();
	pv = (TPaveStats *)grA->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.12);
	pv->SetX2NDC(0.6);
	pv->SetY1NDC(0.6);
	pv->SetY2NDC(0.75);
	pv->SetLineColor(kBlue);
	pv->SetTextColor(kBlue);
	pv->Draw();
	pv = (TPaveStats *)grB->FindObject("stats");
	pv->SetTextSize(0.05);
	pv->SetX1NDC(0.12);
	pv->SetX2NDC(0.6);
	pv->SetY1NDC(0.75);
	pv->SetY2NDC(0.9);
	pv->SetLineColor(kGreen);
	pv->SetTextColor(kGreen);
	lg->AddEntry(grA, "Hit number correction", "pel");
	lg->AddEntry(grB, "Average correction", "pel");
	lg->Draw();
	cv->Update();
}
