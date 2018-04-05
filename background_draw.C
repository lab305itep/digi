#define NHISTS 17

void background_draw(const char *rootname = "background_plots2.root")
{
	char strs[128];
	char strl[1024];
	const char titles[NHISTS][16] = {"gtDiff", "R1", "R2", "RZ", "PX", "PY", "PZ", "NX", "NY", "NZ", "NE", "NH", "PH", "AH", "AE", "AM", "AMO"};
	const char suffix[4][10] = {"A-rand", "A-diff", "B-diff", "C-diff"};
	const Color_t color[4] = {kGreen+2, kBlue, kRed, kOrange};
	const int marker[4] = {kOpenCircle, kFullCircle, kOpenSquare, kOpenStar};
	TH1D *h[NHISTS][4];
	int i, j;
	double hMax;
	int iMax;
	char pdfname[1024];
	char *ptr;

	strcpy(pdfname, rootname);
	ptr = strstr(pdfname, ".root");
	if (ptr) {
		strcpy(ptr, ".pdf");
	} else {
		strcat(pdfname, ".pdf");
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	
	TFile *fRoot = new TFile(rootname);
	if (!fRoot->IsOpen()) {
		printf("root-file %s not found  - run background_calc() first!\n", rootname);
		return;
	}
	for (i=0; i<NHISTS; i++) {
		for (j=0; j<4; j++) {
			sprintf(strs, "h%s%s", titles[i], suffix[j]);
			h[i][j] = (TH1D*) fRoot->Get(strs);
			if (!h[i][j]) {
				printf("%s not found  - rerun background_calc() to create all hists!\n", strs);
				fRoot->Close();
				return;
			}
			h[i][j]->SetLineWidth(2);
			h[i][j]->SetLineColor(color[j]);
			h[i][j]->SetMarkerColor(color[j]);
			h[i][j]->SetMarkerStyle(marker[j]);
			h[i][j]->GetYaxis()->SetLabelSize(0.05);
			h[i][j]->SetMinimum(0);
			h[i][j]->GetYaxis()->SetTitle("");
		}
	}
	
	TCanvas *cv = new TCanvas("CV", "Background plots", 1200, 900);
	TLegend *lg = new TLegend(0.7, 0.8, 0.95, 0.95);
	lg->AddEntry(h[0][0], "Random", "LP");
	lg->AddEntry(h[0][1], "Neutrino", "LP");
	lg->AddEntry(h[0][2], "Cosmic-A", "LP");
	lg->AddEntry(h[0][3], "Cosmic-B", "LP");
	
	sprintf(strl, "%s[", pdfname);
	cv->SaveAs(strl);
	
	for (i=0; i<NHISTS; i++) {
		cv->Clear();

		hMax = 0;
		iMax = 0;
		for (j=0; j<3; j++) if (h[i][j]->GetMaximum() > hMax) {
			hMax = h[i][j]->GetMaximum();
			iMax = j;
		}
		h[i][iMax]->Draw();
		for (j=0; j<4; j++) if (j != iMax) h[i][j]->Draw("same");
		lg->Draw();
		cv->SaveAs(pdfname);
	}
	
	sprintf(strl, "%s]", pdfname);
	cv->SaveAs(strl);
	fRoot->Close();
}
