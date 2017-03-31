void xy_plots_draw(void)
{
	const char position[3][10] = {"DOWN", "MIDDLE", "UP"};
	const Color_t color[3][2] = {{kGreen, kOrange+4}, {kCyan, kViolet}, {kBlue, kRed}};
	char str[128];
	int i, j, k;
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);

	TFile *fRoot = new TFile("xy_plots.root");
	if (!fRoot->IsOpen()) {
		printf("Do .xy_plots.C first\n");
		return;
	}

	TH1D *h[2][3][2];	// X/Y, D/M/U, S/C
	for (i=0; i<2; i++) for (j=0; j<3; j++) for (k=0; k<2; k++) {
		sprintf(str, "h%c%d%c", (i) ? 'Y' : 'X', j, (k) ? 'C' : 'S');
		h[i][j][k] = (TH1D*) fRoot->Get(str);
		if (!h[i][j][k]) {
			printf("%s not found\n", str);
			return;
		}
		h[i][j][k]->SetLineWidth(2);
		h[i][j][k]->Scale(21.6);
		sprintf(str, ";%c, cm;Events/day/cm", (i) ? 'Y' : 'X');
		h[i][j][k]->SetTitle(str);
		h[i][j][k]->GetXaxis()->SetLabelSize(0.06);
		h[i][j][k]->GetXaxis()->SetTitleSize(0.06);
		h[i][j][k]->GetYaxis()->SetLabelSize(0.06);
		h[i][j][k]->GetYaxis()->SetTitleSize(0.06);
		h[i][j][k]->SetLineColor(color[j][k]);
		h[i][j][k]->SetMarkerColor(color[j][k]);
	}
	

	TCanvas *cv = new TCanvas("CV", "XY", 1200, 900);
	TLegend *lg[2];
	for (i=0; i<2; i++) lg[i] = new TLegend(0.3, 0.2, 0.75, 0.35);
	for (j=0; j<3; j++) for (k=0; k<2; k++) {
		sprintf(str, "%s: %s", position[j], (k) ? "Cosmic" : "IBD");
		lg[k]->AddEntry(h[0][j][k], str, "LP");
	}
	TText txt;
	txt.SetTextSize(0.08);
	cv->Divide(2, 1);
	for (i=0; i<2; i++) {
		cv->cd(i+1);
		h[i][2][0]->SetMinimum(0);
		h[i][2][0]->DrawCopy();
		for (j=0; j<3; j++) for (k=0; k<2; k++) h[i][j][k]->DrawCopy("same");
		lg[i]->Draw();
		txt.DrawTextNDC(0.45, 0.92, (i) ? "Y" : "X");
	}
		
	cv->Update();
	fRoot->Close();
}
