void draw_R2(double size = 0)
{
//		Original data
	double r[9] = {10.395, 11.395, 12.395, 10.7, 11.7, 12.7, 11.005, 12.005, 13.005};
	double er[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
//		???
//	double cnt[9] = {1878, 1577, 1361, 1851, 1572, 1321, 1585, 1341, 1140};
//	double ecnt[9] = {14, 17, 15, 14, 17, 14, 12, 15, 12};
//		Oct 16 - May 17
//	double cnt[9] = {1689, 1399, 1198, 1710, 1429, 1201, 1377, 1155, 958};
//	double ecnt[9] = {8, 8, 6, 8, 8, 6, 7, 7, 5};
//		All
	double cnt[9] = {1675, 1391, 1187, 1701, 1424, 1190, 1370, 1147, 952};
	double ecnt[9] = {7, 8, 6, 7, 8, 6, 6, 7, 5};
//		J-M
//	double cnt[9] = {1712, 1394, 1193, 1719, 1416, 1196, 1387, 1154, 961};
//	double ecnt[9] = {12, 11, 8, 12, 11, 8, 10, 9, 7};
//		A-M
//	double cnt[9] = {1669, 1408, 1205, 1685, 1444, 1196, 1378, 1138, 959};
//	double ecnt[9] = {14, 17, 11, 14, 17, 11, 12, 14, 9};
	double dcnt[9];

	TGraphErrors *gr;
	TGraphErrors *grd;
	TH1D *hst;
	TH1D *hstd;
	double par[3];
	int i;
	char str[1024];
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	sprintf(str, "[0]/(x * x - %f)", size);
	TF1 *fR2 = new TF1("fR2", str);
	fR2->SetParNames("Const.");
//		Fit 3 individual distributions
	for (i=0; i<3; i++) {
		gr = new TGraphErrors(3, &r[3*i], &cnt[3*i], &er[3*i], &ecnt[3*i]);
		fR2->SetParameter(0, 200000);
		gr->Fit(fR2);
		par[i] = fR2->GetParameter(0);
	}
//		Renorm top and bottom sections to the middle
	for (i=0; i<3; i++) {
		cnt[i] *= par[1] / par[0];
		cnt[6+i] *= par[1] / par[2];
		ecnt[i] *= par[1] / par[0];
		ecnt[6+i] *= par[1] / par[2];
	}
//	for (i=0; i<9; i++) printf("L = %8.3f    CNT = %6.1f +- %4.1f\n", r[i], cnt[i], ecnt[i]);
	gr = new TGraphErrors(9, r, cnt, er, ecnt);
	gr->SetLineColor(kBlue);
	gr->SetLineWidth(4);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(kBlue);
	gr->SetMarkerSize(2);
	fR2->SetLineColor(kRed);
	fR2->SetLineWidth(3);
	hst = new TH1D("H", ";Distance to reactor core center, m;Events per day", 35, 10, 13.5);
	hst->SetMinimum(900);
	hst->SetMaximum(2000);
	hst->GetXaxis()->SetLabelSize(0.06);
	hst->GetYaxis()->SetLabelSize(0.06);
	hst->GetXaxis()->SetTitleSize(0.06);
	hst->GetYaxis()->SetTitleSize(0.06);
	hst->GetYaxis()->SetTitleOffset(1.25);

	hstd = new TH1D("H", ";Distance to reactor core center, m;Events per day", 35, 10, 13.5);
	hstd->SetMinimum(-50);
	hstd->SetMaximum(50);
	hstd->GetXaxis()->SetLabelSize(0.06);
	hstd->GetYaxis()->SetLabelSize(0.06);
	hstd->GetXaxis()->SetTitleSize(0.06);
	hstd->GetYaxis()->SetTitleSize(0.06);
	hstd->GetYaxis()->SetTitleOffset(1.25);
	
//		Do common fit and draw
	TCanvas *cv = new TCanvas("CV", "R2", 1200, 900);
	cv->SetLeftMargin(0.15);
	cv->SetBottomMargin(0.12);
	hst->Draw();
	fR2->SetParameter(0, 200000);
	gr->Fit(fR2);
	gr->Draw("p");
	TLatex txt;
	txt.SetTextSize(0.07);
	sprintf(str, "#chi^{2}/n.d.f. = %6.2f/6", fR2->GetChisquare());
	txt.DrawLatex(11.3, 1800, str);
//		Draw difference
	for (i=0; i<9; i++) dcnt[i] = cnt[i] - fR2->Eval(r[i]);
	grd = new TGraphErrors(9, r, dcnt, er, ecnt);
	grd->SetLineColor(kBlue);
	grd->SetLineWidth(4);
	grd->SetMarkerStyle(20);
	grd->SetMarkerColor(kBlue);
	grd->SetMarkerSize(2);
	
	cv = new TCanvas("CVD", "R2 difference", 1200, 900);
	cv->SetLeftMargin(0.15);
	cv->SetBottomMargin(0.12);
	hstd->Draw();
	grd->Draw("p");
}

