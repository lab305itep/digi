const struct {
	double r[9] = {10.412, 11.412, 12.412, 10.7, 11.7, 12.7, 10.999, 11.999, 12.999};
//		All but very old
	double cnt[9] = {1681, 1392, 1191, 1702, 1421, 1193, 1371, 1148, 952};
	double ecnt[9] = {8, 8, 6, 8, 8, 6, 7, 7, 5};
} DataArray;

TMinuit *MyMinuit = NULL;
// par[0] - Constant
// par[1,2] - Up and down scaling relative to middle
// par[3] - shift
// par[4] - effective size
static void chi2fun(int &npar, double *gin, double &f, double *par,  int  iflag)
{
	double chi2, fun;
	int i;
	int N;
	const double *x;
	const double *y;
	const double *ey;

	x = DataArray.r;
	y = DataArray.cnt;
	ey = DataArray.ecnt;
	N = 9;
	
	switch (iflag) {
	case 1: // Initialization
		break;
	case 2: //Compute derivatives
		break;
	case 3: // after the fit is finished
		break;
	default:
		chi2 = 0;
		for (i=0; i<N; i++) {
			fun = par[0] / ((x[i] - par[3]) * (x[i] - par[3]) - par[4] * par[4] / 4);
			if ((i/3) == 0) {
				fun *= par[1];
			} else if ((i/3) == 2) {
				fun *= par[2];
			}
			chi2 += (fun - y[i]) * (fun - y[i]) / (ey[i] * ey[i]);
		}
		f = chi2;
	}
}

void draw_R2m(void)
{
	const double er[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	double ccnt[9];
	double eccnt[9];
	double dcnt[9];
	TGraphErrors *gr;
	TGraphErrors *grd;
	TH1D *hst;
	TH1D *hstd;
	double effUp, effDown, C, shift, size;
	double eeffUp, eeffDown, eC, eshift, esize;
	double fmin, fedm, errdef;
	int npari, nparx;
	int i, irc;
	char str[1024];
	TF1 *fR2;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	MyMinuit = new TMinuit(5);
	MyMinuit->SetFCN(chi2fun);
	MyMinuit->DefineParameter(0, "Const", 1000, 10, 0, 1E10);
	MyMinuit->DefineParameter(1, "EffUp", 1, 0.1, 0.5, 1.5);
	MyMinuit->DefineParameter(2, "EffDown", 1, 0.1, 0.5, 1.5);
	MyMinuit->DefineParameter(3, "Shift", 0, 0.1, -10, 10);
	MyMinuit->DefineParameter(4, "Size", 3, 1, 0, 10);
	MyMinuit->FixParameter(3);
	MyMinuit->FixParameter(4);
	MyMinuit->Migrad();
	MyMinuit->GetParameter(0, C, eC);
	MyMinuit->GetParameter(1, effUp, eeffUp);
	MyMinuit->GetParameter(2, effDown, eeffDown);
	MyMinuit->GetParameter(3, shift, eshift);
	MyMinuit->GetParameter(4, size, esize);
	MyMinuit->mnstat(fmin, fedm, errdef, npari, nparx, irc);

//		Renorm top and bottom sections to the middle
	for (i=0; i<3; i++) {
		ccnt[i] = DataArray.cnt[i] / effUp;
		ccnt[i+3] = DataArray.cnt[i+3];
		ccnt[i+6] = DataArray.cnt[i+6] / effDown;
		eccnt[i] = DataArray.ecnt[i] / effUp;
		eccnt[i+3] = DataArray.ecnt[i+3];
		eccnt[i+6] = DataArray.ecnt[i+6] / effDown;
	}
//	for (i=0; i<9; i++) printf("L = %8.3f    CNT = %6.1f +- %4.1f\n", r[i], cnt[i], ecnt[i]);
	gr = new TGraphErrors(9, DataArray.r, ccnt, er, eccnt);
	gr->SetLineColor(kBlue);
	gr->SetLineWidth(4);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(kBlue);
	gr->SetMarkerSize(2);
	
	fR2 = new TF1("fR2", "[0] / ((x - [1]) * (x - [1]) - [2] * [2] / 4.0)", 1, 100);
	fR2->SetParameter(0, C);
	fR2->SetParameter(1, shift);
	fR2->SetParameter(2, size);
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
	gr->Draw("p");
	fR2->Draw("same");
	TLatex txt;
	txt.SetTextSize(0.07);
	sprintf(str, "#chi^{2}/n.d.f. = %6.2f/5", fmin);
	txt.DrawLatex(11.3, 1800, str);
//		Draw difference
	for (i=0; i<9; i++) dcnt[i] = ccnt[i] - fR2->Eval(DataArray.r[i]);
	grd = new TGraphErrors(9, DataArray.r, dcnt, er, eccnt);
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

