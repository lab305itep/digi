struct {
	double *r;
	double *cnt;
	double *ecnt;
	int N;
} DataArray;

TMinuit *MyMinuit = NULL;
// par[0] - Constant
// par[1] - shift
// par[2] - effective size
// par[3..] - scaling relative to the first (lowest) section
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
	N = 3 * DataArray.N;
	
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
			fun = par[0] / ((x[i] - par[1]) * (x[i] - par[1]) - par[2] * par[2] / 4);
			if ((i/3) != 0) {
				fun *= par[i/3 + 2];
			}
			chi2 += (fun - y[i]) * (fun - y[i]) / (ey[i] * ey[i]);
		}
		f = chi2;
	}
}

void fit_and_draw(void)
{
	double *er;
	double *ccnt;
	double *eccnt;
	double *dcnt;
	TGraphErrors *gr;
	TGraphErrors *grd;
	TH1D *hst;
	TH1D *hstd;
	double C, shift, size;
	double eC, eshift, esize;
	double *eff;
	double *eeff;
	double fmin, fedm, errdef;
	int npari, nparx;
	int i, irc;
	char str[1024];
	TF1 *fR2;
	TVirtualPad *pd;
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	MyMinuit = new TMinuit(DataArray.N + 2);
	MyMinuit->SetFCN(chi2fun);
	MyMinuit->DefineParameter(0, "Const", 1000, 10, 0, 1E10);
	MyMinuit->DefineParameter(1, "Shift", 0, 0.1, -10, 10);
	MyMinuit->DefineParameter(2, "Size", 3, 1, 0, 10);
	for (i=0; i<DataArray.N - 1; i++) {
		sprintf(str, "Eff%d", i+2);
		MyMinuit->DefineParameter(i+3, str, 1, 0.1, 0.5, 1.5);
	}
	MyMinuit->FixParameter(1);
	MyMinuit->FixParameter(2);
	MyMinuit->Migrad();
	MyMinuit->GetParameter(0, C, eC);
	MyMinuit->GetParameter(1, shift, eshift);
	MyMinuit->GetParameter(2, size, esize);
	eff  = (double *) malloc((DataArray.N - 1) * sizeof(double));
	eeff = (double *) malloc((DataArray.N - 1) * sizeof(double));
	for (i=0; i<DataArray.N - 1; i++) MyMinuit->GetParameter(i+3, eff[i], eeff[i]);
	MyMinuit->mnstat(fmin, fedm, errdef, npari, nparx, irc);

	er    = (double *) malloc(3 * DataArray.N * sizeof(double));
	memset(er, 0, 3 * DataArray.N * sizeof(double));
	ccnt  = (double *) malloc(3 * DataArray.N * sizeof(double));
	eccnt = (double *) malloc(3 * DataArray.N * sizeof(double));
	dcnt  = (double *) malloc(3 * DataArray.N * sizeof(double));

//		Renorm top and bottom sections to the middle
	for (i=0; i < 3 * DataArray.N; i++) {
		if (i < 3) {
			ccnt[i]  = DataArray.cnt[i];
			eccnt[i] = DataArray.ecnt[i];
		} else {
			ccnt[i]  = DataArray.cnt[i]  / eff[i/3-1];
			eccnt[i] = DataArray.ecnt[i] / eff[i/3-1];
		}
	}
//	for (i=0; i<9; i++) printf("L = %8.3f    CNT = %6.1f +- %4.1f\n", r[i], cnt[i], ecnt[i]);
	gr = new TGraphErrors(3 * DataArray.N, DataArray.r, ccnt, er, eccnt);
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
	hst->SetMinimum(DataArray.cnt[0] * 0.6);
	hst->SetMaximum(DataArray.cnt[0] * 1.2);
	hst->GetXaxis()->SetLabelSize(0.06);
	hst->GetYaxis()->SetLabelSize(0.06);
	hst->GetXaxis()->SetTitleSize(0.06);
	hst->GetYaxis()->SetTitleSize(0.06);
	hst->GetYaxis()->SetTitleOffset(1.25);

	hstd = new TH1D("H", ";Distance to reactor core center, m;#Delta Events per day", 35, 10, 13.5);
	hstd->SetMinimum(-25);
	hstd->SetMaximum(25);
	hstd->GetXaxis()->SetLabelSize(0.06);
	hstd->GetYaxis()->SetLabelSize(0.06);
	hstd->GetXaxis()->SetTitleSize(0.06);
	hstd->GetYaxis()->SetTitleSize(0.06);
	hstd->GetYaxis()->SetTitleOffset(1.25);
	
	TCanvas *cv = new TCanvas("CV", "R2", 800, 1200);
	cv->Divide(1, 2);
//		Do common fit and draw
	pd = cv->cd(1);
	pd->SetLeftMargin(0.15);
	pd->SetBottomMargin(0.15);
	pd->SetTopMargin(0.03);
	hst->Draw();
	gr->Draw("p");
	fR2->Draw("same");
	TLatex txt;
	txt.SetTextSize(0.07);
	sprintf(str, "#chi^{2}/n.d.f. = %6.2f/%d", fmin, 2 * DataArray.N);
	txt.DrawLatex(11.3, DataArray.cnt[0] * 1.05, str);
//		Draw difference
	for (i=0; i < 3 * DataArray.N; i++) dcnt[i] = ccnt[i] - fR2->Eval(DataArray.r[i]);
	grd = new TGraphErrors(3 * DataArray.N, DataArray.r, dcnt, er, eccnt);
	grd->SetLineColor(kBlue);
	grd->SetLineWidth(4);
	grd->SetMarkerStyle(20);
	grd->SetMarkerColor(kBlue);
	grd->SetMarkerSize(2);
	
	pd = cv->cd(2);
	pd->SetLeftMargin(0.15);
	pd->SetBottomMargin(0.15);
	pd->SetTopMargin(0.03);
	hstd->Draw();
	grd->Draw("p");
}

void fill_and_fit(int nSect, double e_min, double e_max, int mask)
{
	const char name_pattern[] = "danss_report_v4_sep17-sect%d_of_%d-calc.root";
	const char pos[3][20] = {"hUp_%d", "hMid_%d", "hDown_%d"};
	char fname[1024];
	char cpos[32];
	int i, j;
	TFile *f;
	TH1* h;
	double val, err, r;
	const double bz[3] = {11.2, 12.2, 13.2};
	const double dz3[3] = {0.212, 0.502, 0.799};
	const double dz5[5] = {0.141, 0.317, 0.501, 0.685, 0.866};
//	const double dz3[3] = {0.195, 0.5, 0.805};
//	const double dz5[5] = {0.132, 0.316, 0.5, 0.684, 0.868};
	
	if (nSect != 3 && nSect != 5) {
		printf("nSect = %d is not supported. 3 and 5 only for now.\n", nSect);
	}
	
	DataArray.N = nSect;
	DataArray.r    = (double *) malloc(3 * DataArray.N * sizeof(double));
	DataArray.cnt  = (double *) malloc(3 * DataArray.N * sizeof(double));
	DataArray.ecnt = (double *) malloc(3 * DataArray.N * sizeof(double));
	
	for (i=0; i<nSect; i++) {
		sprintf(fname, name_pattern, i+1, nSect);
		f = new TFile(fname);
		if (!f->IsOpen()) return;
		for (j=0; j<3; j++) {
			switch (nSect) {
			case 3:
				r = bz[j] - dz3[i];
				break;
			case 5:
				r = bz[j] - dz5[i];
				break;
			}
			sprintf(cpos, pos[j], mask);
			h = (TH1*) f->Get(cpos);
			if (!h) {
				printf("Something is wrong: %d %d %s.\n", i, j, cpos);
				return;
			}
			val = h->IntegralAndError(h->FindBin(e_min), h->FindBin(e_max), err);
			DataArray.r[3*i+j] = r;
			DataArray.cnt[3*i+j] = val;
			DataArray.ecnt[3*i+j] = err;
			printf("%d %d r = %f: %f +- %f\n", i, j, r, val, err);
		}
		f->Close();
	}
	fit_and_draw();
}
