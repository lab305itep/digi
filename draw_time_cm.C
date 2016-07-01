void draw_time_cm(char *fname)
{
	const float tm = 2701;
	TFile *f;
	TTree *t;
	TH1D  *h;
	int i;
	char strs[64];

	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	
	f = new TFile(fname);
	t = (TTree *) f->Get("DanssCm");
	h = new TH1D("hTcm", "Time to the first Gd neutron capture for 248Cm source at 25cm;T, us;Hz", 50, 0, 50);
	h->SetLineWidth(4);
//	t->Project("hTcm", "gtDiff[1]/125", "N>1 && NeutronEnergy[1]>3 && PositronX[0][0] > 40 && PositronX[0][0] < 60 && PositronX[0][1] > 10 && PositronX[0][1] < 40 && PositronX[0][2] > 40 && PositronX[0][2] < 60");	
	t->Project("hTcm", "gtDiff[1]/125", "N>1 && NeutronEnergy[1]>3 && PositronX[0][0] > 40 && PositronX[0][0] < 60 && PositronX[0][1] < 15 && PositronX[0][2] > 40 && PositronX[0][2] < 60");
//	t->Project("hTcm", "gtDiff[1]/125", "N>1 && NeutronEnergy[1]>3");
	h->Sumw2();
	h->Scale(1/tm);
	
	TF1 *fdec3 = new TF1("FDEC3", "[0]*(exp(-x/[1]) - exp(-x/[2]))*pow(([1]*exp(-x/[1]) - [2]*exp(-x/[2])), [3]-1)", 4);
	fdec3->SetParNames("Const", "t_{CAPTURE}", "t_{THERM}", "K");
	fdec3->SetParameters(h->Integral()/4, 15, 5, 3);
	fdec3->FixParameter(3, 3);
	h->Fit("FDEC3", "", "", 1, 50);	
}

