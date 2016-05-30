TH1D *ddiff(char *what, char *cut, TTree *tA, TTree *tB, int N=20, double amin=0, double amax=10, double ratio=1.0)
{
	TH1D *h = gROOT->FindObject("hDiff");
	if (h) delete h;
	h = new TH1D("hDiff", what, N, amin, amax);
	TH1D *hTmpA = new TH1D("_hTmpA", "Tmp A", N, amin, amax);
	TH1D *hTmpB = new TH1D("_hTmpB", "Tmp B", N, amin, amax);
	tA->Project("_hTmpA", what, cut);
	tB->Project("_hTmpB", what, cut);
	hTmpA->Sumw2();
	hTmpB->Sumw2();
	h->Add(hTmpA, hTmpB, 1, -ratio);
	delete hTmpA;
	delete hTmpB;
	return h;
}

