void histdump(const char *fname, const char *hname)
{
	int i;
	
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) return;
	TH1 *h = (TH1 *)f->Get(hname);
	if (!h) {
		printf("%s:%s not found!\n", fname, hname);
		return;
	}
	
	printf("%s:%s (%s)\n%s\t\t%s\n", fname, hname, h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
	for (i=0; i<h->GetNbinsX(); i++) printf("%f-%f:  %f +- %f\n",
		h->GetBinLowEdge(i+1), h->GetBinLowEdge(i+1) + h->GetBinWidth(i+1), h->GetBinContent(i+1), h->GetBinError(i+1));
	f->Close();
}
