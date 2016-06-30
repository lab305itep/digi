void hits5m(int first, int last)
{
	int i;
	char str[1024];

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(10);

	TChain *ch = new TChain("DanssEvent");
	for (i = first; i <= last; i++) {
		sprintf(str, "danss_root/danss_%6.6d.root", i);
		ch->AddFile(str);
	}
	TH2D *hXZ = new TH2D("hXZ", "XZ-distribution of positron vertices with energy above 5 MeV;X, cm;Z, cm", 25, 0, 100, 100, 0, 100);
	TH2D *hYZ = new TH2D("hYZ", "YZ-distribution of positron vertices with energy above 5 MeV;Y, cm;Z, cm", 25, 0, 100, 100, 0, 100);

	ch->Project("hXZ", "PositronX[2]:PositronX[0]", "PositronSiPmEnergy > 5 && PositronX[0] >= 0");
	ch->Project("hYZ", "PositronX[2]:PositronX[1]", "PositronSiPmEnergy > 5 && PositronX[1] >= 0");

	TCanvas *cv = new TCanvas("CV", "CV", 1200, 600);
	cv->Divide(2,1);
	cv->cd(1);
	hXZ->Draw("colz");
	cv->cd(2);
	hYZ->Draw("colz");
}

