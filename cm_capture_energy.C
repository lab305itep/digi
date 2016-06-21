void cm_capture_energy(char *fname)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	TH1D *hcm = new TH1D("HDC", "Energy detected in neutron capture;E, MeV;N", 60, 0, 12);
	hcm->SetLineWidth(4);

	TFile f(fname);
	if (!f.IsOpen()) return;
	TTree *t = f.Get("DanssCm");
	if (!t) return;
	gROOT->cd();
	t->Project("HDC", "(SiPmCleanEnergy[1]+PmtCleanEnergy[1])/2", "N>1 && gtDiff[1]/125<50");
	
	new TCanvas("CV", "Neutron Capture", 600, 800);
	hcm->DrawCopy("errors");
	TText txt;
	txt.DrawText(2, 200, "Curium n source");
	f.Close();
}

