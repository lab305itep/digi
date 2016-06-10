void draw_cme(void)
{
	const char *fname[6] = {"cm_05cm.root", "cm_25cm.root", "cm_45cm.root", "cm_65cm.root", "cm_85cm.root", "cm_back.root"};
	const char *label[6] = {"L =  5 cm", "L = 25 cm", "L = 45 cm", "L = 65 cm", "L = 85 cm", "No 248Cm"};
	const float tm[6] = {2787, 2701, 2847, 2738, 2766, 4760};
	TFile *f[6];
	TTree *t[6];
	TH1D  *h[6];
	int i;
	char strs[64];
	TLegend *lg;

	gStyle->SetOptStat(0);
	for (i=0; i<6; i++) {
		f[i] = new TFile(fname[i]);
		t[i] = (TTree *) f[i]->Get("DanssCm");
		sprintf(strs, "h%d", i+1);
		h[i] = new TH1D(strs, "Gd neutron capture energy mesured for different 248Cm source position;MeV;Hz", 50, 0, 10);
		h[i]->SetLineColor(i+1);
		t[i]->Project(strs, "(SiPmCleanEnergy[1] + PmtCleanEnergy[1])/2", "N>1");
		h[i]->Sumw2();
		h[i]->Scale(1/tm[i]);
		printf("%s - %s\n", fname[i], label[i]);
	}
	
	h[2]->DrawCopy();
	for (i=0; i<6; i++) if (i != 2) h[i]->DrawCopy("same");

	lg = new TLegend(0.15, 0.6, 0.3, 0.9);
	for (i=0; i<6; i++) lg->AddEntry(h[i], label[i], "l");
	lg->Draw();
	gPad->Update();
	for (i=0; i<6; i++) f[i]->Close();
}

void draw_cmn(void)
{
	const char *fname[6] = {"cm_05cm.root", "cm_25cm.root", "cm_45cm.root", "cm_65cm.root", "cm_85cm.root", "cm_back.root"};
	const char *label[6] = {"L =  5 cm", "L = 25 cm", "L = 45 cm", "L = 65 cm", "L = 85 cm", "No 248Cm"};
	const float tm[6] = {2787, 2701, 2847, 2738, 2766, 4760};
	TFile *f[6];
	TTree *t[6];
	TH1D  *h[6];
	int i;
	char strs[64];
	TLegend *lg;

	gStyle->SetOptStat(0);
	for (i=0; i<6; i++) {
		f[i] = new TFile(fname[i]);
		t[i] = (TTree *) f[i]->Get("DanssCm");
		sprintf(strs, "h%d", i+1);
		h[i] = new TH1D(strs, "Gd neutron capture number of hits for different 248Cm source position;;Hz", 20, 0, 20);
		h[i]->SetLineColor(i+1);
		t[i]->Project(strs, "Hits[1]", "N>1");
		h[i]->Sumw2();
		h[i]->Scale(1/tm[i]);
		printf("%s - %s\n", fname[i], label[i]);
	}
	
	h[2]->DrawCopy();
	for (i=0; i<6; i++) if (i != 2) h[i]->DrawCopy("same");

	lg = new TLegend(0.15, 0.6, 0.3, 0.9);
	for (i=0; i<6; i++) lg->AddEntry(h[i], label[i], "l");
	lg->Draw();
	gPad->Update();
	for (i=0; i<6; i++) f[i]->Close();
}

void draw_cmt(void)
{
	const char *fname[6] = {"cm_05cm.root", "cm_25cm.root", "cm_45cm.root", "cm_65cm.root", "cm_85cm.root", "cm_back.root"};
	const char *label[6] = {"L =  5 cm", "L = 25 cm", "L = 45 cm", "L = 65 cm", "L = 85 cm", "No 248Cm"};
	const float tm[6] = {2787, 2701, 2847, 2738, 2766, 4760};
	TFile *f[6];
	TTree *t[6];
	TH1D  *h[6];
	int i;
	char strs[64];
	TLegend *lg;

	gStyle->SetOptStat(0);
	for (i=0; i<6; i++) {
		f[i] = new TFile(fname[i]);
		t[i] = (TTree *) f[i]->Get("DanssCm");
		sprintf(strs, "h%d", i+1);
		h[i] = new TH1D(strs, "Time to the first Gd neutron capture for different 248Cm source position;us;Hz", 70, 0, 70);
		h[i]->SetLineColor(i+1);
		t[i]->Project(strs, "gtDiff[1]/125", "N>1");
		h[i]->Sumw2();
		h[i]->Scale(1/tm[i]);
		printf("%s - %s\n", fname[i], label[i]);
	}
	
	h[2]->DrawCopy();
	for (i=0; i<6; i++) if (i != 2) h[i]->DrawCopy("same");

	lg = new TLegend(0.15, 0.6, 0.3, 0.9);
	for (i=0; i<6; i++) lg->AddEntry(h[i], label[i], "l");
	lg->Draw();
	gPad->Update();
	for (i=0; i<6; i++) f[i]->Close();
}

void draw_cmnn(void)
{
	const char *fname[6] = {"cm_05cm.root", "cm_25cm.root", "cm_45cm.root", "cm_65cm.root", "cm_85cm.root", "cm_back.root"};
	const char *label[6] = {"L =  5 cm", "L = 25 cm", "L = 45 cm", "L = 65 cm", "L = 85 cm", "No 248Cm"};
	const float tm[6] = {2787, 2701, 2847, 2738, 2766, 4760};
	TFile *f[6];
	TTree *t[6];
	TH1D  *h[6];
	int i;
	char strs[64];
	TLegend *lg;

	gStyle->SetOptStat(0);
	for (i=0; i<6; i++) {
		f[i] = new TFile(fname[i]);
		t[i] = (TTree *) f[i]->Get("DanssCm");
		sprintf(strs, "h%d", i+1);
		h[i] = new TH1D(strs, "Number of Gd neutron capture events from one fission for different 248Cm source position;;Hz", 10, 0, 10);
		h[i]->SetLineColor(i+1);
		t[i]->Project(strs, "N-1");
		h[i]->Sumw2();
		h[i]->Scale(1/tm[i]);
		printf("%s - %s\n", fname[i], label[i]);
	}
	
	h[2]->DrawCopy();
	for (i=0; i<6; i++) if (i != 2) h[i]->DrawCopy("same");

	lg = new TLegend(0.15, 0.6, 0.3, 0.9);
	for (i=0; i<6; i++) lg->AddEntry(h[i], label[i], "l");
	lg->Draw();
	gPad->Update();
	for (i=0; i<6; i++) f[i]->Close();
}

void draw_cmxyz(char *fname, int num=1)
{
	const char *label[3] = {"X", "Y", "Z"};	
	TFile *f;
	TTree *t;
	TH1D  *h[3];
	int i;
	char strs[64], strl[256], str[256];
	TLegend *lg;

	gStyle->SetOptStat(1);
	f = new TFile(fname);
	if (!f->IsOpen()) {
		printf("File not opened: %s\n", fname);
		return;
	}
	t = (TTree *) f->Get("DanssCm");
	if (!t) {
		printf("No tree in the file\n");
		return;
	}
	for (i=0; i<3; i++) {
		sprintf(strs, "h%s", label[i]);
		h[i] = new TH1D(strs, "Neutron center position", 25, 0, 100);
		h[i]->SetLineColor(i+2);
		sprintf(strl, "NeutronX[%d][%d]", num, i);
		sprintf(str, "N>%d && %s > 0", num, strl);
		t->Project(strs, strl, str);
		h[i]->Sumw2();
	}
	h[2]->DrawCopy();
	h[0]->DrawCopy("sames");
	h[1]->DrawCopy("sames");

	lg = new TLegend(0.15, 0.6, 0.3, 0.9);
	for (i=0; i<3; i++) lg->AddEntry(h[i], label[i], "l");
	lg->Draw();
	gPad->Update();
	f->Close();
}


