#define NRANDOM	16

class HPainter {
private:
	TFile *fSig;
	TFile *fRand;
	TTree *tSig;
	TTree *tRand;
	float upTime;
public:
	HPainter(char *base);
	HPainter(char *sname, char *rname);
	~HPainter(void);
	void Init(char *sname, char *rname);
	inline int IsOpen(void) { return tSig && tRand; };
	void Project(TH1 *hist, char *what, TCut &cut);
};

HPainter::HPainter(char *base)
{
	char strs[2048], strr[2048];

	snprintf(strs, sizeof(strs), "sel_%s.root", base);
	snprintf(strr, sizeof(strr), "sel_%s_random.root", base);

	Init(strs, strr);
}

HPainter::HPainter(char *sname, char *rname)
{
	Init(sname, rname);
}

HPainter::Init(char *sname, char *rname)
{
	TTree *info;

	tSig = tRand = NULL;
	upTime = 0;
	
	fSig = new TFile(sname);
	if (!fSig->IsOpen()) {
		printf("No file %s\n", sname);
		return;
	}
	
	fRand = new TFile(rname);
	if (!fRand->IsOpen()) {
		printf("No file %s\n", rname);
		return;
	}
	
	tSig = (TTree *) fSig->Get("DanssPair");
	tRand = (TTree *) fRand->Get("DanssPair");
	
	info = (TTree *) fSig->Get("SumInfo");
	info->GetEntry(0);
	upTime = info->GetLeaf("gTime")->GetValue() / 125E6;
	gROOT->cd();
}

HPainter::~HPainter(void)
{
	fSig->Close();
	fRand->Close();
}

void HPainter::Project(TH1 *hist, char *what, TCut &cut)
{
	TH1 *hSig;
	TH1 *hRand;
	
	if (!IsOpen()) return;
	
	hSig = (TH1 *) hist->Clone("_sigtmp");
	hRand = (TH1 *) hist->Clone("_randtmp");
	
	tSig->Project("_sigtmp", what, cut);
	tRand->Project("_randtmp", what, cut);
	
	hSig->Sumw2();
	hRand->Sumw2();
	hRand->Scale(1.0/NRANDOM);
	
	hist->Add(hSig, hRand, 1.0, -1.0);
	hist->Scale(1000.0 / upTime);	// mHz
	
	delete hSig;
	delete hRand;
}
