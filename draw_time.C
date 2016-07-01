void draw_time(char *base, TCut cut)
{
	TH1D *h[2];
	char str[64];
	int i;
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);
	gStyle->SetOptFit();

	HPainter *p = new HPainter(base);
	if (!p->IsOpen()) {
		printf("Something wrong with base=%s.\n", base);
		return;
	}

	for (i=0; i<2; i++) {
		sprintf(str, "hT%d", i);
		h[i] = new TH1D(str, "Time between positron and neutron;T, us;mHz", 50, 0, 50);
		h[i]->SetLineWidth(4);
	}
	
	printf("Histograms are created\n");
	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
	// cut edges - more muon background there
	TCut cex("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cey("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cez("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut chot("!(PositronX[2] >= 80 && PositronX[2] < 81 && PositronX[1] >= 24 && PositronX[1] < 28)");
	// Distance cut - nothing beyond 100 cm in R and +-40cm in RZ
	TCut cr("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cg06("gtDiff > 0.6");
	TCut cgamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");

	cf = cg06 && cr && cgamma && cex && cey && cez && chot && cut;
	
	p->Project(h[0], "gtDiff", cs && cf);
	p->Project(h[1], "gtDiff", (!cs) && cf);
	
	printf("Projections are done\n");
	
	h[0]->Add(h[1], -0.05);
	
	TF1 *fdec = new TF1("FDEC", "[0]*(exp(-x/[1]) - exp(-x/[2]))", 3);
	fdec->SetParNames("Const", "t_{CAPTURE}", "t_{THERM}");
	fdec->SetParameters(h[0]->Integral()/4, 15, 5);
	h[0]->Fit("FDEC", "", "", 1, 50);
}

