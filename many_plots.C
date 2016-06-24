void many_plots(char *base)
{
	TH1D *h[12][2];
	int i, j;
	char str[64];
	TCanvas *cv[2];
	TPaveStats *st;
	float y, dy;
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);
	gStyle->SetOptFit();

	HPainter *p = new HPainter(base);
	if (!p->IsOpen()) {
		printf("Something wrong with base=%s.\n", base);
		return;
	}
	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
	TCut cg06("gtDiff > 0.6");									// dead time
	// cut edges - more muon background there
	TCut cx("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cy("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cz("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut cg20("gtDiff > 2 && gtDiff < 20");									// strong cut
	TCut cxy("PositronX[0]>=0 && PositronX[1]>=0");
	TCut cn5("NeutronEnergy > 5 && NeutronHits >= 5");
	TCut cgamma("AnnihilationEnergy < 1.2 && AnnihilationGammas > 0");
	TF1 *fdec = new TF1("FDEC", "expo(0) - expo(2)", 4);

	for (i=0; i<2; i++) {
		sprintf(str, "hR%d", i);
		h[0][i] = new TH1D(str, "Distance between positron and neutron;R, cm;mHz", 40, 0, 160);
		sprintf(str, "hRZ%d", i);
		h[1][i] = new TH1D(str, "Distance between positron and neutron, projection Z;R_{z}, cm;mHz", 100, -100, 100);
		sprintf(str, "hT%d", i);
		h[2][i] = new TH1D(str, "Time between positron and neutron;T, us;mHz", 50, 0, 50);
		sprintf(str, "hX%d", i);
		h[3][i] = new TH1D(str, "Positron vertex X;X, cm;mHz", 25, 0, 100);
		sprintf(str, "hY%d", i);
		h[4][i] = new TH1D(str, "Positron vertex Y;Y, cm;mHz", 25, 0, 100);
		sprintf(str, "hZ%d", i);
		h[5][i] = new TH1D(str, "Positron vertex Z;Z, cm;mHz", 25, 0, 100);
		sprintf(str, "hNE%d", i);
		h[6][i] = new TH1D(str, "Energy detected in neutron capture;E_{n}, MeV;mHz", 50, 0, 10);
		sprintf(str, "hNN%d", i);
		h[7][i] = new TH1D(str, "Number of hits in SiPM for neutron capture;N_{n};mHz", 20, 0, 20);
		sprintf(str, "hGE%d", i);
		h[8][i] = new TH1D(str, "Energy beyond positron cluster;E_{#gamma}, MeV;mHz", 15, 0, 3);
		sprintf(str, "hGN%d", i);
		h[9][i] = new TH1D(str, "Number of SiPM hits out of positron cluster;N_{#gamma};mHz", 10, 0, 10);
		sprintf(str, "hE%d", i);
		h[10][i] = new TH1D(str, "Positorn kinetic energy;E, MeV;mHz", 40, 0, 8);
		sprintf(str, "hEC%d", i);
		h[11][i] = new TH1D(str, "Positorn kinetic energy, strong background cuts;E, MeV;mHz", 40, 0, 8);
	}
	
	printf("Histograms are created\n");
	
	p->Project(h[0][0], "Distance", cs && cg06);
	p->Project(h[0][1], "Distance", (!cs) && cg06);
	printf("Distance.\n");
	p->Project(h[1][0], "DistanceZ", cs && cg06);
	p->Project(h[1][1], "DistanceZ", (!cs) && cg06);
	printf("DistanceZ.\n");
	p->Project(h[2][0], "gtDiff", cs && cg06);
	p->Project(h[2][1], "gtDiff", (!cs) && cg06);
	printf("gtDiff.\n");
	p->Project(h[3][0], "PositronX[0]", cs && cg06);
	p->Project(h[3][1], "PositronX[0]", (!cs) && cg06);
	printf("X.\n");
	p->Project(h[4][0], "PositronX[1]", cs && cg06);
	p->Project(h[4][1], "PositronX[1]", (!cs) && cg06);
	printf("Y.\n");
	p->Project(h[5][0], "PositronX[2]", cs && cg06);
	p->Project(h[5][1], "PositronX[2]", (!cs) && cg06);
	printf("Z.\n");
	p->Project(h[6][0], "NeutronEnergy", cs && cg06);
	p->Project(h[6][1], "NeutronEnergy", (!cs) && cg06);
	printf("NE.\n");
	p->Project(h[7][0], "NeutronHits", cs && cg06);
	p->Project(h[7][1], "NeutronHits", (!cs) && cg06);
	printf("NN.\n");
	p->Project(h[8][0], "AnnihilationEnergy", cs && cg06);
	p->Project(h[8][1], "AnnihilationEnergy", (!cs) && cg06);
	printf("AE.\n");
	p->Project(h[9][0], "AnnihilationGammas", cs && cg06);
	p->Project(h[9][1], "AnnihilationGammas", (!cs) && cg06);
	printf("AG.\n");
	p->Project(h[10][0], "PositronEnergy", cs && cg06);
	p->Project(h[10][1], "PositronEnergy", (!cs) && cg06);
	p->Project(h[11][0], "PositronEnergy", cs && cg20 && cn5 && cgamma);
	p->Project(h[11][1], "PositronEnergy", (!cs) && cg20 && cn5 && cgamma);
	
	printf("Projections are done\n");
	
	for (i=0; i<12; i++) {
		for (j=0; j<2; j++) h[i][j]->SetLineWidth(4);
		h[i][0]->SetLineColor(kGreen);
		h[i][1]->SetLineColor(kRed);
		h[i][1]->SetMinimum(0);
	}

	h[1][0]->Fit("gaus", "", "0");
	h[1][1]->Fit("gaus", "", "0");
	h[2][0]->Fit("FDEC", "", "0");

	for (i=0; i<2; i++) {
		sprintf(str, "CV%d", i);
		cv[i] = new TCanvas(str, "Plots", 1800, 1200);
		cv[i]->Divide(3, 2);
		for (j=0; j<6; j++) {
			cv[i]->cd(j+1);
			h[6*i+j][1]->Draw();
			h[6*i+j][0]->Draw("sames");
			gPad->Update();
			if (i == 0 && j == 1) {
				st = (TPaveStats *) h[6*i+j][0]->FindObject("stats");
				st->SetOptStat(1100);
				st->SetLineColor(kGreen);
				st->SetTextColor(kGreen);
				y = st->GetY1NDC();
				dy = st->GetY2NDC() - y;
				st->SetX1NDC(0.72);
				st = (TPaveStats *) h[6*i+j][1]->FindObject("stats");
				st->SetOptStat(1100);
				st->SetLineColor(kRed);
				st->SetTextColor(kRed);
				st->SetX1NDC(0.72);
				st->SetY1NDC(y - dy);
			} else {
				st = (TPaveStats *) h[6*i+j][0]->FindObject("stats");
				st->SetLineColor(kGreen);
				st->SetTextColor(kGreen);
				y = st->GetY1NDC();
				dy = st->GetY2NDC() - y;
				st->SetX1NDC(0.72);
				st->SetY1NDC(y + dy/2);
				st = (TPaveStats *) h[6*i+j][1]->FindObject("stats");
				st->SetLineColor(kRed);
				st->SetTextColor(kRed);
				st->SetX1NDC(0.72);
				st->SetY2NDC(y + dy/2);
			}
			h[6*i+j][1]->Draw();
			h[6*i+j][0]->Draw("sames");
			gPad->Update();
		}
		cv[i]->Update();
	}
}

