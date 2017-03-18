void many_plots2(int first = 5808, int last = 11688, int mask = 0x801E)
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
	
	TFile *fRoot = new TFile("many_plots.root", "RECREATE");

	HPainter *p = new HPainter(mask, first, last);
	p->SetFile(fRoot);
	if (!p->IsOpen()) {
		printf("Something wrong with data files.\n");
		return;
	}
	
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cXYZ = cX && cY && cZ;
	TCut cR("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cT10("gtDiff > 1");
	TCut cT20("gtDiff > 2");
	TCut cT200("gtDiff < 20");						// strong cut
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cXY("PositronX[0]>=0 && PositronX[1]>=0");
	TCut cN4("NeutronEnergy > 4 && NeutronHits >= 5");

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
		h[5][i] = new TH1D(str, "Positron vertex Z;Z, cm;mHz", 100, 0, 100);
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
	
	p->Project(h[0][0], "Distance", cVeto && cIso && cT20 && cGamma && cPe && cXYZ);
	p->Project(h[0][1], "Distance", !cVeto && cIso && cT20 && cGamma && cPe && cXYZ);
	printf("Distance.\n");
	p->Project(h[1][0], "DistanceZ", cVeto && cIso && cT20 && cGamma && cPe && cXYZ);
	p->Project(h[1][1], "DistanceZ", !cVeto && cIso && cT20 && cGamma && cPe && cXYZ);
	printf("DistanceZ.\n");
	p->Project(h[2][0], "gtDiff", cVeto && cIso && cT10 && cGamma && cPe && cXYZ && cR);
	p->Project(h[2][1], "gtDiff", !cVeto && cIso && cT10 && cGamma && cPe && cXYZ && cR);
	printf("gtDiff.\n");
	p->Project(h[3][0], "PositronX[0]+2", cVeto && cIso && cT20 && cGamma && cPe && cY && cZ && cR && "PositronX[0] >= 0");
	p->Project(h[3][1], "PositronX[0]+2", !cVeto && cIso && cT20 && cGamma && cPe && cY && cZ && cR && "PositronX[0] >= 0");
	printf("X.\n");
	p->Project(h[4][0], "PositronX[1]+2", cVeto && cIso && cT20 && cGamma && cPe && cX && cZ && cR && "PositronX[1] >= 0");
	p->Project(h[4][1], "PositronX[1]+2", !cVeto && cIso && cT20 && cGamma && cPe && cX && cZ && cR && "PositronX[1] >= 0");
	printf("Y.\n");
	p->Project(h[5][0], "PositronX[2]+0.5", cVeto && cIso && cT20 && cGamma && cPe && cX && cY && cR);
	p->Project(h[5][1], "PositronX[2]+0.5", !cVeto && cIso && cT20 && cGamma && cPe && cX && cY && cR);
	printf("Z.\n");
	p->Project(h[6][0], "NeutronEnergy", cVeto && cIso && cT20 && cGamma && cPe && cY && cZ && cR);
	p->Project(h[6][1], "NeutronEnergy", !cVeto && cIso && cT20 && cGamma && cPe && cY && cZ && cR);
	printf("NE.\n");
	p->Project(h[7][0], "NeutronHits", cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR);
	p->Project(h[7][1], "NeutronHits", !cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR);
	printf("NN.\n");
	p->Project(h[8][0], "AnnihilationEnergy", cVeto && cIso && cT20 && cPe && cXYZ && cR);
	p->Project(h[8][1], "AnnihilationEnergy", !cVeto && cIso && cT20 && cPe && cXYZ && cR);
	printf("AE.\n");
	p->Project(h[9][0], "AnnihilationGammas", cVeto && cIso && cT20 && cPe && cXYZ && cR);
	p->Project(h[9][1], "AnnihilationGammas", !cVeto && cIso && cT20 && cPe && cXYZ && cR);
	printf("AG.\n");
	p->Project(h[10][0], "PositronEnergy", cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR);
	p->Project(h[10][1], "PositronEnergy", !cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR);
	p->Project(h[11][0], "PositronEnergy", cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR && cN4 && cT200);
	p->Project(h[11][1], "PositronEnergy", !cVeto && cIso && cT20 && cGamma && cPe && cXYZ && cR && cN4 && cT200);
	
	printf("Projections are done\n");
	
	for (i=0; i<12; i++) {
		for (j=0; j<2; j++) h[i][j]->SetLineWidth(4);
		h[i][0]->SetLineColor(kGreen);
		h[i][1]->SetLineColor(kRed);
		h[i][0]->SetMinimum(0);
		h[i][1]->SetMinimum(0);
	}

//	h[1][0]->Fit("gaus", "0");
//	h[1][1]->Fit("gaus", "0");
	
//	TF1 *fdec = new TF1("FDEC", "[0]*(exp(-x/[1]) - exp(-x/[2]))", 3);
//	fdec->SetParNames("Const", "t_{CAPTURE}", "t_{THERM}");
//	fdec->SetParameters(h[2][0]->Integral()/4, 15, 5);
//	h[2][0]->Fit("FDEC", "0", "", 2, 50);

	for (i=0; i<2; i++) {
		sprintf(str, "CV%d", i);
		cv[i] = new TCanvas(str, "Plots", 1800, 1200);
		cv[i]->Divide(3, 2);
		for (j=0; j<6; j++) {
			cv[i]->cd(j+1);
			h[6*i+j][0]->Draw();
			h[6*i+j][1]->Draw("sames");
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
				st->SetY2NDC(y);
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
			if (h[6*i+j][0]->GetMaximum() > h[6*i+j][1]->GetMaximum()) {
				h[6*i+j][0]->Draw();
				h[6*i+j][1]->Draw("sames");
			} else {
				h[6*i+j][1]->Draw();
				h[6*i+j][0]->Draw("sames");
			}
			gPad->Update();
		}
		cv[i]->Update();
	}

	fRoot->cd();
	for (i=0; i<12; i++) for (j=0; j<2; j++) h[i][j]->Write();
	fRoot->Close();
}
