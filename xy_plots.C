#include "HPainter.h"

void xy_plots(int first = 5808, int last = 11688)
{
	TH1D *h[2][3][2];
	int i, j, k;
	char str[64];
	TCanvas *cv[2];
	TPaveStats *st;
	float y, dy;
	int mask[3] = {2, 4, 8};	// down, middle, up
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1000000);
	gStyle->SetOptFit();
	
	TFile *fRoot = new TFile("xy_plots.root", "RECREATE");

	HPainter *p[3];
	for (i=0; i<3; i++) p[i] = new HPainter(mask[i], first, last);
	if (!p[0]->IsOpen() || !p[1]->IsOpen() || !p[2]->IsOpen()) {
		printf("Something wrong with data files.\n");
		return;
	}
	
	for (i=0; i<3; i++) p[i]->SetFile(fRoot);
	
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ = "PositronX[2] > 3.5 && PositronX[2] < 95.5";
	TCut cR("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cT20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
        TCut cSig = cIso && cR && cT20 && cGamma && cPe;

	fRoot->cd();
	for (i=0; i<2; i++) for (j=0; j<3; j++) {
		sprintf(str, "hX%d%c", j, (i) ? 'C' : 'S');
		h[0][j][i] = new TH1D(str, "Positron vertex X;X, cm;mHz", 25, 0, 100);
		sprintf(str, "hY%d%c", j, (i) ? 'C' : 'S');
		h[1][j][i] = new TH1D(str, "Positron vertex Y;Y, cm;mHz", 25, 0, 100);
	}
	
	printf("Histograms are created\n");
	for (i=0; i<3; i++) {
		p[i]->Project(h[0][i][0], "PositronX[0]+2", cVeto && cSig && cY && "PositronX[0] >= 0" && cZ);
		p[i]->Project(h[0][i][1], "PositronX[0]+2", !cVeto && cSig && cY && "PositronX[0] >= 0" && cZ);
		printf("X%d\n", i);
		p[i]->Project(h[1][i][0], "PositronX[1]+2", cVeto && cSig && cX && "PositronX[1] >= 0" && cZ);
		p[i]->Project(h[1][i][1], "PositronX[1]+2", !cVeto && cSig && cX && "PositronX[1] >= 0" && cZ);
		printf("Y%d\n", i);
	}
	printf("Projections are done\n");

	fRoot->cd();
	for (i=0; i<2; i++) for (j=0; j<3; j++) for (k=0; k<2; k++) h[i][j][k]->Write();
	fRoot->Close();
}
