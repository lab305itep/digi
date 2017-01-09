#include <stdio.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>

#include "HPainter.h"

TH1D *spectr_sect(int sect, const char *base1, float bgnd1 = 0.05, const char *base2 = NULL, float bgnd2 = 0.05)
{
	TH1D *h[2][2];	// spectrum
	TH1D *hSum;
	int i;
	char str[64], strl[1024];
	HPainter *p[2];

	p[0] = new HPainter(base1);
	if (!p[0]->IsOpen()) {
		printf("Something wrong with base=%s.\n", base1);
		return NULL;
	}
	if (base2 && strlen(base2)) {
		p[1] = new HPainter(base2);
		if (!p[1]->IsOpen()) {
			printf("Something wrong with base=%s.\n", base2);
			return NULL;
		}
	} else {
		p[1] = NULL;
	}

	TCut cs("gtFromVeto > 100 && EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100");	// standard Veto
	// cut edges - more muon background there
	TCut cex("(PositronX[0]>=4 && PositronX[0]<=96) || PositronX[0] < 0");
	TCut cey("(PositronX[1]>=4 && PositronX[1]<=96) || PositronX[1] < 0");
	TCut cez("(PositronX[2]>=4 && PositronX[2]<=96) || PositronX[2] < 0");
	TCut chot("!(PositronX[2] >= 80 && PositronX[2] < 81 && PositronX[1] >= 24 && PositronX[1] < 28)");
	// Distance cut - nothing beyond 100 cm in R and +-40cm in RZ
	TCut cr("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut cg20("gtDiff > 2");
	TCut cgamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
	TCut cpe("PositronEnergy>1");

	TCut cf = cg20 && cr && cgamma && cex && cey && cez && chot && cpe;
	if (sect>0) {
		sprintf(strl, "PositronX[2]>=%6.1f && PositronX[2]<%6.1f", 20.0*(sect-1), 20.0*sect);
		TCut csect(strl);
		cf = cf && csect;
	}

	for (i=0; i<2; i++) {
		sprintf(str, "hS%s%d", base1, i);
		sprintf(strl, "Positron spectrum, DANSS position = %s;E, MeV;mHz", base1);
		h[0][i] = new TH1D(str, strl, 40, 0, 8);
		if (p[1]) {
			sprintf(str, "hS%s%d", base2, i);
			sprintf(strl, "Positron spectrum, DANSS position = %s;E, MeV;mHz", base2);
			h[1][i] = new TH1D(str, strl, 40, 0, 8);			
		}
	}
	printf("Histograms are created\n");
	
	for (i=0; i<2; i++) if (p[i]) {
		p[i]->Project(h[i][0], "PositronEnergy", cs && cf);
		p[i]->Project(h[i][1], "PositronEnergy", (!cs) && cf);
		h[i][0]->Add(h[i][1], -((i) ? bgnd2 : bgnd1));	// background
	}
	
	printf("Projections are done\n");
	
	sprintf(str, "h%s_%d", base1, sect);
	hSum = (TH1D *) h[0][0]->Clone(str);
	if (p[1]) {
		hSum->SetBit(TH1::kIsAverage);
		h[1][0]->SetBit(TH1::kIsAverage);
		hSum->Add(h[1][0]);
	}
	for (i=0; i<2; i++) if (p[i]) delete p[i];
	
	return hSum;
}

void spectr2(void)
{
	const char *base1[3] = {"upa", "downa", "stucka"};
	const char *base2[3] = {"upb", "downb", ""};
	const char *sname[3] = {"UP", "DOWN", "MID"};
	const char *rname[3] = {"DOWN_UP", "MID_UP", "DOWN_MID"};
	const int rats[3][2] = {{1, 0}, {2, 0}, {1, 2}};
	const float bgnd1 = 0.05;
	const float bgnd2 = 0.10;
	int i; int j;
	char strs[64], strl[1024];
	TH1D *h[3][6];
	TH1D *hr[3][6];
	
	TFile f("spectr2.root", "RECREATE");
//		Spectra	
	for (i=0; i<3; i++) for (j=0; j<6; j++) {
		h[i][j] = spectr_sect(j, base1[i], bgnd1, base2[i], bgnd2);
		if (j) {
			sprintf(strs, "spec_%s_sect%d", sname[i], j-1);
			sprintf(strl, "Positron kinetic energy. Position %s. Section %d.;E_{KIN}, MeV;F, mHz", sname[i], j-1);
		} else {
			sprintf(strs, "spec_%s", sname[i]);
			sprintf(strl, "Positron kinetic energy. Position %s.;E_{KIN}, MeV;F, mHz", sname[i]);
		}
		h[i][j]->SetName(strs);
		h[i][j]->SetTitle(strl);
		f.cd();
		h[i][j]->Write();
	}
//		Ratios	
	for (i=0; i<3; i++) for (j=0; j<6; j++) {
		if (j) {
			sprintf(strs, "ratio_%s_sect%d", rname[i], j-1);
			sprintf(strl, "Ratio %s. Section %d.;E_{KIN}, MeV;", rname[i], j-1);
		} else {
			sprintf(strs, "ratio_%s", rname[i]);
			sprintf(strl, "Ratio %s.;E_{KIN}, MeV;", rname[i]);
		}
		hr[i][j] = new TH1D(strs, strl, 40, 0, 8);
		hr[i][j]->Divide(h[rats[i][0]][j], h[rats[i][1]][j]);
		f.cd();
		hr[i][j]->Write();
	}
	f.Close();
}

