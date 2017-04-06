#include "HPainter.h"
#define NHISTS 18
void background_calc(const char *fname = "background_plots.root", int run_first = 5808, int run_last = 11688)
{
	char strs[128];
	char strl[1024];
	const char titles[NHISTS][16] = {"gtDiff", "R1", "R2", "RZ", "PX", "PY", "PZ", "NX", "NY", "NZ", "NE", "NH", "NR", "PH", "PMaxE", "AH", "AE", "PF"};
	const char titlel[NHISTS][256] = {"Time from positron to neutron;us", "Distance between positron and neutron, 2D;cm", 
		"Distance between positron and neutron, 3D;cm", "Distance between positron and neutron, Z;cm", 
		"Positron vertex X;cm", "Positron vertex Y;cm", "Positron vertex Z;cm", 
		"Neutron vertex X;cm", "Neutron vertex Y;cm", "Neutron vertex Z;cm", "Neutron capture energy", "Neutron capture SiPM hits",
		"Neutron capture average photon flight",
		"Hits in positron cluster", "Maximum energy in one hit", "Number of hits out of cluster", "Energy out of cluster;MeV",
		"Edge hit in prompt event"};
	TH1D *h[NHISTS][2];
	int i, j;
	
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR2("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
        TCut cR("Distance < 60 && fabs(DistanceZ) < 40");
        TCut ct;

	TFile *fRoot = new TFile(fname, "RECREATE");
	for (i=0; i<NHISTS; i++) for (j=0; j<2; j++) {
		sprintf(strs, "h%s%c", titles[i], (j) ? 'C' : 'N');
		sprintf(strl, "%s: %s", (j) ? "Cosmic" : "Neutrino", titlel[i]);
		switch(i) {
		case 0:		// gtDiff
			h[i][j] = new TH1D(strs, strl, 50, 0, 50.0);
			break;
		case 1:		// R1, R2
		case 2:
			h[i][j] = new TH1D(strs, strl, 40, 0, 160.0);
			break;
		case 3:		// RZ
			h[i][j] = new TH1D(strs, strl, 100, -50.0, 50.0);
			break;
		case 4:		// PX, PY, NX, NY
		case 5:
		case 7:
		case 8:
			h[i][j] = new TH1D(strs, strl, 25, 0, 100.0);
			break;
		case 6:		// PZ, NZ
		case 9:
			h[i][j] = new TH1D(strs, strl, 100, 0, 100.0);
			break;
		case 10:	// NE
			h[i][j] = new TH1D(strs, strl, 45, 3.0, 12.0);
			break;
		case 11:	// NH
			h[i][j] = new TH1D(strs, strl, 20, 0, 20.0);
			break;
		case 12:	// NR
			h[i][j] = new TH1D(strs, strl, 15, 0, 60.0);
			break;
		case 13:	// PH
			h[i][j] = new TH1D(strs, strl, 10, 0, 10.0);
			break;
		case 14:	// PMaxE
			h[i][j] = new TH1D(strs, strl, 40, 0, 8.0);
			break;
		case 15:	// AH
			h[i][j] = new TH1D(strs, strl, 20, 0, 20.0);
			break;
		case 16:	// AE
			h[i][j] = new TH1D(strs, strl, 20, 0, 4.0);
			break;
		case 17:	// PF
			h[i][j] = new TH1D(strs, strl, 2, 0, 2);
			break;
		}
	}

	HPainter *hp = new HPainter(0x801E, run_first, run_last);
	hp->SetFile(fRoot);

	for (j=0; j<2; j++) {
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[0][j], "gtDiff", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma && !cR2;
		hp->Project(h[1][j], "Distance", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma && cR2;
		hp->Project(h[2][j], "Distance", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma;
		hp->Project(h[3][j], "DistanceZ", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cY && cZ && cR && cPe && cGamma && "PositronX[0] >= 0";
		hp->Project(h[4][j], "PositronX[0] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cZ && cR && cPe && cGamma && "PositronX[1] >= 0";
		hp->Project(h[5][j], "PositronX[1] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cR && cPe && cGamma && "PositronX[2] >= 0";
		hp->Project(h[6][j], "PositronX[2] + 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cY && cZ && cR && cPe && cGamma && "NeutronX[0] >= 0";
		hp->Project(h[7][j], "NeutronX[0] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cZ && cR && cPe && cGamma && "NeutronX[1] >= 0";
		hp->Project(h[8][j], "NeutronX[1] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cR && cPe && cGamma && "NeutronX[2] >= 0";
		hp->Project(h[9][j], "NeutronX[2] + 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[10][j], "NeutronEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[11][j], "NeutronHits", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[12][j], "NeutronRadius", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[13][j], "PositronHits", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[14][j], "MaxHitEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe;
		hp->Project(h[15][j], "AnnihilationGammas", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe;
		hp->Project(h[16][j], "AnnihilationEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[17][j], "(PositronFlags[0] && 0x3FF0000) ? 1.5 : 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
	}
	
	fRoot->cd();
	for (i=0; i<NHISTS; i++) for (j=0; j<2; j++) h[i][j]->Write();
	delete hp;
	fRoot->Close();
}
