#include "HPainter.h"

TH1D *spectr5(const char *prefix, int mask, int run_from, int run_to, double bgnd)
{
	char str[256];
	TFile *fRoot;
	
	gStyle->SetOptStat(1001100);
//		Set cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cPe;
	TCut cSig = cSel && cVeto && cIso;
	TCut cBgnd = cSel && (!cVeto);

	fRoot = new TFile ("danss_report.root", "UPDATE");

	sprintf(str, "%s_hSig", prefix);
	TH1D *hSig  = new TH1D(str,  "Positron Energy;MeV;mHz", 35, 1, 8);
	sprintf(str, "%s_hCosm", prefix);
	TH1D *hBgnd = new TH1D(str, "Positron Energy;MeV;mHz", 35, 1, 8);
	sprintf(str, "%s_hRes", prefix);
	TH1D *hRes  = new TH1D(str,  "Positron Energy;MeV;mHz", 35, 1, 8);
	
	HPainter *ptr = new HPainter(mask, run_from, run_to);
	ptr->SetFile(fRoot);
	ptr->Project(hSig,  "PositronEnergy", cSig);
	ptr->Project(hBgnd, "PositronEnergy", cBgnd);

	hRes->Add(hSig, hBgnd, 1, -bgnd);
	fRoot->cd();
	hSig->Write();
	hBgnd->Write();
	hRes->Write();
	fRoot->Close();
	delete ptr;
	return hRes;
}

void spectr_all(void)
{
	const int mask = 0x801E;
	const struct {
		char name[32];
		int first;
		int last;
	} positions[] = {
		{ "raised_30.09.16", 5540, 5807},
		{ "raised_04.10.16", 5808, 5903},	// veto corners on
		{ "mid_05.10.16",    5907, 5995},
		{ "up_10.10.16",     6007, 6130},
		{ "mid_12.10.16",    6136, 6275},
		{ "down_14.10.16",   6278, 6467},
		{ "up_17.10.16",     6469, 6570},
		{ "mid_21.10.16",    6573, 6582},
		{ "down_21.10.16",   6587, 6745},
		{ "up_24.10.16",     6757, 6815},
		{ "mid_27.10.16",    6842, 6923},
		{ "down_28.10.16",   6926, 7095},
		{ "up_31.10.16",     7106, 7364},
		{ "mid_11.11.16",    7387, 7406},
		{ "down_11.11.16",   7418, 7458},
		{ "up_14.11.16",     7478, 7579},
		{ "mid_16.11.16",    7581, 7717},
		{ "down_18.11.16",   7727, 7913},
		{ "up_21.11.16",     7922, 8042},
		{ "mid_23.11.16",    8048, 8179},
		{ "down_25.11.16",   8185, 8353},
		{ "up_28.11.16",     8357, 8430},
		{ "mid_01.12.16",    8470, 8571},
		{ "down_02.12.16",   8574, 8738},
		{ "up_05.12.16",     8741, 8869},
		{ "mid_07.12.16",    8873, 9009},
		{ "up_12.12.16",     9012, 9112},
		{ "mid_14.12.16",    9116, 9245},
		{ "down_16.12.16",   9253, 9470},
		{ "up_19.12.16",     9475, 9600},
		{ "mid_21.12.16",    9603, 9712},
		{ "down_23.12.16",   9715, 9869},
		{ "up_26.12.16",     9871, 10019},
		{ "mid_28.12.16",    10021, 10171},
		{ "down_30.12.16",   10175, 10307},
		{ "down_02.01.17",   10308, 10356},
		{ "down_03.01.17",   10357, 10424},
		{ "up_04.01.17",     10433, 10832},
		{ "mid_11.01.17",    10834, 10973},
		{ "down_13.01.17",   10979, 11147},
		{ "up_16.01.17",     11150, 11267},
		{ "mid_18.01.17",    11271, 11401},
		{ "down_20.01.17",   11404, 11563},
		{ "up_23.01.17",     11570, 11619}
	};
	int i;
	int N;
	
	N = sizeof(positions) / sizeof(positions[0]);
	for (i=0; i<N; i++) 
		spectr5(positions[i].name, mask, positions[i].first, positions[i].last, (!i) ? 0.05 : 0.025);
}

