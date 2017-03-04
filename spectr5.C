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
	TCut cBgnd = cSel && (!cVeto) && (cIso || "gtFromPrevious == gtFromVeto");
#ifdef STRONG_CUTS
        TCut cStrong("NeutronEnergy > 4 && gtDiff < 20");
        cSig = cSig && cStrong;
        cBgnd = cBgnd && cStrong;
	fRoot = new TFile ("danss_report_strong.root", "UPDATE");
#else
	fRoot = new TFile ("danss_report.root", "UPDATE");
#endif

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
		double bgnd;
	} positions[] = {
		{ "down_21.04.16",   2307, 2361, 0.05},
		{ "up_22.04.16",     2366, 2387, 0.05},
		{ "down_23.04.16",   2399, 2445, 0.05},
		{ "mid_24.04.16",    2449, 2512, 0.05},
		{ "up_25.04.16",     2514, 2562, 0.05},
		{ "down_26.04.16",   2564, 2609, 0.05},
		{ "mid_27.04.16",    2620, 2685, 0.05},
		{ "up_28.04.16",     2687, 2730, 0.05},
		{ "down_29.04.16",   2733, 2788, 0.05},
		{ "mid_30.04.16",    2791, 2832, 0.05},
		{ "stuck_01.05.16",  2836, 3399, 0.05},
		{ "stuck_10.05.16",  3400, 3788, 0.05},
		{ "up_03.06.16",     4261, 4400, 0.1},
		{ "down_06.06.16",   4403, 4439, 0.1},
		{ "up_08.06.16",     4508, 4601, 0.1},
		{ "raised_30.09.16", 5540, 5807, 0.05},
		{ "raised_04.10.16", 5808, 5903, 0.025},	// veto corners on
		{ "mid_05.10.16",    5907, 5995, 0.025},
		{ "up_10.10.16",     6007, 6130, 0.025},
		{ "mid_12.10.16",    6136, 6275, 0.025},
		{ "down_14.10.16",   6278, 6467, 0.025},
		{ "up_17.10.16",     6469, 6570, 0.025},
		{ "mid_21.10.16",    6573, 6582, 0.025},
		{ "down_21.10.16",   6587, 6745, 0.025},
		{ "up_24.10.16",     6757, 6815, 0.025},
		{ "mid_27.10.16",    6842, 6909, 0.025},	// was 6923
		{ "down_28.10.16",   6926, 7095, 0.025},
		{ "up_31.10.16",     7106, 7364, 0.025},
		{ "mid_11.11.16",    7387, 7406, 0.025},
		{ "down_11.11.16",   7418, 7458, 0.025},
		{ "up_14.11.16",     7478, 7579, 0.025},
		{ "mid_16.11.16",    7581, 7717, 0.025},
		{ "down_18.11.16",   7727, 7913, 0.025},
		{ "up_21.11.16",     7922, 8042, 0.025},
		{ "mid_23.11.16",    8048, 8179, 0.025},
		{ "down_25.11.16",   8185, 8353, 0.025},
		{ "up_28.11.16",     8357, 8430, 0.025},
		{ "mid_01.12.16",    8470, 8571, 0.025},
		{ "down_02.12.16",   8574, 8738, 0.025},
		{ "up_05.12.16",     8741, 8869, 0.025},
		{ "mid_07.12.16",    8873, 9009, 0.025},
		{ "up_12.12.16",     9012, 9112, 0.025},
		{ "mid_14.12.16",    9116, 9245, 0.025},
		{ "down_16.12.16",   9253, 9470, 0.025},
		{ "up_19.12.16",     9475, 9600, 0.025},
		{ "mid_21.12.16",    9603, 9712, 0.025},
		{ "down_23.12.16",   9715, 9869, 0.025},
		{ "up_26.12.16",     9871, 10019, 0.025},
		{ "mid_28.12.16",    10021, 10171, 0.025},
		{ "down_30.12.16",   10175, 10307, 0.025},
		{ "down_02.01.17",   10308, 10356, 0.025},
		{ "down_03.01.17",   10357, 10424, 0.025},
		{ "up_04.01.17",     10433, 10832, 0.025},
		{ "mid_11.01.17",    10834, 10973, 0.025},
		{ "down_13.01.17",   10979, 11147, 0.025},
		{ "up_16.01.17",     11150, 11267, 0.025},
		{ "mid_18.01.17",    11271, 11401, 0.025},
		{ "down_20.01.17",   11404, 11563, 0.025},
		{ "up_23.01.17",     11570, 11619, 0.025}
	};
	int i;
	int N;
	
	N = sizeof(positions) / sizeof(positions[0]);
	for (i=0; i<N; i++) 
		spectr5(positions[i].name, mask, positions[i].first, positions[i].last, positions[i].bgnd);
}
