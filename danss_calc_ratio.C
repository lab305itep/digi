void danss_calc_ratio(void)
{
#ifdef STRONG_CUTS
	const char fname[] = "danss_report_strong.root";
	const char pname[] = "danss_ratio_strong.pdf";
#else
	const char fname[] = "danss_report.root";
	const char pname[] = "danss_ratio.pdf";
#endif
	const struct {
		char name[32];
		int first;
		int last;
	} positions[] = {
		{ "down_21.04.16",   2307, 2361},	// [0] - 0-
		{ "up_22.04.16",     2366, 2387},	// [1] - 0+
		{ "down_23.04.16",   2399, 2445},	// [2] - 0-
		{ "mid_24.04.16",    2449, 2512},
		{ "up_25.04.16",     2514, 2562},	// [4] - 1+
		{ "down_26.04.16",   2564, 2609},	// [5] - 1-
		{ "mid_27.04.16",    2620, 2685},
		{ "up_28.04.16",     2687, 2730},	// [7] - 2+
		{ "down_29.04.16",   2733, 2788},	// [8] - 2-
		{ "mid_30.04.16",    2791, 2832},
		{ "stuck_01.05.16",  2836, 3399},
		{ "stuck_10.05.16",  3400, 3788},
		{ "up_03.06.16",     4261, 4400},	// [12] - 3+
		{ "down_06.06.16",   4403, 4439},	// [13] - 3-
		{ "up_08.06.16",     4508, 4601},	// [14] - 3+
		{ "raised_30.09.16", 5540, 5807},
		{ "raised_04.10.16", 5808, 5903},	// veto corners on
		{ "mid_05.10.16",    5907, 5995},	// [17]
		{ "up_10.10.16",     6007, 6130},	// [18] - 5+
//		{ "mid_12.10.16",    6136, 6275},
//		{ "down_14.10.16",   6278, 6467},
//		{ "up_17.10.16",     6469, 6570},
//		{ "mid_21.10.16",    6573, 6582},
//		{ "down_21.10.16",   6587, 6745},
//		{ "up_24.10.16",     6757, 6815},
//		{ "mid_27.10.16",    6842, 6923},
//		{ "down_28.10.16",   6926, 7095},
		{ "up_31.10.16",     7106, 7364},	// [19] - 5+
		{ "mid_11.11.16",    7387, 7406},	// [20]
		{ "down_11.11.16",   7418, 7458},	// [21] - 5-
		{ "up_14.11.16",     7478, 7579},	// [22] - 6+
//		{ "mid_16.11.16",    7581, 7717},
		{ "down_18.11.16",   7727, 7913},	// [23] - 6-
		{ "up_21.11.16",     7922, 8042},	// [24] - 7+
		{ "mid_23.11.16",    8048, 8179},	// [25]
		{ "down_25.11.16",   8185, 8353},	// [26] - 7-
		{ "up_28.11.16",     8357, 8430},	// [27] - 8+
		{ "mid_01.12.16",    8470, 8571},	// [28]
		{ "down_02.12.16",   8574, 8738},	// [29] - 8-
		{ "up_05.12.16",     8741, 8869},	// [30] - 8+
		{ "mid_07.12.16",    8873, 9009},	// [31]
		{ "up_12.12.16",     9012, 9112},	// [32] - 9+
		{ "mid_14.12.16",    9116, 9245},	// [33]
		{ "down_16.12.16",   9253, 9470},	// [34] - 9-
//		{ "up_19.12.16",     9475, 9600},
//		{ "mid_21.12.16",    9603, 9712},
//		{ "down_23.12.16",   9715, 9869},
//		{ "up_26.12.16",     9871, 10019}
//		{ "mid_28.12.16",    10021, 10171},
//		{ "down_30.12.16",   10175, 10307},
//		{ "down_02.01.17",   10308, 10356},
//		{ "down_03.01.17",   10357, 10424},
		{ "up_04.01.17",     10433, 10832},	// [35] - 10+
		{ "mid_11.01.17",    10834, 10973},	// [36]
		{ "down_13.01.17",   10979, 11147},	// [37] - 10-
		{ "up_16.01.17",     11150, 11267},	// [38] - 11+
		{ "mid_18.01.17",    11271, 11401},	// [39]
		{ "down_20.01.17",   11404, 11563},	// [40] - 11-
		{ "up_23.01.17",     11570, 11619}	// [41] - 11+
	};
	int i, j;
	int N;
	TCanvas *cv;
	TFile *f;
	TH1 *h[sizeof(positions) / sizeof(positions[0])];
	TH1 *hSum[12];
	TH1 *hRatio[7];
	TH1 *hTmp[3];
	char str[1024];
	TLatex *txt;
	double val, err;
	TVirtualPad *pd;
	const char periods[3][30] = {"April-June", "October-November", "December-January"};
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	N = sizeof(positions) / sizeof(positions[0]);
	cv = new TCanvas("CV", "Results", 2400, 1200);
	f = new TFile(fname, "UPDATE");
	if (!f->IsOpen()) return;
	sprintf(str, "%s[", pname);
	cv->Print(str);
	txt = new TLatex();
	
	for (i=0; i<N; i++) {
		sprintf(str, "%s_hRes", positions[i].name);
		h[i] = (TH1 *) f->Get(str);
		if (!h[i]) {
			printf("%s not found in %s\n", positions[i].name, fname);
			return;
		}
		h[i]->SetBit(TH1::kIsAverage);
	}
	
	hSum[0] = (TH1 *) h[0]->Clone("hUp");
	hSum[1] = (TH1 *) h[0]->Clone("hDown");
	hSum[2] = (TH1 *) h[0]->Clone("hMid");
	hSum[3] = (TH1 *) h[0]->Clone("hSum");
	hSum[4] = (TH1 *) h[0]->Clone("hSum_first");
	hSum[5] = (TH1 *) h[0]->Clone("hSum_last");
	hSum[6] = (TH1 *) h[0]->Clone("hUpAprilJune");
	hSum[7] = (TH1 *) h[0]->Clone("hDownAprilJune");
	hSum[8] = (TH1 *) h[0]->Clone("hUpOctoberNovember");
	hSum[9] = (TH1 *) h[0]->Clone("hDownOctoberNovember");
	hSum[10] = (TH1 *) h[0]->Clone("hUpDecemberJanuary");
	hSum[11] = (TH1 *) h[0]->Clone("hDownDecemberJanuary");
	hRatio[0] = (TH1 *) h[0]->Clone("hRatio");
	hRatio[1] = (TH1 *) h[0]->Clone("hMix");
	hRatio[2] = (TH1 *) h[0]->Clone("hRatio2");
	hRatio[3] = (TH1 *) h[0]->Clone("hRatioMidUp");
	hRatio[4] = (TH1 *) h[0]->Clone("hRatioAprilJune");
	hRatio[5] = (TH1 *) h[0]->Clone("hRatioOctoberNovember");
	hRatio[6] = (TH1 *) h[0]->Clone("hRatioDecemberJanuary");
	hTmp[0] = (TH1 *) h[0]->Clone("hTmpUp");
	hTmp[1] = (TH1 *) h[0]->Clone("hTmpDown");
	hTmp[2] = (TH1 *) h[0]->Clone("hTmpRatio");
	for (i=0; i<12; i++) hSum[i]->Reset();
	for (i=0; i<12; i++) hSum[i]->SetLineWidth(2);
	
	for (i=0; i<N; i++) {
		switch (positions[i].name[0]) {
			case 'u' : j = 0; break;
			case 'd' : j = 1; break;
			case 'm' : j = 2; break;
			default : j = -1; break;
		}
		if (j >= 0) hSum[j]->Add(h[i]);
		hSum[(i < N/2) ? 4 : 5]->Add(h[i]);
		hSum[3]->Add(h[i]);
	}
//	for (i=0; i<3; i++) hSum[3]->Add(hSum[i]);
	hRatio[0]->Divide(hSum[1], hSum[0]);
	hRatio[1]->Divide(hSum[5], hSum[4]);
	hRatio[3]->Divide(hSum[2], hSum[0]);
//	Ratio from short ratios
	hRatio[2]->Reset();
	hRatio[2]->SetBit(TH1::kIsAverage);
//	0+2 / 1
	hTmp[0]->Add(h[0], h[2]);
	hTmp[2]->Divide(hTmp[0], h[1]);
	hRatio[2]->Add(hTmp[2]);
//	5 / 4
	hTmp[2]->Divide(h[5], h[4]);
	hRatio[2]->Add(hTmp[2]);
//	8 / 7
	hTmp[2]->Divide(h[8], h[7]);
	hRatio[2]->Add(hTmp[2]);
//	13 / 12 + 14
	hTmp[0]->Add(h[12], h[14]);
	hTmp[2]->Divide(h[13], hTmp[0]);
	hRatio[2]->Add(hTmp[2]);
//	21 / 18+19
	hTmp[0]->Add(h[18], h[19]);
	hTmp[2]->Divide(h[21], hTmp[0]);
	hRatio[2]->Add(hTmp[2]);
//	23 / 22
	hTmp[2]->Divide(h[23], h[22]);
	hRatio[2]->Add(hTmp[2]);
//	26 / 24
	hTmp[2]->Divide(h[26], h[24]);
	hRatio[2]->Add(hTmp[2]);
//	29 / 27 + 30
	hTmp[0]->Add(h[27], h[30]);
	hTmp[2]->Divide(h[29], hTmp[0]);
	hRatio[2]->Add(hTmp[2]);
//	34 / 32
	hTmp[2]->Divide(h[34], h[32]);
	hRatio[2]->Add(hTmp[2]);
//	37 / 35
	hTmp[2]->Divide(h[37], h[35]);
	hRatio[2]->Add(hTmp[2]);
//	40 / 38 + 41
	hTmp[0]->Add(h[38], h[41]);
	hTmp[2]->Divide(h[40], hTmp[0]);
	hRatio[2]->Add(hTmp[2]);
//	AprilJune (0 + 2 + 5 + 8 + 13) / (1 + 4 + 7 + 12 + 14)
	hSum[6]->Add(h[1]);
	hSum[6]->Add(h[4]);
	hSum[6]->Add(h[7]);
	hSum[6]->Add(h[12]);
	hSum[6]->Add(h[14]);
	hSum[7]->Add(h[0]);
	hSum[7]->Add(h[2]);
	hSum[7]->Add(h[5]);
	hSum[7]->Add(h[8]);
	hSum[7]->Add(h[13]);
	hRatio[4]->Divide(hSum[7], hSum[6]);
//	OctoberNovember (21 + 23 + 26 + 29) / (18 + 19 + 22 + 24 + 27 + 30)
	hSum[8]->Add(h[18]);
	hSum[8]->Add(h[19]);
	hSum[8]->Add(h[22]);
	hSum[8]->Add(h[24]);
	hSum[8]->Add(h[27]);
	hSum[8]->Add(h[30]);
	hSum[9]->Add(h[21]);
	hSum[9]->Add(h[23]);
	hSum[9]->Add(h[26]);
	hSum[9]->Add(h[29]);
	hRatio[5]->Divide(hSum[9], hSum[8]);
//	DecemberJanuary (34 + 37 + 40) / (32 + 35 + 38 + 41)
	hSum[10]->Add(h[32]);
	hSum[10]->Add(h[35]);
	hSum[10]->Add(h[38]);
	hSum[10]->Add(h[41]);
	hSum[11]->Add(h[34]);
	hSum[11]->Add(h[37]);
	hSum[11]->Add(h[40]);
	hRatio[6]->Divide(hSum[11], hSum[10]);
	cv->Divide(2, 1);

	pd = cv->cd(1);
	pd->Divide(1, 2);
	pd->cd(1);
	hSum[0]->SetTitle("Up - Middle - Down;MeV;mHz");
	hSum[0]->SetLineColor(kRed);
	hSum[2]->SetLineColor(kGreen);
	hSum[1]->SetLineColor(kBlue);
	for (i=0; i<4; i++) hSum[i]->SetLineWidth(2);
	hSum[0]->Draw();
	hSum[2]->Draw("same");
	hSum[1]->Draw("same");
	val = hSum[0]->IntegralAndError(1, 35, err);
	sprintf(str, "Up = %5.2f+-%4.2f", val, err);
	txt->SetTextColor(kRed);
	txt->DrawTextNDC(0.4, 0.85, str);
	val = hSum[2]->IntegralAndError(1, 35, err);
	sprintf(str, "Mid = %5.2f+-%4.2f", val, err);
	txt->SetTextColor(kGreen);
	txt->DrawTextNDC(0.4, 0.8, str);
	val = hSum[1]->IntegralAndError(1, 35, err);
	sprintf(str, "Down = %5.2f+-%4.2f", val, err);
	txt->SetTextColor(kBlue);
	txt->DrawTextNDC(0.4, 0.75, str);
	
	pd->cd(2);
	hSum[3]->SetTitle("Average of all positions;MeV;mHz");
	hSum[3]->SetLineColor(kBlack);
	hSum[3]->Draw();
	val = hSum[3]->IntegralAndError(1, 35, err);
	sprintf(str, "Avr = %5.2f+-%4.2f", val, err);
	txt->SetTextColor(kBlack);
	txt->DrawTextNDC(0.4, 0.8, str);

	pd = cv->cd(2);
	pd->Divide(1, 4);
	pd->cd(1);
	hRatio[0]->SetTitle("Down/Up all by all;MeV;");
	hRatio[0]->SetLineColor(kBlack);
	hRatio[0]->SetLineWidth(2);
	hRatio[0]->SetMinimum(0.5);
	hRatio[0]->SetMaximum(1.0);
	hRatio[0]->Fit("pol0", "", "", 1, 7);
	
	pd->cd(2);
	hRatio[2]->SetTitle("Down/Up pair by pair;MeV;");
	hRatio[2]->SetLineColor(kBlack);
	hRatio[2]->SetLineWidth(2);
	hRatio[2]->SetMinimum(0.5);
	hRatio[2]->SetMaximum(1.0);
	hRatio[2]->Fit("pol0", "", "", 1, 7);

	pd->cd(3);
	hRatio[3]->SetTitle("Mid/Up all by all;MeV;");
	hRatio[3]->SetLineColor(kBlack);
	hRatio[3]->SetLineWidth(2);
	hRatio[3]->SetMinimum(0.65);
	hRatio[3]->SetMaximum(1.15);
	hRatio[3]->Fit("pol0", "", "", 1, 7);

	pd->cd(4);
	hRatio[1]->SetTitle("First half/Last half;MeV;");
	hRatio[1]->SetLineColor(kBlack);
	hRatio[1]->SetLineWidth(2);
	hRatio[1]->SetMinimum(0.75);
	hRatio[1]->SetMaximum(1.25);
	hRatio[1]->Fit("pol0", "", "", 1, 7);

	cv->Update();
	cv->Print(pname);

	cv->Clear();
	cv->Divide(2, 2);
	for (i=0; i<3; i++) {
		cv->cd(i+1);
		sprintf(str, "Down/Up all by all for %s;MeV;", periods[i]);
		hRatio[4+i]->SetTitle(str);
		hRatio[4+i]->SetLineColor(kBlack);
		hRatio[4+i]->SetLineWidth(2);
		hRatio[4+i]->SetMinimum(0.5);
		hRatio[4+i]->SetMaximum(1.0);
		hRatio[4+i]->Fit("pol0", "", "", 1, 7);
	}
	cv->Update();
	cv->Print(pname);

	sprintf(str, "%s]", pname);
	cv->Print(str);
	for (i=0; i<12; i++) hSum[i]->Write();
	for (i=0; i<7; i++) hRatio[i]->Write();
	f->Close();
	delete cv;
}
