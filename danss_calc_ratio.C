void danss_calc_ratio(void)
{
	const char fname[] = "danss_report.root";
	const char pname[] = "danss_ratio.pdf";
	const struct {
		char name[32];
		int first;
		int last;
	} positions[] = {
//		{ "raised_30.09.16", 5540, 5807},
//		{ "raised_04.10.16", 5808, 5903},	// veto corners on
		{ "mid_05.10.16",    5907, 5995},	// [0]
		{ "up_10.10.16",     6007, 6130},	// [1] - 0+
//		{ "mid_12.10.16",    6136, 6275},
//		{ "down_14.10.16",   6278, 6467},
//		{ "up_17.10.16",     6469, 6570},
//		{ "mid_21.10.16",    6573, 6582},
//		{ "down_21.10.16",   6587, 6745},
//		{ "up_24.10.16",     6757, 6815},
//		{ "mid_27.10.16",    6842, 6923},
//		{ "down_28.10.16",   6926, 7095},
		{ "up_31.10.16",     7106, 7364},	// [2] - 0+
		{ "mid_11.11.16",    7387, 7406},	// [3]
		{ "down_11.11.16",   7418, 7458},	// [4] - 0-
		{ "up_14.11.16",     7478, 7579},	// [5] - 1+
//		{ "mid_16.11.16",    7581, 7717},
		{ "down_18.11.16",   7727, 7913},	// [6] - 1-
		{ "up_21.11.16",     7922, 8042},	// [7] - 2+
		{ "mid_23.11.16",    8048, 8179},	// [8]
		{ "down_25.11.16",   8185, 8353},	// [9] - 2-
		{ "up_28.11.16",     8357, 8430},	// [10] - 3+
		{ "mid_01.12.16",    8470, 8571},	// [11]
		{ "down_02.12.16",   8574, 8738},	// [12] - 3-
		{ "up_05.12.16",     8741, 8869},	// [13] - 4+
		{ "mid_07.12.16",    8873, 9009},	// [14]
		{ "up_12.12.16",     9012, 9112},	// [15] - 4+
		{ "mid_14.12.16",    9116, 9245},	// [16]
		{ "down_16.12.16",   9253, 9470},	// [17] - 4-
//		{ "up_19.12.16",     9475, 9600},
//		{ "mid_21.12.16",    9603, 9712},
//		{ "down_23.12.16",   9715, 9869},
//		{ "up_26.12.16",     9871, 10019}
//		{ "mid_28.12.16",    10021, 10171},
//		{ "down_30.12.16",   10175, 10307},
//		{ "down_02.01.17",   10308, 10356},
//		{ "down_03.01.17",   10357, 10424},
		{ "up_04.01.17",     10433, 10832},	// [18] - 5+
		{ "mid_11.01.17",    10834, 10973},	// [19]
		{ "down_13.01.17",   10979, 11147},	// [20] - 5-
		{ "up_16.01.17",     11150, 11267},	// [21] - 6+
		{ "mid_18.01.17",    11271, 11401},	// [22]
		{ "down_20.01.17",   11404, 11563},	// [23] - 6-
		{ "up_23.01.17",     11570, 11619}	// [24] - 6+
	};
	int i, j;
	int N;
	TCanvas *cv;
	TFile *f;
	TH1 *h[sizeof(positions) / sizeof(positions[0])];
	TH1 *hSum[6];
	TH1 *hRatio[3];
	TH1 *hTmp[3];
	char str[1024];
	TLatex *txt;
	double val, err;
	TVirtualPad *pd;
	
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
	hRatio[0] = (TH1 *) h[0]->Clone("hRatio");
	hRatio[1] = (TH1 *) h[0]->Clone("hMix");
	hRatio[2] = (TH1 *) h[0]->Clone("hRatio2");
	hTmp[0] = (TH1 *) h[0]->Clone("hTmpUp");
	hTmp[1] = (TH1 *) h[0]->Clone("hTmpDown");
	hTmp[2] = (TH1 *) h[0]->Clone("hTmpRatio");
	for (i=0; i<4; i++) hSum[i]->Reset();
	for (i=0; i<4; i++) hSum[i]->SetLineWidth(2);
	
	for (i=0; i<N; i++) {
		switch (positions[i].name[0]) {
			case 'u' : j = 0; break;
			case 'd' : j = 1; break;
			case 'm' : j = 2; break;
			default : j = -1; break;
		}
		if (j >= 0) hSum[j]->Add(h[i]);
		hSum[(i < N/2) ? 4 : 5]->Add(h[i]);
	}
	for (i=0; i<3; i++) hSum[3]->Add(hSum[i]);
	hRatio[0]->Divide(hSum[1], hSum[0]);
	hRatio[1]->Divide(hSum[5], hSum[4]);
	
	hRatio2->Reset();
	hRatio2->SetBit(TH1::kIsAverage);

	hTmp[0]->Add(h[1], h[2]);
//	hTmp[1]->Add(h[4], h[4], 1, 0);
	hTmp[2]->Divide(h[4], hTmp[0]);
	hRatio2->Add(hTmp[2]);
	
//	hTmp[0]->Add(h[5], h[5], 1, 0);
//	hTmp[1]->Add(h[6], h[6], 1, 0);
	hTmp[2]->Divide(h[6], h[5]);
	hRatio2->Add(hTmp[2]);
	
//	hTmp[0]->Add(h[7], h[7], 1, 0);
//	hTmp[1]->Add(h[9], h[9], 1, 0);
	hTmp[2]->Divide(h[9], h[7]);
	hRatio2->Add(hTmp[2]);
	
//	hTmp[0]->Add(h[10], h[10], 1, 0);
//	hTmp[1]->Add(h[12], h[12], 1, 0);
	hTmp[2]->Divide(h[12], h[10]);
	hRatio2->Add(hTmp[2]);
	
	hTmp[0]->Add(h[13], h[15]);
//	hTmp[1]->Add(h[17], h[17], 1, 0);
	hTmp[2]->Divide(h[17], hTmp[0]);
	hRatio2->Add(hTmp[2]);
	
//	hTmp[0]->Add(h[18], h[18], 1, 0);
//	hTmp[1]->Add(h[20], h[20], 1, 0);
	hTmp[2]->Divide(h[20], h[18]);
	hRatio2->Add(hTmp[2]);
	
	hTmp[0]->Add(h[21], h[24]);
//	hTmp[1]->Add(h[23], h[23], 1, 0);
	hTmp[2]->Divide(h[23], hTmp[0]);
	hRatio2->Add(hTmp[2]);
	
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
	pd->Divide(1, 3);
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
	hRatio[1]->SetTitle("First half/Last half;MeV;");
	hRatio[1]->SetLineColor(kBlack);
	hRatio[1]->SetLineWidth(2);
	hRatio[1]->SetMinimum(0.75);
	hRatio[1]->SetMaximum(1.25);
	hRatio[1]->Fit("pol0", "", "", 1, 7);

	cv->Update();
	cv->Print(pname);

	sprintf(str, "%s]", pname);
	cv->Print(str);
	for (i=0; i<6; i++) hSum[i]->Write();
	for (i=0; i<3; i++) hRatio[i]->Write();
	f->Close();
	delete cv;
}
