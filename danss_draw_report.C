void danss_draw_report(void)
{
	const char fname[] = "danss_report.root";
	const char pname[] = "danss_report.pdf";
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
		{ "stuck_10.05.16",  3400, 3788, 0.1},
		{ "up_03.06.16",     4261, 4400, 0.05},
		{ "down_06.06.16",   4403, 4439, 0.05},
		{ "up_08.06.16",     4508, 4601, 0.05},
		{ "down_10.06.16",   4604, 4643, 0.05},
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
	TCanvas *cv;
	TFile *f;
	TH1 *h;
	char str[1024];
	TVirtualPad *pd;
	TLatex *txt;
	double val, err;
	double evnts, rate;
	
	gStyle->SetOptStat(0);
	N = sizeof(positions) / sizeof(positions[0]);
	cv = new TCanvas("CV", "Results", 1600, 800);
	f = new TFile(fname);
	if (!f->IsOpen()) return;
	sprintf(str, "%s[", pname);
	cv->Print(str);
	txt = new TLatex();
	txt->SetTextSize(0.05);
	
	for (i=0; i<N; i++) {
		cv->Clear();
		cv->Divide(2, 1);
		pd = cv->cd(1);
		pd->Divide(1, 2);
		pd->cd(1);

		sprintf(str, "%s_hSig-sig", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) {
			printf("%s not found in %s\n", positions[i].name, fname);
			continue;
		}
		h->SetLineWidth(2);
		h->SetLineColor(kBlue);
		h->SetTitle("Neutrino events after veto");
		h->Draw();
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Signal = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kBlue);
		txt->DrawTextNDC(0.4, 0.8, str);
		sprintf(str, "%s_hSig-rand", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kBlack);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Random = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kBlack);
		txt->DrawTextNDC(0.4, 0.7, str);
		sprintf(str, "%s_hSig-diff", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kRed);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Result = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kRed);
		txt->DrawTextNDC(0.4, 0.6, str);
		evnts = val;
		
		pd->cd(2);
		sprintf(str, "%s_hCosm-sig", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kBlue);
		h->SetTitle("Cosmic background");
		h->Draw();
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Signal = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kBlue);
		txt->DrawTextNDC(0.4, 0.8, str);
		sprintf(str, "%s_hCosm-rand", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kBlack);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Random = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kBlack);
		txt->DrawTextNDC(0.4, 0.7, str);
		sprintf(str, "%s_hCosm-diff", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kRed);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Result = %5.0f+-%4.0f", val, err);
		txt->SetTextColor(kRed);
		txt->DrawTextNDC(0.4, 0.6, str);
		
		cv->cd(2);
		sprintf(str, "%s_hSig", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		rate = h->Integral(1, 35);
		sprintf(str, "%s_hRes", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kRed);
		sprintf(str, "%s (%d - %d)", positions[i].name, positions[i].first, positions[i].last);
		h->SetTitle(str);
		h->Draw();
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Neutrio = %5.2f+-%4.2f mHz", val, err);
		txt->SetTextColor(kRed);
		txt->DrawTextNDC(0.3, 0.85, str);
		sprintf(str, "%s_hCosm", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kBlue);
		h->SetFillStyle(kSolid);
		h->SetFillColor(kBlue);
		h->Scale(positions[i].bgnd);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Cosm = %5.2f+-%4.2f mHz (%4.1f%%)", val, err, 100*positions[i].bgnd);
		txt->SetTextColor(kBlue);
		txt->DrawTextNDC(0.3, 0.8, str);
		sprintf(str, "Time = %6.1f*10^3 s", evnts / rate);
		txt->SetTextColor(kBlack);
		txt->DrawTextNDC(0.3, 0.75, str);
		cv->Update();
//		getchar();
		cv->Print(pname);
	}
	sprintf(str, "%s]", pname);
	cv->Print(str);
	f->Close();
	delete cv;
}
