void danss_draw_report(void)
{
	const char fname[] = "danss_report.root";
	const char pname[] = "danss_report.pdf";
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
	TCanvas *cv;
	TFile *f;
	TH1 *h;
	char str[1024];
	TVirtualPad *pd;
	TLatex *txt;
	double val, err;
	
	gStyle->SetOptStat(0);
	N = sizeof(positions) / sizeof(positions[0]);
	cv = new TCanvas("CV", "Results", 1600, 800);
	f = new TFile(fname);
	if (!f->IsOpen()) return;
	sprintf(str, "%s[", pname);
	cv->Print(str);
	txt = new TLatex();
	
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
		sprintf(str, "%s_hRes", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kRed);
		h->SetTitle(positions[i].name);
		h->Draw();
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Neutrio = %5.1f+-%4.1f mHz", val, err);
		txt->SetTextColor(kRed);
		txt->DrawTextNDC(0.3, 0.85, str);
		sprintf(str, "%s_hCosm", positions[i].name);
		h = (TH1 *) f->Get(str);
		if (!h) continue;
		h->SetLineWidth(2);
		h->SetLineColor(kBlue);
		h->SetFillStyle(kSolid);
		h->SetFillColor(kBlue);
		h->Scale((i) ? 0.025 : 0.05);
		h->Draw("same");
		val = h->IntegralAndError(1, 35, err);
		sprintf(str, "Cosmic = %5.1f+-%4.1f mHz", val, err);
		txt->SetTextColor(kBlue);
		txt->DrawTextNDC(0.3, 0.8, str);
		cv->Update();
//		getchar();
		cv->Print(pname);
	}
	sprintf(str, "%s]", pname);
	cv->Print(str);
	f->Close();
	delete cv;
}
