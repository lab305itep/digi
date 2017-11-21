void OnOffProject(TH1D* hist, const char *what, TCut cut)
{
	double tSig;
	double tBgnd;
	int i;
//	const int range[] = {18684, 23041, 23082, 24706};
	const int range[] = {18684, 18690, 23082, 23090};

	TChain *chSig  = new TChain("DanssEvent", "Signal");
	TChain *chBgnd = new TChain("DanssEvent", "Background");
	TChain *chSigInfo  = new TChain("DanssInfo", "SignalInfo");
	TChain *chBgndInfo = new TChain("DanssInfo", "BackgrondInfo");
	
	Add2Chain(chSig,  range[0], range[1], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chSigInfo,  range[0], range[1], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chBgnd, range[2], range[3], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chBgndInfo, range[2], range[3], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	
	TH1D *hSig = (TH1D*) hist->Clone("_hSig");
	TH1D *hBgnd = (TH1D*) hist->Clone("_hBgnd");
	chSig->Project(hSig->GetName(), what, cut);
	chBgnd->Project(hBgnd->GetName(), what, cut);
	hSig->Sumw2();
	hBgnd->Sumw2();
	
	tSig = 0;
	for (i=0; i<chSigInfo->GetEntries(); i++) {
		chSigInfo->GetEntry(i);
		tSig += chSigInfo->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
	}
	tSig /= 86400.0;	// s/day

	tBgnd = 0;
	for (i=0; i<chBgndInfo->GetEntries(); i++) {
		chBgndInfo->GetEntry(i);
		tBgnd += chBgndInfo->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
	}
	
	tBgnd /= 86400.0;	// s/day
	
	printf("Signal: %5.2f days / %8.0f events    Background: %5.2f days / %8.0f events\n",
		tSig, hSig->Integral(), tBgnd, hBgnd->Integral());
	hSig->Scale(1.0/tSig);
	hBgnd->Scale(1.0/tBgnd);
	hist->Add(hSig, hBgnd, 1.0, -1.0);
	
	hSig->SetLineColor(kRed);
	hBgnd->SetLineColor(kBlue);
	hist->SetLineColor(kGreen);
	
	hSig->DrawCopy();
	hBgnd->DrawCopy("same");
	hist->DrawCopy("same");
	
	delete hSig;
	delete hBgnd;
	delete chSig;
	delete chBgnd;
	delete chSigInfo;
	delete chBgndInfo;
}

void look4Fe(void)
{
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	
	TH1D *h = new TH1D("Diff", "Clean energy, Reactor on - off", 100, 0, 10);
	OnOffProject(h, "(SiPmCleanEnergy + PmtCleanEnergy)/2", cX && cY && cZ);
//	h->Draw();
}
