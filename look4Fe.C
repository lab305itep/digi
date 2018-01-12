void OnOffProject(TH1D* hist, TH1D **histS, TH1D **histB, const char *what, TCut cut, int MaxFiles = 10)
{
	double tSig;
	double tBgnd;
	int i;
	const int range[4] = {18684, 23041, 23082, 24706};
	int ActiveRange[2];
	char str[128];
	TObject *tmp;
//	const int range[] = {18684, 18690, 23082, 23090};

	for (i=0; i<2; i++) ActiveRange[i] = (range[2*i + 1] - range[2*i] + 1 > MaxFiles) ? 
		range[2*i] + MaxFiles - 1 : range[2*i + 1];

	TChain *chSig  = new TChain("DanssEvent", "Signal");
	TChain *chBgnd = new TChain("DanssEvent", "Background");
	TChain *chSigInfo  = new TChain("DanssInfo", "SignalInfo");
	TChain *chBgndInfo = new TChain("DanssInfo", "BackgrondInfo");
	
	Add2Chain(chSig,  range[0], ActiveRange[0], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chSigInfo,  range[0], ActiveRange[0], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chBgnd, range[2], ActiveRange[1], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	Add2Chain(chBgndInfo, range[2], ActiveRange[1], 0xE, "/mnt/root0/danss_root4/danss_%6.6d.root");
	
	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetTitleSize(0.05);
	
	sprintf(str, "%s_Sig", hist->GetName());
	tmp = gROOT->FindObject(str);
	if (tmp) delete tmp;
	TH1D *hSig = (TH1D*) hist->Clone(str);
	sprintf(str, "%s_Bgnd", hist->GetName());
	tmp = gROOT->FindObject(str);
	if (tmp) delete tmp;
	TH1D *hBgnd = (TH1D*) hist->Clone(str);
	chSig->Project(hSig->GetName(), what, cut);
	chBgnd->Project(hBgnd->GetName(), what, cut);
	hSig->Sumw2();
	hBgnd->Sumw2();
	
	tSig = 0;
	for (i = 0; i < chSigInfo->GetEntries() && i < MaxFiles; i++) {
		chSigInfo->GetEntry(i);
		tSig += chSigInfo->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
	}
	tSig /= 86400.0;	// s/day

	tBgnd = 0;
	for (i = 0; i < chBgndInfo->GetEntries() && i < MaxFiles; i++) {
		chBgndInfo->GetEntry(i);
		tBgnd += chBgndInfo->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
	}
	
	tBgnd /= 86400.0;	// s/day
	
	printf("Signal: %d files / %5.2f days / %8.0f events    Background: %d files / %5.2f days / %8.0f events\n",
		(int) chSigInfo->GetEntries(), tSig, hSig->Integral(), 
		(int) chBgndInfo->GetEntries(), tBgnd, hBgnd->Integral());
	hSig->Scale(1.0/tSig);
	hBgnd->Scale(1.0/tBgnd);
	hist->Add(hSig, hBgnd, 1.0, -1.0);
	
	hSig->SetLineColor(kRed);
	hBgnd->SetLineColor(kBlue);
	hist->SetLineColor(kGreen);
	
	*histS = hSig;
	*histB = hBgnd;
	
	delete chSig;
	delete chBgnd;
	delete chSigInfo;
	delete chBgndInfo;
}

void look4Fe(int MaxFiles = 10)
{
	TH1D *hClust[3];
	TH1D *hSum[3];
	TVirtualPad *pd;
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	hClust[2] = new TH1D("hClust", "Positron energy, Reactor on - off;MeV;Events/day/100 keV", 90, 1, 10);
	hSum[2] = new TH1D("hSum", "Total energy, Reactor on - off;MeV;Events/day/100 keV", 90, 1, 10);
	OnOffProject(hClust[2], &hClust[0], &hClust[1], "PositronEnergy", cX && cY && cZ, MaxFiles);
	OnOffProject(hSum[2], &hSum[0], &hSum[1], "(SiPmCleanEnergy + PmtCleanEnergy)/2", cX && cY && cZ, MaxFiles);
	
	TLegend *lg = new TLegend(0.6, 0.65, 0.85, 0.8);
	lg->AddEntry(hClust[0], "Reactor ON", "LE");
	lg->AddEntry(hClust[1], "Reactor OFF", "LE");
	lg->AddEntry(hClust[2], "Difference", "LE");
	
	TCanvas *cv = new TCanvas("CV", "CV", 800, 1200);
	cv->SaveAs("look4Fe.pdf[");
	
	cv->Divide(1, 2);
	pd = cv->cd(1);
	pd->SetLogy();
	hClust[0]->Draw();
	hClust[1]->Draw("same");
	hClust[2]->Draw("same");
	lg->Draw();
	pd = cv->cd(2);
	pd->SetLogy();
	hSum[0]->Draw();
	hSum[1]->Draw("same");
	hSum[2]->Draw("same");
	lg->Draw();
	cv->SaveAs("look4Fe.pdf");
	
	cv->Clear();
	cv->Divide(1, 2);
	pd = cv->cd(1);
	hClust[0]->SetTitle("Positron energy, Reactor ON");
	hClust[0]->Draw();
	pd = cv->cd(2);
	hSum[0]->SetTitle("Total energy, Reactor ON");
	hSum[0]->Draw();
	cv->SaveAs("look4Fe.pdf");

	cv->Clear();
	cv->Divide(1, 2);
	pd = cv->cd(1);
	hClust[2]->SetTitle("Positron energy, ON-OFF");
	hClust[2]->Draw();
	pd = cv->cd(2);
	hSum[2]->SetTitle("Total energy, ON-OFF");
	hSum[2]->Draw();
	cv->SaveAs("look4Fe.pdf");
	
	cv->SaveAs("look4Fe.pdf]");
}

void DrawOnOff(const char *fname)
{
	TH1D *h[6];
	TVirtualPad *pd;
	const char hname[6][20] = {"hSum", "hSum_Sig", "hSum_Bgnd", "hClust", "hClust_Sig", "hClust_Bgnd"};
	int i;
	
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) {
		printf("Can not open %s.\n", fname);
		return;
	}
	
	for(i=0; i<6; i++) {
		h[i] = (TH1D*) f->Get(hname[i]);
		if (!h[i]) {
			printf("%s not in flie %s\n", hname[i], fname);
			return;
		}
	}
	
	h[0]->Add(h[1], h[2], 1, -1);
	h[3]->Add(h[4], h[5], 1, -1);
	
	TLegend *lg = new TLegend(0.6, 0.65, 0.85, 0.8);
	lg->AddEntry(h[1], "Reactor ON", "LE");
	lg->AddEntry(h[2], "Reactor OFF", "LE");
	lg->AddEntry(h[0], "Difference", "LE");

	TCanvas *cv = new TCanvas("CV", "CV", 800, 1200);
	cv->SaveAs("look4Fe.pdf[");
	
	cv->Divide(1, 2);
	pd = cv->cd(1);
	pd->SetLogy();
	h[4]->Draw();
	h[3]->Draw("same");
	h[5]->Draw("same");
	lg->Draw();
	pd = cv->cd(2);
	pd->SetLogy();
	h[1]->Draw();
	h[0]->Draw("same");
	h[2]->Draw("same");
	lg->Draw();
	cv->SaveAs("look4Fe.pdf");
	
	cv->Clear();
	cv->Divide(1, 2);
	pd = cv->cd(1);
	h[4]->SetTitle("Positron energy, Reactor ON");
	h[4]->Draw();
	pd = cv->cd(2);
	h[1]->SetTitle("Total energy, Reactor ON");
	h[1]->Draw();
	cv->SaveAs("look4Fe.pdf");

	cv->Clear();
	cv->Divide(1, 2);
	pd = cv->cd(1);
	h[3]->SetTitle("Positron energy, ON-OFF");
	h[3]->Draw();
	pd = cv->cd(2);
	h[0]->SetTitle("Total energy, ON-OFF");
	h[0]->Draw();
	cv->SaveAs("look4Fe.pdf");
	
	cv->SaveAs("look4Fe.pdf]");
	
}
