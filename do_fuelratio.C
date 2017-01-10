void do_fuelratio(void)
{
	TH1D *hApr;
	TH1D *hMay;
	TH1D *hAprMay;
	TH1D *hOct;
	TH1D *rAprMay;
	TH1D *rOctAprMay;
	
	rAprMay = new TH1D("rAprMay", "April to May spectrum ratio;MeV", 35, 1, 8); 
	rOctAprMay = new TH1D("rOctAprMay", "(April + May) to October spectrum ratio;MeV", 35, 1, 8); 
	
	hApr = spectr3("root/2306_2834.root", "root/2306_2834_r.root", "PositronSpectrumApril", 1461229103, 1462074649, 0.05);
	hMay = spectr3("root/2834_3778.root", "root/2834_3778_r.root", "PositronSpectrumMay", 1462074649, 1463848610, 0.05);
	hOct = spectr3("root/5469_6694.root", "root/5469_6694_r.root", "PositronSpectrumOctober", 1473359327, 1476300960, 0.025);
	
	hAprMay = (TH1D *) hApr->Clone("hAprMay");
	hMay->SetBit(TH1::kIsAverage);
	hAprMay->SetBit(TH1::kIsAverage);
	hAprMay->Add(hMay);
	
	rAprMay->Divide(hApr, hMay);
	rOctAprMay->Divide(hAprMay, hOct);
	
	TFile *fOut = new TFile("apr_oct_may.root", "recreate");
	hApr->Write();
	hMay->Write();
	hOct->Write();
	hAprMay->Write();
	rAprMay->Write();
	rOctAprMay->Write();
	fOut->Close();
}
