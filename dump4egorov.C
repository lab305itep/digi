void dump4egorov(const char *fname)
{
	int i;
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) return;
	
	TH1 *hUp = f->Get("hUp_28");
	TH1 *hMid = f->Get("hMid_28");
	TH1 *hDown = f->Get("hDown_28");
	
	if (!hUp || !hMid || !hDown) {
		printf("Something is wrong.\n");
		return;
	}
	
	printf(" Energy         Up            Middle        Down\n");
//              123456789  123456789012  123456789012  123456789012  
	for (i=0; i<28; i++) printf("%4.2f-%4.2f  %6.2f+-%4.2f  %6.2f+-%4.2f  %6.2f+-%4.2f\n",
		0.25*i + 1, 0.25*i + 1.25,
		hUp->GetBinContent(i+1), hUp->GetBinError(i+1),
		hMid->GetBinContent(i+1), hMid->GetBinError(i+1),
		hDown->GetBinContent(i+1), hDown->GetBinError(i+1));
	
	f->Close();
}
