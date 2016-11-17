int UnitNumber;
int ChannelNumber;
TFile *RootFile;
TCanvas *cv;

void Open(const char *name)
{
	RootFile = new TFile(name);
	UnitNumber = 1;
	ChannelNumber = 0;
}

void Draw(int update = 1)
{
	char str[128];
	int i;
	TH1 *h;
	
	if (!RootFile->IsOpen()) return;
	cv = (TCanvas *) gROOT->FindObject("CVT");
	if (!cv) cv = new TCanvas("CVT", "Canvas", 1500, 1200);
	cv->Clear();
	cv->Divide(4, 4);
	for (i=0; i<16; i++) {
		cv->cd(i+1);
		sprintf(str, "hDT%2.2dc%2.2d", UnitNumber, ChannelNumber + i);
		h = (TH1 *) RootFile->Get(str);
		if (h) h->Draw();
	}
	if (update) cv->Update();
}

void Next(void)
{
	Draw();
	ChannelNumber += 16;
	if (ChannelNumber >= 64) {
		ChannelNumber  = 0;
		UnitNumber++;
		if (UnitNumber >= 50) UnitNumber = 0;
	}
}

void Prev(void)
{
	Draw();
	ChannelNumber -= 16;
	if (ChannelNumber < 0) {
		ChannelNumber  = 48;
		UnitNumber--;
		if (UnitNumber < 1) UnitNumber = 47;
	}
}

void Print(const char *name)
{
	char str[1024];
	char *ptr;
	strncpy(str, name, sizeof(str)-10);
	ptr = strrchr(str, '.');
	if (!ptr) ptr = strchr(str, '\0');
	strcpy(ptr, ".pdf");
	ptr += strlen(".pdf");
	
	cv = new TCanvas("CVT", "CVT", 1500, 1200);
	strcpy(ptr, "[");
	cv->Print(str);
	*ptr = '\0';
	for (UnitNumber = 1; UnitNumber < 48; UnitNumber++) for (ChannelNumber = 0; ChannelNumber<64; ChannelNumber += 16) {
		Draw(1);
		cv->Print(str);
	}
	strcpy(ptr, "]");
	cv->Print(str);
}

void drawtimehists(const char *name)
{
	Open(name);
	Print(name);
}
