void digi_statone(int num)
{
	char fname[1024];
	double tm;
	char start[64], stop[64];
	
	sprintf(fname, "danss_root/danss_%6.6d.root", num);
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) {
		printf("File %s not found\n", fname);
		return;
	}

	TTree *info = f->Get("DanssInfo");
	if (!info) {
		printf("File %s - no info\n", fname);
		return;
	}
	info->GetEntry(0);
	tm = info->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
//	time_t tStart = info->GetLeaf("startTime")->GetValue();
	time_t tStop = info->GetLeaf("stopTime")->GetValue();
	
	TTree *evt = f->Get("DanssEvent");
	if (!evt) {
		printf("File %s - no events\n", fname);
		return;
	}

	evt->GetEntry(0);
	time_t tStart = evt->GetLeaf("unixTime")->GetValue();

	int N = evt->GetEntries();
	
	int veto = evt->GetEntries("VetoCleanHits > 1 || VetoCleanEnergy > 4 || PmtCleanEnergy + SiPmCleanEnergy > 40");
	int vetoOnly = evt->GetEntries("(VetoCleanHits > 1 || VetoCleanEnergy > 4) && !(PmtCleanEnergy + SiPmCleanEnergy > 40)");
	int danssOnly = evt->GetEntries("!(VetoCleanHits > 1 || VetoCleanEnergy > 4) && (PmtCleanEnergy + SiPmCleanEnergy > 40)");

	int gt1MeV = evt->GetEntries("PmtCleanEnergy + SiPmCleanEnergy > 2");
	int gt3MeV = evt->GetEntries("PmtCleanEnergy + SiPmCleanEnergy > 6");
	int gt20MeV = evt->GetEntries("PmtCleanEnergy + SiPmCleanEnergy > 40");
	int positrons = evt->GetEntries("PositronPmtEnergy + PositronSiPmEnergy > 2 && PositronPmtEnergy + PositronSiPmEnergy < 16");
	int neutrons = evt->GetEntries("PmtCleanEnergy + SiPmCleanEnergy > 6 && PmtCleanEnergy + SiPmCleanEnergy < 20 && SiPmCleanHits > 2");
//	
	strftime(start, sizeof(start), "%F %R", localtime(&tStart));
	strftime(stop , sizeof(stop) , "%R", localtime(&tStop));
	printf("%4d  %s %s %5.0f   %8d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", 
		num, start, stop, tm, N, N/tm, veto/tm, vetoOnly/tm, danssOnly/tm, gt1MeV/tm, gt3MeV/tm, gt20MeV/tm, positrons/tm, neutrons/tm);

	f->Close();
}

void digi_stattitle(void)
{
	printf("Run    Start           Stop   len, s   Events   Trig   Veto  NoPMT NoVeto  >1MeV  >3MeV >20MeV     e+     n\n");
}

void digi_statlist(void)
{
	const int runs[] = {2320, 2400, 2600, 2750, 2375, 2550, 2700, 2500, 2650, 2800, 3001, 3500};
	int i;
	
	digi_stattitle();
	for (i=0; i<sizeof(runs)/sizeof(runs[0]); i++) digi_statone(runs[i]);
}

void digi_stat(int first, int last)
{
	int i;
	
	digi_stattitle();
	for (i=first; i<=last; i++) digi_statone(i);
}

