void digi_statone(int num, const char *root_dir)
{
	char fname[1024];
	double tm;
	char start[64], stop[64];
	
	sprintf(fname, "%s/danss_%6.6d.root", root_dir, num);
	TFile *f = new TFile(fname);
	if (!f->IsOpen()) {
		printf("%d 0 file %s not found\n", num, fname);
		return;
	}

	TTree *info = (TTree *) f->Get("DanssInfo");
	if (!info) {
		printf("%d 0 file %s - no info\n", num, fname);
		return;
	}
	info->GetEntry(0);
	tm = info->GetLeaf("gTime")->GetValueLong64() / 125.0E6;
//	time_t tStart = info->GetLeaf("startTime")->GetValue();
	time_t tStop = info->GetLeaf("stopTime")->GetValue();
	int ipos = info->GetLeaf("position")->GetValue();
	
	TTree *evt = (TTree *) f->Get("DanssEvent");
	if (!evt) {
		printf("%d 0 file %s - no events\n", num, fname);
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
	int positrons = evt->GetEntries("PositronEnergy > 1 && PositronEnergy < 8");
	int neutrons = evt->GetEntries("PmtCleanEnergy + SiPmCleanEnergy > 6 && PmtCleanEnergy + SiPmCleanEnergy < 20 && SiPmCleanHits > 2");
//	
	strftime(start, sizeof(start), "%F %R", localtime(&tStart));
	strftime(stop , sizeof(stop) , "%R", localtime(&tStop));
	printf("%6d %2d %s %s %5.0f   %8d %6.1f %6.1f %6.1f %6.2f %6.1f %6.1f %6.1f %6.1f %6.1f\n", 
		num, ipos, start, stop, tm, N, N/tm, veto/tm, vetoOnly/tm, danssOnly/tm, gt1MeV/tm, gt3MeV/tm, gt20MeV/tm, positrons/tm, neutrons/tm);

	f->Close();
}

void digi_stattitle(void)
{
	printf("Run   Pos Start           Stop   len, s   Events   Trig   Veto  NoPMT NoVeto  >1MeV  >3MeV >20MeV     e+     n\n");
}

void digi_statlist(void)
{
	const int runs[] = {2320, 2400, 2600, 2750, 2375, 2550, 2700, 2500, 2650, 2800, 3001, 3500};
	int i;
	
	digi_stattitle();
	for (i=0; i<sizeof(runs)/sizeof(runs[0]); i++) digi_statone(runs[i], "/mnt/root1/danss_root5");
}

void digi_stat(int first, int last, const char *root_dir = "/mnt/root1/danss_root5")
{
	int i;
	
	digi_stattitle();
	for (i=first; i<=last; i++) digi_statone(i, root_dir);
}
