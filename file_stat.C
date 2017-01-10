#include <time.h>

void file_stat(int num)
{
	char str[1024];
	TFile *f;
	struct {
		long long	upTime;			// running time in terms of 125 MHz
		int		runNumber;		// the run number
		int		startTime;		// linux start time, seconds
		int		stopTime;		// linux stop time, seconds
		int		events;			// number of events
	} info;
	TTree *t;
	long tm;
	
	sprintf(str, "danss_root3/danss_%6.6d.root", num);
	f = new TFile(str);
	if (!f->IsOpen()) {
		printf("File %s not found.\n", str);
		return;
	}
	t = (TTree *) f->Get("DanssInfo");
	if (!t) {
		printf("DanssInfo not found in %s.\n", str);
		return;
	}
	t->SetBranchAddress("Info", &info);
	t->GetEntry(0);
	printf("File:\t%s\n", str);
	printf("Run:\t%d\n", info.runNumber);
	tm = info.startTime;
	printf("Begin:\t%s", ctime(&tm));
	tm = info.stopTime;
	printf("End:\t%s", ctime(&tm));
	printf("Duration:\t%8.1f\n", info.upTime / 125.0E6);
	printf("Events:\t%d\n", info.events);
}
