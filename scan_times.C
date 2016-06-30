#include "evtbuilder.h"

void scan_times(char *flist)
{
	FILE *fl;
	char str[2048];
	TFile *f;
	TTree *ti;
	TTree *te;
	int S[3];
	int t[2];
	float t125[2];
	float gt;	

	fl = fopen(flist, "rt");
	if (!fl) {
		printf("File not found %s\n", flist);
		return;
	}

	printf("Run                          S1      S2      S3\n");
	//      danss_root/danss_123456.root 123456781234567812345678
	memset(S, 0, sizeof(S));
	for(;!feof(fl);) {
		if (!fgets(str, sizeof(str), fl)) break;
		if (!strlen(str)) continue;
		str[strlen(str) - 1] = '\0';
		f = new TFile(str);
		if (!f->IsOpen()) {
			printf("File not found %s\n", str);
			continue;
		}
		ti = (TTree *) f->Get("DanssInfo");
		te = (TTree *) f->Get("DanssEvent");
		if (!ti || !te) {
			printf("Wrong file %s\n", str);
			delete f;
			continue;
		}

		ti->GetEntry(0);
		gt = ti->GetLeaf("gTime")->GetValueLong64() / GLOBALFREQ;

		te->GetEntry(0);
		t[0] = te->GetLeaf("unixTime")->GetValue();
		t125[0] = te->GetLeaf("globalTime")->GetValueLong64() / GLOBALFREQ;
		te->GetEntry(te->GetEntries() - 1);
		t[1] = te->GetLeaf("unixTime")->GetValue();
		t125[1] = te->GetLeaf("globalTime")->GetValueLong64() / GLOBALFREQ;

		printf("%s %7.0f %7.0f %7d\n", str, gt, t125[1] - t125[0], t[1] - t[0]);
		S[0] += gt;
		S[1] += t125[1] - t125[0];
		S[2] += t[1] - t[0];
		delete f;
	}

	printf("TOTAL sums: %d, %d, %d\n", S[0], S[1], S[2]);

	fclose(fl);
}

