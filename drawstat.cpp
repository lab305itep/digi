#include <string.h>
#include <stdlib.h>


void drawstat(char * fname = "stat_all.txt", int startfile=5540, int col = 8, int mask = 0xFFFFFFFF, char * listpref = "danss_root/") {
    FILE *f, *fl;
    char line[1024];
    char * ptr;
    int i, first = -1, last, n, type;
    double val, min=1e10, max=-1e10;
    const int mycol[] = {
	(int)kBlack,
	(int)kWhite,
	(int)kRed,
	(int)kGreen,
	(int)kBlue,
	(int)kYellow,
	(int)kMagenta,
	(int)kCyan,
	(int)kGray+2,
	(int)kRed+2,
	(int)kGreen+2,
	(int)kBlue+2,
	(int)kYellow+2,
	(int)kMagenta+2,
	(int)kCyan+2,
	(int)kGray,
	(int)kRed-6,
	(int)kGreen-6,
	(int)kBlue-6,
	(int)kYellow-6,
	(int)kMagenta-6,
	(int)kCyan-6,
	(int)kGray-6
    };

    const char tit[][50] = {
	"Run duration, s",
	"Events in the run",
	"Trigger frequency, Hz",
	"Veto frequency, Hz",
	"No PMT vetos, Hz",
	"No veto >20 MeV, Hz",
	">1 MeV, Hz",
	">3 MeV, Hz",
	">20 MeV, Hz",
	"Positron, Hz",
	"Neutron, Hz"
    };
    TLegend * l;
    const char leg[][20] = {
	"Bad",
	"",
	"Undefined",
	"Down",
	"Middle",
	"Up",
	"Stuck",
	"Old Undef power",
	"Old OFF",
	"Old 30% power",
	"Old 60% power",
	"Old ON",
	"Old Cm",
	"Source Cm",
	"Source Na",
	"Source Co",
	"Test LED",
	"Raised 23 cm",
	"Source Cs"
    };
    
    printf("Usage: drawstat(<statfile>=\"stat_all.txt\", <startfile>=5540, <column>=8(>20Mev), <mask>=0xFFFFFFFF, <listpref>=\"danss_root/\")\n");

    if (col < 0 || col >= sizeof(tit)/sizeof(tit[0])) return;
    
    TGraph * gr[sizeof(leg)/sizeof(leg[0])];
    for (i=0; i<sizeof(gr)/sizeof(gr[0]); i++) {
	gr[i] = new TGraph();
	gr[i]->SetMarkerStyle(20);
	gr[i]->SetMarkerColor(mycol[i]);
    }

    f = fopen(fname, "rt");
    if (!f) return;
    fl = fopen("selected.list", "wt");
    if (!fl) return;
    
    for (;;) {
	fgets(line, sizeof(line), f);
	if (feof(f)) break;
	if (strlen(line) < 2 || line[0] == '*') continue;
	// run #
	ptr = strtok(line, " \t");
	if (!ptr) continue;
	n = strtol(ptr, NULL, 10);
	if (n < startfile) continue;
	// type
	ptr = strtok(NULL, " \t");
	type = strtol(ptr, NULL, 10);
	if (type == 0) continue;	// no file
	if (type < 0) type = -1;
	if (!( ((type == -1) && (mask & 1)) || ((type > 0) && ((1 << type) & mask)) )) continue;
	if (first < 0) first = n;
	// skip times
	ptr = strtok(NULL, " \t");
	ptr = strtok(NULL, " \t");
	ptr = strtok(NULL, " \t");
	// skip unwanted
	for (i=0; i<col; i++) ptr = strtok(NULL, " \t");
	// need this column
	ptr = strtok(NULL, " \t");
	if (ptr) val = strtod(ptr, NULL); else val=0;
	if (val < 0) val = 0;
	if (val > max) max = val;
	if (val < min) min = val;
	// fill graph
	gr[type+1]->SetPoint(gr[type+1]->GetN(), (double)n, val);
	// write list
	fprintf(fl, "%sdanss_%06d.root\n", listpref, n);
	last = n;
    }

    l = new TLegend(0.8, 0.8, 0.9, 1.0);
    for (i=0; i<sizeof(gr)/sizeof(gr[0]); i++) {
	if (i != 1 && ((1 << (i-1)) & mask)) l->AddEntry(gr[i], leg[i], "p");
    }
    
    
    TH1D h("hdummy", tit[col], last-first+2, first-1, last+1);
    h.SetMinimum(min - (max-min)/20.);
    h.SetMaximum(max + (max-min)/20.);
    h.GetXaxis()->SetLabelSize(0.055);
    h.GetYaxis()->SetLabelSize(0.055);
    gStyle->SetOptStat(0);
    h.DrawCopy();
    for (i=0; i<sizeof(gr)/sizeof(gr[0]); i++) {
	if (gr[i]->GetN() > 0) gr[i]->Draw("p");
    }
    l->Draw();
    
    fclose(fl);
    fclose(f);
}