#include <string.h>
#include <stdlib.h>

void power(void)
{

    int i, j, n, sn;
    struct {
	unsigned long time;
	double power;
    } pow[550];
    FILE * f;
    char str[256];
    char * tok;
    double pnorm, s;


    unsigned long periods[][2] = {
	{1475218800, 1475398300},
	{1475398300, 1475577800},
	{1475578800, 1475851400},
//	{1475852400, 1476094100},
	{1476095100, 1476260900},
	{1476261900, 1476457100},
	{1476458100, 1476696500},
	{1476697500, 1476875600}};


    
    int np = sizeof(periods)/sizeof(periods[0]);
    
    struct {
	double power;
	int n;
    } avpow[20];

    memset(avpow, '\0', sizeof(avpow));
    
/*
    double data[][2] = {
         {44.32, 0.91},
         {44.32, 0.91},
         {42.72, 0.90},
//         {48.8,  5.5},
         {59.97, 1.00},
         {27.02, 0.79},
         {3.74,  0.60},
         {2.90,  0.74}};
*/

    double data[][2] = {
	{43.6, 0.9},
	{43.5, 0.9},
	{41.4, 0.9},
	{59.3, 1.0},
	{25.2, 0.8},
	{ 2.1, 0.6},
	{ 1.8, 0.7}};

         
    double const up = 1.0;
    double const mid = (11.7/10.7)*(11.7/10.7);
    double const down = (12.7/10.7)*(12.7/10.7);
    double const raised = (12.47/10.7)*(12.47/10.7);

    double norm[] = { raised, raised, down, up, mid, down, up };
    
    for (j=0; j < np; j++) {
	data[j][0] *= norm[j];
	data[j][1] *= norm[j];
    }
    
    f = fopen("b4power.txt", "rt");
    if (!f) return;
    fgets(str, sizeof(str), f);
    
    for (i=0, n=0, j=0; i<sizeof(pow)/sizeof(pow[0]); i++) {
        fgets(str, sizeof(str), f);
	if (feof(f)) break;
        tok = strtok(str, "\t ");
        pow[n].time = strtol(tok, NULL, 0) ;	// sutract 3 hours to convert to GMT
        tok = strtok(NULL, "\t ");
        pow[n].power = strtod(tok, NULL); 
	
	if (j < np && pow[n].time > periods[j][0]) {
	    if (pow[n].time < periods[j][1]) {
		avpow[j].power += pow[n].power;
		avpow[j].n ++;
	    } else {
		avpow[j].power /= (double)avpow[j].n;
		j++;
	    }
	}
        n++;
    }

    fclose(f);

//   normalizing by first 4 points 
    // average power
    for (j=0, s=0.0, sn=0; j<4; j++) {
	s += avpow[j].power*avpow[j].n;
	sn += avpow[j].n;
    }
    pnorm = s / (double)sn;
    for (j=0, s=0.0, sn=0; j<4; j++) {
	s += data[j][0];
	sn ++;
    }    
    pnorm /= (s / (double)sn);
    
    printf("Normalization=%.5f mHz/MW (%.2f MW/mHz)\n", 1./pnorm, pnorm);

    TH1D* hpow = new TH1D("hpow", "Block 4 power", n, pow[0].time, pow[n-1].time);
    TGraphErrors* apow = new TGraphErrors(np);
    TGraphErrors* dpow = new TGraphErrors(np);

    for (i=0; i<n; i++) hpow->SetBinContent(i+1, pow[i].power);
    for (j=0; j < np; j++) {
	apow->SetPoint(j, (periods[j][0] + periods[j][1])/2., avpow[j].power);
	apow->SetPointError(j, (periods[j][1] - periods[j][0])/2., 0);
    }
    for (j=0; j < np; j++) {
	dpow->SetPoint(j, (periods[j][0] + periods[j][1])/2., data[j][0]*pnorm);
	dpow->SetPointError(j, 0., data[j][1]*pnorm);
    }

    gStyle->SetTimeOffset(10800);	// -3 hours to convert back from GMT
    gStyle->SetEndErrorSize(7);
    gStyle->SetOptStat(0);
    hpow->GetXaxis()->SetTimeDisplay(1);
    hpow->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    hpow->GetXaxis()->SetNdivisions(710);
    hpow->GetXaxis()->SetTitle("Date day/month hr:min");
    hpow->SetLineColor(kBlue);
    hpow->SetLineWidth(2);
    apow->SetLineColor(kRed);
    apow->SetLineWidth(3);
    dpow->SetLineColor(kGreen+2);
    dpow->SetLineWidth(3);
    dpow->SetMarkerColor(kGreen+2);
    dpow->SetMarkerStyle(20);
    dpow->SetMarkerSize(1.5);
    

    
    hpow->Draw();
    apow->Draw("P");
    dpow->Draw("P");


}