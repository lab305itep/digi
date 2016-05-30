void make_tcalib(char *fname)
{
	TFile *fIn;
	TH1D *hDt;
	FILE *fOut;
	TH1D *hEntries[50];
	TH1D *hMean[50];
	TH1D *hSigma[50];
	char outName[1024];
	char *ptr;
	char str[1024];
	char strs[128];
	int i, j;
	double mean, sigma, entries;
	double left, right;
	
	fIn = new TFile(fname, "UPDATE");
	if (!fIn->IsOpen()) {
		printf("Can not open file %s\n", fname);
		return;
	}
	strncpy(outName, fname, sizeof(outName)-10);
	ptr = strstr(outName, ".root");
	if (!ptr) ptr = &outName[strlen(OutName)];
	*ptr = '\0';
	strcat(outName, ".calib");
	fOut = fopen(outName, "wt");
	if (!fOut) goto fin;

	for (i=0; i<50; i++) {
		sprintf(strs, "hEntries%2.2d", i);
		sprintf(str, "Number of entries per channel. Module %d", i);
		hEntries[i] = new TH1D(strs, str, 64, 0, 64);
		sprintf(strs, "hMean%2.2d", i);
		sprintf(str, "Mean time per channel. Module %d", i);
		hMean[i] = new TH1D(strs, str, 64, 0, 64);
		sprintf(strs, "hSigma%2.2d", i);
		sprintf(str, "RMS per channel. Module %d", i);
		hSigma[i] = new TH1D(strs, str, 64, 0, 64);
		for (j=0; j<64; j++) {
			sprintf(str, "hDt%2.2dx%2.2d", i, j);
			hDt = (TH1D*) fIn->Get(str);
			if (!hDt) continue;
			entries = hDt->GetEntries();
			hEntries[i]->SetBinContent(j+1, entries);
			if (entries < 1000) continue;
			mean = hDt->GetMean();
			sigma = hDt->GetRMS();
			if (sigma > 15) continue;
			if (fabs(mean) > 150) continue;
			left = mean - 5*sigma;
			if (left < -200) left = -200;
			right = mean + 5*sigma;
			if (right > 200) right = 200;
			hDt->GetXaxis()->SetRangeUser(left, right);
			mean = hDt->GetMean();
			sigma = hDt->GetRMS();
			hMean[i]->SetBinContent(j+1, mean);
			hSigma[i]->SetBinContent(j+1, sigma);
			fprintf(fOut, "Channel=%2.2d.%2.2d  DT=%6.1f\n", i, j, mean);
		}
		hEntries[i]->Write();
		hMean[i]->Write();
		hSigma[i]->Write();
	}
	fIn->Write();
fin:
	fIn->Close();
	fclose(fOut);
}

void sum_tcalib(char *fna, char *fnb, char *fnc)
{
	FILE *fa;
	FILE *fb;
	FILE *fc;
	char *stra[1024];
	char *strb[1024];
	char *ptra;
	char *ptrb;
	int ia, ib, k;
	double dta, dtb;

	fa = fopen(fna, "rt");
	fb = fopen(fnb, "rt");
	fc = fopen(fnc, "wt");
	if (!fa || !fb || !fc) {
		printf("Can not open all requested files\n");
		goto fin;
	}
	
	k = 0;
	for(;;) {
		ptra = fgets(stra, sizeof(stra), fa);
		ptrb = fgets(strb, sizeof(strb), fb);
		if (!ptra || !ptrb) break;
		if (stra[0] == '*') continue;	// Comment ?

		ptra = strstr(stra, "Channel=");
		ptrb = strstr(strb, "Channel=");
		if (!ptra || !ptrb) continue;		// strange string
		ptra += strlen("Channel=");
		ptrb += strlen("Channel=");
		ia = 100 * (strtod(ptra, 0) + 0.002);
		ib = 100 * (strtod(ptrb, 0) + 0.002);
		if (ia != ib) {
			printf("Bad synchronization of files: %5.2f(a) != %5.2f(b)\n", ia / 100.0, ib / 100.0);
			break;
		}
		ptra = strstr(stra, "DT=");
		ptrb = strstr(strb, "DT=");
		if (!ptra || !ptrb) continue;		// strange string
		ptra += strlen("DT=");
		ptrb += strlen("DT=");
		dta = strtod(ptra, 0);
		dtb = strtod(ptrb, 0);
		fprintf(fc, "Channel=%2.2d.%2.2d  DT=%6.1f\n", ia / 100, ia % 100, dta + dtb);
		k++;

	}
	printf("Done: %d channels\n", k);
fin:
	fclose(fa);
	fclose(fb);
	fclose(fc);
}

