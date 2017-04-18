#include "HPainter.h"
#define NHISTS 18
void background_calc(const char *fname = "background_plots.root", int run_first = 5808, int run_last = 15028)
{
	char strs[128];
	char strl[1024];
	const char titles[NHISTS][16] = {"gtDiff", "R1", "R2", "RZ", "PX", "PY", "PZ", "NX", "NY", "NZ", "NE", "NH", "NR", "PH", "PMaxE", "AH", "AE", "PF"};
	const char titlel[NHISTS][256] = {"Time from positron to neutron;us", "Distance between positron and neutron, 2D;cm", 
		"Distance between positron and neutron, 3D;cm", "Distance between positron and neutron, Z;cm", 
		"Positron vertex X;cm", "Positron vertex Y;cm", "Positron vertex Z;cm", 
		"Neutron vertex X;cm", "Neutron vertex Y;cm", "Neutron vertex Z;cm", "Neutron capture energy", "Neutron capture SiPM hits",
		"Neutron capture average photon flight",
		"Hits in positron cluster", "Maximum energy in one hit", "Number of hits out of cluster", "Energy out of cluster;MeV",
		"Edge hit in prompt event"};
	TH1D *h[NHISTS][2];
	int i, j;
	
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct;

	TFile *fRoot = new TFile(fname, "RECREATE");
	for (i=0; i<NHISTS; i++) for (j=0; j<2; j++) {
		sprintf(strs, "h%s%c", titles[i], (j) ? 'C' : 'N');
		sprintf(strl, "%s: %s", (j) ? "Cosmic" : "Neutrino", titlel[i]);
		switch(i) {
		case 0:		// gtDiff
			h[i][j] = new TH1D(strs, strl, 50, 0, 50.0);
			break;
		case 1:		// R1, R2
		case 2:
			h[i][j] = new TH1D(strs, strl, 40, 0, 160.0);
			break;
		case 3:		// RZ
			h[i][j] = new TH1D(strs, strl, 100, -50.0, 50.0);
			break;
		case 4:		// PX, PY, NX, NY
		case 5:
		case 7:
		case 8:
			h[i][j] = new TH1D(strs, strl, 25, 0, 100.0);
			break;
		case 6:		// PZ, NZ
		case 9:
			h[i][j] = new TH1D(strs, strl, 100, 0, 100.0);
			break;
		case 10:	// NE
			h[i][j] = new TH1D(strs, strl, 45, 3.0, 12.0);
			break;
		case 11:	// NH
			h[i][j] = new TH1D(strs, strl, 20, 0, 20.0);
			break;
		case 12:	// NR
			h[i][j] = new TH1D(strs, strl, 15, 0, 60.0);
			break;
		case 13:	// PH
			h[i][j] = new TH1D(strs, strl, 10, 0, 10.0);
			break;
		case 14:	// PMaxE
			h[i][j] = new TH1D(strs, strl, 40, 0, 8.0);
			break;
		case 15:	// AH
			h[i][j] = new TH1D(strs, strl, 20, 0, 20.0);
			break;
		case 16:	// AE
			h[i][j] = new TH1D(strs, strl, 20, 0, 4.0);
			break;
		case 17:	// PF
			h[i][j] = new TH1D(strs, strl, 2, 0, 2);
			break;
		}
	}

	HPainter *hp = new HPainter(0x801E, run_first, run_last);
	hp->SetFile(fRoot);

	for (j=0; j<2; j++) {
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN;
		hp->Project(h[0][j], "gtDiff", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma && cN && !cRXY;
		hp->Project(h[1][j], "Distance", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma && cN && cRXY;
		hp->Project(h[2][j], "Distance", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cPe && cGamma && cN;
		hp->Project(h[3][j], "DistanceZ", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cY && cZ && cR && cPe && cGamma && cN && "PositronX[0] >= 0";
		hp->Project(h[4][j], "PositronX[0] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cZ && cR && cPe && cGamma && cN && "PositronX[1] >= 0";
		hp->Project(h[5][j], "PositronX[1] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cR && cPe && cGamma && cN && "PositronX[2] >= 0";
		hp->Project(h[6][j], "PositronX[2] + 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cY && cZ && cR && cPe && cGamma && cN && "NeutronX[0] >= 0";
		hp->Project(h[7][j], "NeutronX[0] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cZ && cR && cPe && cGamma && cN && "NeutronX[1] >= 0";
		hp->Project(h[8][j], "NeutronX[1] + 2.0", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cR && cPe && cGamma && cN && "NeutronX[2] >= 0";
		hp->Project(h[9][j], "NeutronX[2] + 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[10][j], "NeutronEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[11][j], "NeutronHits", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma;
		hp->Project(h[12][j], "NeutronRadius", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN;
		hp->Project(h[13][j], "PositronHits", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN;
		hp->Project(h[14][j], "MaxHitEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cN;
		hp->Project(h[15][j], "AnnihilationGammas", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cN;
		hp->Project(h[16][j], "AnnihilationEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN;
		hp->Project(h[17][j], "(PositronFlags[0] && 0x3FF0000) ? 1.5 : 0.5", (j) ? (ct && !cVeto) : (ct && cVeto));
	}
	
	fRoot->cd();
	for (i=0; i<NHISTS; i++) for (j=0; j<2; j++) h[i][j]->Write();
	delete hp;
	fRoot->Close();
}


void background_draw_all(const char *rootname = "background_plots.root")
{
	char strs[128];
	char strl[1024];
	const char titles[NHISTS][16] = {"gtDiff", "R1", "R2", "RZ", "PX", "PY", "PZ", "NX", "NY", "NZ", "NE", "NH", "NR", "PH", "PMaxE", "AH", "AE", "PF"};
	const char suffix[3][10] = {"N-rand", "N-diff", "C-diff"};
	const Color_t color[3] = {kGreen+2, kBlue, kRed};
	const int marker[3] = {kFullStar, kFullCircle, kFullCross};
	TH1D *h[NHISTS][3];
	TH1D *hf[NHISTS][3];
	int i, j, k, kl, ku;
	double hMax;
	int iMax;
	double total, totale;
	double frac, frace;
	char pdfname[1024];
	char *ptr;

	strcpy(pdfname, rootname);
	ptr = strstr(pdfname, ".root");
	if (ptr) {
		strcpy(ptr, ".pdf");
	} else {
		strcat(pdfname, ".pdf");
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	
	TFile *fRoot = new TFile(rootname);
	if (!fRoot->IsOpen()) {
		printf("root-file not found  - run background_calc() first!\n");
		return;
	}
	for (i=0; i<NHISTS; i++) {
		for (j=0; j<3; j++) {
			sprintf(strs, "h%s%s", titles[i], suffix[j]);
			h[i][j] = (TH1D*) fRoot->Get(strs);
			if (!h[i][j]) {
				printf("%s not found  - rerun background_calc() to create all hists!\n", rootname);
				fRoot->Close();
				return;
			}
			h[i][j]->SetLineWidth(2);
			h[i][j]->SetLineColor(color[j]);
			h[i][j]->SetMarkerColor(color[j]);
			h[i][j]->SetMarkerStyle(marker[j]);
			h[i][j]->GetYaxis()->SetLabelSize(0.05);
			h[i][j]->SetMinimum(0);
			h[i][j]->GetYaxis()->SetTitle("");
//			Fractions
			sprintf(strs, "h%s%s-f", titles[i], suffix[j]);
			hf[i][j] = (TH1D *) h[i][j]->Clone(strs);
			hf[i][j]->Clear();
			switch(i) {
			case 0:		// from top
			case 1:
			case 2:
			case 12:
			case 15:
			case 16:
			case 17:
				hf[i][j]->SetTitle("Cutting from top");
				total = h[i][j]->IntegralAndError(1, hf[i][j]->GetNbinsX(), totale);
				for (k=0; k<hf[i][j]->GetNbinsX(); k++) {
					frac = h[i][j]->IntegralAndError(1, k+1, frace);
					hf[i][j]->SetBinContent(k+1, frac/total);
					hf[i][j]->SetBinError(k+1, (frac/total) * 
						TMath::Sqrt((frace*frace/(frac*frac) + totale*totale/(total*total))));
				}
				break;
			case 3:		// from both sides
			case 4:
			case 5:
			case 6:
			case 7:
			case 8:
			case 9:
				hf[i][j]->SetTitle("Cutting from both sides");
				total = h[i][j]->IntegralAndError(1, hf[i][j]->GetNbinsX(), totale);
				for (k=0; k<hf[i][j]->GetNbinsX(); k++) {
					if (k+1 <= hf[i][j]->GetNbinsX()/2) {
						kl = k+1;
						ku = hf[i][j]->GetNbinsX() - k;
					} else {
						ku = k+1;
						kl = hf[i][j]->GetNbinsX() - k;
					}
					frac = h[i][j]->IntegralAndError(kl, ku, frace);
					hf[i][j]->SetBinContent(k+1, frac/total);
					hf[i][j]->SetBinError(k+1, (frac/total) * 
						TMath::Sqrt((frace*frace/(frac*frac) + totale*totale/(total*total))));
				}
				break;
			case 10: 	// from bottom
			case 11:
			case 13:
			case 14:
				hf[i][j]->SetTitle("Cutting from bottom");
				total = h[i][j]->IntegralAndError(1, hf[i][j]->GetNbinsX(), totale);
				for (k=0; k<hf[i][j]->GetNbinsX(); k++) {
					frac = h[i][j]->IntegralAndError(k+1, hf[i][j]->GetNbinsX(), frace);
					hf[i][j]->SetBinContent(k+1, frac/total);
					hf[i][j]->SetBinError(k+1, (frac/total) * 
						TMath::Sqrt((frace*frace/(frac*frac) + totale*totale/(total*total))));
				}
				break;
			}
			hf[i][j]->SetMaximum(1.0);
			hf[i][j]->SetMinimum(0.0);
		}
	}
	
	TCanvas *cv = new TCanvas("CV", "Background plots", 1200, 900);
	TLegend *lg = new TLegend(0.7, 0.8, 0.95, 0.95);
	lg->AddEntry(h[0][0], "Random", "LP");
	lg->AddEntry(h[0][1], "Neutrino", "LP");
	lg->AddEntry(h[0][2], "Cosmic", "LP");
	
	sprintf(strl, "%s[", pdfname);
	cv->SaveAs(strl);
	
	for (i=0; i<NHISTS; i++) {
		cv->Clear();
		cv->Divide(2, 1);

		cv->cd(1);
		hMax = 0;
		iMax = 0;
		for (j=0; j<3; j++) if (h[i][j]->GetMaximum() > hMax) {
			hMax = h[i][j]->GetMaximum();
			iMax = j;
		}
		h[i][iMax]->Draw();
		for (j=0; j<3; j++) if (j != iMax) h[i][j]->Draw("same");
		lg->Draw();
		
		cv->cd(2);
		for (j=0; j<3; j++) hf[i][j]->Draw((j) ? "same" : "");
		
		cv->SaveAs(pdfname);
	}
	
	sprintf(strl, "%s]", pdfname);
	cv->SaveAs(strl);
	fRoot->Close();
}

void background_calcgt(const char *fname = "background_plotsgt.root", int run_first = 5808, int run_last = 15028)
{
	char strs[128];
	char strl[1024];
	TH1D *h[3][2];
	int i, j;
	
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe("PositronEnergy > 1");
        TCut cR60("Distance < 60");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut cN3("NeutronEnergy > 3.0");
        TCut ct;

	TFile *fRoot = new TFile(fname, "RECREATE");
	for (i=0; i<3; i++) for (j=0; j<2; j++) {
		sprintf(strs, "hgtDiff%c%c", 'A'+i, (j) ? 'C' : 'N');
		sprintf(strl, "Time from positron to neutron: %s;us", (j) ? "Cosmic" : "Neutrino");
		h[i][j] = new TH1D(strs, strl, 50, 0, 50.0);
	}

	HPainter *hp = new HPainter(0x801E, run_first, run_last);
	hp->SetFile(fRoot);

	for (j=0; j<2; j++) {
		ct = cIso && cX && cY && cZ && cR60 && cRZ && cPe && cGamma && cN3;
		hp->Project(h[0][j], "gtDiff", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN3;
		hp->Project(h[1][j], "gtDiff", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cR && cPe && cGamma && cN;
		hp->Project(h[2][j], "gtDiff", (j) ? (ct && !cVeto) : (ct && cVeto));
	}
	
	fRoot->cd();
	for (i=0; i<3; i++) for (j=0; j<2; j++) h[i][j]->Write();
	delete hp;
	fRoot->Close();
}

void background_draw_gt(const char *rootname = "background_plotsgt.root")
{
	char strs[128];
	char strl[1024];
	char title[3][128] = {
		"Old cuts: R < 60 cm && NeutronEnergy > 3.0 MeV",
		"Diatance < 45 cm for 2D and Distance < 55 cm for 3 D, NeutronEnergy > 3.0 MeV",
		"Diatance < 45 cm for 2D and Distance < 55 cm for 3 D, NeutronEnergy > 3.5 MeV"
	};
	const char suffix[3][10] = {"N-rand", "N-diff", "C-diff"};
	const Color_t color[3] = {kGreen+2, kBlue, kRed};
	const int marker[3] = {kFullStar, kFullCircle, kFullCross};
	TH1D *h[3][3];
	int i, j, k, kl, ku;
	double hMax;
	int iMax;
	double total, totale;
	double frac, frace;
	char pdfname[1024];
	char *ptr;

	strcpy(pdfname, rootname);
	ptr = strstr(pdfname, ".root");
	if (ptr) {
		strcpy(ptr, ".pdf");
	} else {
		strcat(pdfname, ".pdf");
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetLabelSize(0.05);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	
	TFile *fRoot = new TFile(rootname);
	if (!fRoot->IsOpen()) {
		printf("root-file not found  - run background_calcgt() first!\n");
		return;
	}
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			sprintf(strs, "hgtDiff%c%s", 'A'+i, suffix[j]);
			h[i][j] = (TH1D*) fRoot->Get(strs);
			if (!h[i][j]) {
				printf("%s not found  - rerun background_calc() to create all hists!\n", rootname);
				fRoot->Close();
				return;
			}
			h[i][j]->SetLineWidth(2);
			h[i][j]->SetLineColor(color[j]);
			h[i][j]->SetMarkerColor(color[j]);
			h[i][j]->SetMarkerStyle(marker[j]);
			h[i][j]->GetYaxis()->SetLabelSize(0.05);
			h[i][j]->SetMinimum(0);
			h[i][j]->GetYaxis()->SetTitle("");
			h[i][j]->SetTitle(title[i]);
		}
	}
	TCanvas *cv = new TCanvas("CV", "Background plots", 1200, 900);
	TLegend *lg = new TLegend(0.7, 0.8, 0.95, 0.95);
	lg->AddEntry(h[0][0], "Random", "LP");
	lg->AddEntry(h[0][1], "Neutrino", "LP");
	lg->AddEntry(h[0][2], "Cosmic", "LP");
	
	sprintf(strl, "%s[", pdfname);
	cv->SaveAs(strl);
	
	for (i=0; i<3; i++) {
		cv->Clear();
		hMax = 0;
		iMax = 0;
		for (j=0; j<3; j++) if (h[i][j]->GetMaximum() > hMax) {
			hMax = h[i][j]->GetMaximum();
			iMax = j;
		}
		h[i][iMax]->Draw();
		for (j=0; j<3; j++) if (j != iMax) h[i][j]->Draw("same");
		lg->Draw();
		
		cv->SaveAs(pdfname);
	}
	
	sprintf(strl, "%s]", pdfname);
	cv->SaveAs(strl);
	fRoot->Close();
}

void background_calcpe(const char *fname = "background_plotspe.root", int run_first = 5808, int run_last = 15028)
{
	char strs[128];
	char strl[1024];
	TH1D *h[3][2];
	int i, j;
	
//		Main cuts
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut cT20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
        TCut cPe("PositronEnergy > 1");
        TCut cR60("Distance < 60");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut cN3("NeutronEnergy > 3.0");
//        TCut cPMT("fabs(SiPmCleanEnergy[0] - PmtCleanEnergy[0]) < 5");
        TCut ct;

	TFile *fRoot = new TFile(fname, "RECREATE");
	for (i=0; i<3; i++) for (j=0; j<2; j++) {
		sprintf(strs, "hPosEnergy%c%c", 'A'+i, (j) ? 'C' : 'N');
		sprintf(strl, "Positron Energy %s;MeV", (j) ? "Cosmic" : "Neutrino");
		h[i][j] = new TH1D(strs, strl, 44, 1.0, 12.0);
	}

	HPainter *hp = new HPainter(0x801E, run_first, run_last);
	hp->SetFile(fRoot);

	for (j=0; j<2; j++) {
		ct = cIso && cX && cY && cZ && cT20 && cR60 && cRZ && cPe && cGamma && cN3 && cPMT;
		hp->Project(h[0][j], "PositronEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cT20 && cR && cPe && cGamma && cN3 && cPMT;
		hp->Project(h[1][j], "PositronEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
		ct = cIso && cX && cY && cZ && cT20 && cR && cPe && cGamma && cN && cPMT;
		hp->Project(h[2][j], "PositronEnergy", (j) ? (ct && !cVeto) : (ct && cVeto));
	}
	
	fRoot->cd();
	for (i=0; i<3; i++) for (j=0; j<2; j++) h[i][j]->Write();
	delete hp;
	fRoot->Close();
}

void background_draw_pe(const char *rootname = "background_plotspe.root")
{
	char strs[128];
	char strl[1024];
	char title[3][128] = {
		"Old cuts: R < 60 cm && NeutronEnergy > 3.0 MeV",
		"Diatance < 45 cm for 2D and Distance < 55 cm for 3 D, NeutronEnergy > 3.0 MeV",
		"Diatance < 45 cm for 2D and Distance < 55 cm for 3 D, NeutronEnergy > 3.5 MeV"
	};
	const char suffix[3][10] = {"N-rand", "N-diff", "C-diff"};
	const Color_t color[3] = {kGreen+2, kBlue, kRed};
	const int marker[3] = {kFullStar, kFullCircle, kFullCross};
	TH1D *h[3][3];
	TH1D *hz[3][3];
	TH1D *hr[3][2];
	int i, j, k, kl, ku;
	double hMax;
	int iMax;
	double total, totale;
	double frac, frace;
	char pdfname[1024];
	char *ptr;
	TPad *pd;
	TVirtualPad *pv;

	strcpy(pdfname, rootname);
	ptr = strstr(pdfname, ".root");
	if (ptr) {
		strcpy(ptr, ".pdf");
	} else {
		strcat(pdfname, ".pdf");
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	
	TFile *fRoot = new TFile(rootname);
	if (!fRoot->IsOpen()) {
		printf("root-file not found  - run background_calcgt() first!\n");
		return;
	}
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			sprintf(strs, "hPosEnergy%c%s", 'A'+i, suffix[j]);
			h[i][j] = (TH1D*) fRoot->Get(strs);
			if (!h[i][j]) {
				printf("%s not found  - rerun background_calcpe() to create all hists!\n", rootname);
				fRoot->Close();
				return;
			}
			h[i][j]->SetLineWidth(2);
			h[i][j]->SetLineColor(color[j]);
			h[i][j]->SetMarkerColor(color[j]);
			h[i][j]->SetMarkerStyle(marker[j]);
			h[i][j]->GetXaxis()->SetTitleSize(0.06);
			h[i][j]->GetYaxis()->SetTitleSize(0.06);
			h[i][j]->GetXaxis()->SetLabelSize(0.06);
			h[i][j]->GetYaxis()->SetLabelSize(0.06);
			h[i][j]->SetMinimum(0);
			h[i][j]->GetYaxis()->SetTitle("");
			h[i][j]->SetTitle(title[i]);
			sprintf(strs, "%s_z", h[i][j]->GetName());
			hz[i][j] = (TH1D*)h[i][j]->Clone(strs);
			hz[i][j]->SetTitle("");
			hz[i][j]->GetXaxis()->SetRange(29, 44);
		}
		for (j=0; j<2; j++) {
			sprintf(strs, "hRPES%c%c", (j) ? 'C': 'R', 'A'+i);
			hr[i][j] = (TH1D*)h[i][0]->Clone(strs);
			hr[i][j]->Clear();
			hr[i][j]->Divide(h[i][2*j], h[i][1]);
			hr[i][j]->SetLineWidth(2);
			hr[i][j]->SetLineColor(color[j]);
			hr[i][j]->SetMarkerColor(color[j]);
			hr[i][j]->SetMarkerStyle(marker[j]);
			hr[i][j]->GetXaxis()->SetTitleSize(0.06);
			hr[i][j]->GetYaxis()->SetTitleSize(0.06);
			hr[i][j]->GetXaxis()->SetLabelSize(0.06);
			hr[i][j]->GetYaxis()->SetLabelSize(0.06);
			hr[i][j]->SetMinimum(0);
			hr[i][j]->SetMaximum(5);
			hr[i][j]->GetYaxis()->SetTitle("");
			hr[i][j]->SetTitle(title[i]);
		}
	}
	TCanvas *cv = new TCanvas("CV", "Background plots", 1200, 900);
	TLegend *lg1 = new TLegend(0.7, 0.8, 0.95, 0.95);
	lg1->AddEntry(h[0][0], "Random", "LPE");
	lg1->AddEntry(h[0][1], "Neutrino", "LPE");
	lg1->AddEntry(h[0][2], "Cosmic", "LPE");
	lg1->SetTextSize(0.04);
	TLegend *lg2 = new TLegend(0.25, 0.8, 0.65, 0.9);
	lg2->AddEntry(hr[0][0], "Random/Neutrino", "LPE");
	lg2->AddEntry(hr[0][1], "Cosmic/Neutrino", "LPE");
	lg2->SetTextSize(0.035);
	
	sprintf(strl, "%s[", pdfname);
	cv->SaveAs(strl);
	
	for (i=0; i<3; i++) {
		cv->Clear();
		cv->Divide(2, 1);
		
		pv = cv->cd(1);
		pv->SetFillStyle(4000);

		hMax = 0;
		iMax = 0;
		for (j=0; j<3; j++) if (h[i][j]->GetMaximum() > hMax) {
			hMax = h[i][j]->GetMaximum();
			iMax = j;
		}
		h[i][iMax]->Draw();
		for (j=0; j<3; j++) if (j != iMax) h[i][j]->Draw("same");
		lg1->Draw();

		pv->cd();
		pd = new TPad("PD", "", 0.3, 0.35, 0.99, 0.79);
		pd->SetTopMargin(0);
		pd->SetRightMargin(0);
		pd->Draw();
		pd->cd();
		hMax = 0;
		iMax = 0;
		for (j=0; j<3; j++) if (hz[i][j]->GetMaximum() > hMax) {
			hMax = hz[i][j]->GetMaximum();
			iMax = j;
		}
		hz[i][iMax]->Draw();
		for (j=0; j<3; j++) if (j != iMax) hz[i][j]->Draw("same");
		pv->SetFillStyle(4000);
		pd->Draw();

		cv->cd(2);
		hr[i][0]->Draw();
		hr[i][1]->Draw("same");
		lg2->Draw();
		
		cv->SaveAs(pdfname);
	}
	
	sprintf(strl, "%s]", pdfname);
	cv->SaveAs(strl);
	fRoot->Close();
}

