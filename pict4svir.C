#include "HPainter2.h"
TH1D *pict4svir_fromVeto(int mask, int run_from, int run_to)
{
//		Main cuts
//	TCut cVeto("gtFromVeto > 60");
//	TCut cMuonA("gtFromVeto == 0");
//	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
//	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct = cIso && cX && cY && cZ && c20 && cGamma && cGammaMax && cPe && cR && cN;

        TH1D *h = new TH1D("hFromVeto", "Time from muon;us", 100, 0, 100);
        HPainter2 *hp = new HPainter2(mask, run_from, run_to);
        hp->Project(h, "gtFromVeto", ct);
        h->Scale(hp->GetUpTime()/1000.0);
        return h;
}

TH1D *pict4svir_posiFast(int mask, int run_from, int run_to)
{
//		Main cuts
//	TCut cVeto("gtFromVeto > 60");
	TCut cMuonA("gtFromVeto == 0");
//	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
//	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct = cIso && cX && cY && cZ && c20 && cGamma && cGammaMax && cPe && cR && cN && cMuonA;

        TH1D *h = new TH1D("hPosiFast", "Fast muon background;MeV", 40, 0, 10);
        HPainter2 *hp = new HPainter2(mask, run_from, run_to);
        hp->Project(h, "PositronEnergy", ct);
        h->Scale(hp->GetUpTime()/1000.0);
        return h;
}

TH1D *pict4svir_posiSlow(int mask, int run_from, int run_to)
{
//		Main cuts
//	TCut cVeto("gtFromVeto > 60");
//	TCut cMuonA("gtFromVeto == 0");
	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
//	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct = cIso && cX && cY && cZ && c20 && cGamma && cGammaMax && cPe && cR && cN && cMuonB;

        TH1D *h = new TH1D("hPosiSlow", "Slow muon background;MeV", 40, 0, 10);
        HPainter2 *hp = new HPainter2(mask, run_from, run_to);
        hp->Project(h, "PositronEnergy", ct);
        h->Scale(hp->GetUpTime()/1000.0);
        return h;
}

TH1D *pict4svir_sumFast(int mask, int run_from, int run_to)
{
//		Main cuts
//	TCut cVeto("gtFromVeto > 60");
	TCut cMuonA("gtFromVeto == 0");
//	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
//	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct = cIso && cX && cY && cZ && c20 && cGamma && cGammaMax && cPe && cR && cN && cMuonA;

        TH1D *h = new TH1D("hSumFast", "Fast muon background;MeV", 40, 0, 10);
        HPainter2 *hp = new HPainter2(mask, run_from, run_to);
        hp->Project(h, "(SiPmCleanEnergy[0]+PmtCleanEnergy[0])/2", ct);
        h->Scale(hp->GetUpTime()/1000.0);
        return h;
}

TH1D *pict4svir_sumSlow(int mask, int run_from, int run_to)
{
//		Main cuts
//	TCut cVeto("gtFromVeto > 60");
//	TCut cMuonA("gtFromVeto == 0");
	TCut cMuonB("gtFromVeto > 0 && gtFromVeto <= 60");
	TCut cIso("(gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0) || (gtFromPrevious == gtFromVeto)");
//	TCut cShower("gtFromVeto > 200 || DanssEnergy < 300");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cRXY("PositronX[0] >= 0 && PositronX[1] >= 0 && NeutronX[0] >= 0 && NeutronX[1] >= 0");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.8 && AnnihilationGammas <= 10");
	TCut cGammaMax("AnnihilationMax < 0.8");
        TCut cPe("PositronEnergy > 1");
        TCut cR1("Distance < 45");
        TCut cR2("Distance < 55");
        TCut cRZ("fabs(DistanceZ) < 40");
        TCut cR = cR2 && (cRXY || cR1) && cRZ;
        TCut cN("NeutronEnergy > 3.5");
        TCut ct = cIso && cX && cY && cZ && c20 && cGamma && cGammaMax && cPe && cR && cN && cMuonB;

        TH1D *h = new TH1D("hSumSlow", "Slow muon background;MeV", 40, 0, 10);
        HPainter2 *hp = new HPainter2(mask, run_from, run_to);
        hp->Project(h, "(SiPmCleanEnergy[0]+PmtCleanEnergy[0])/2", ct);
        h->Scale(hp->GetUpTime()/1000.0);
        return h;
}


void pict4svir(int mask, int run_from, int run_to)
{
	TH1D *h;
	TFile *f = new TFile("pict4svir.root", "RECREATE");
	
	h = pict4svir_fromVeto(mask, run_from, run_to);
	f->cd();
	h->Write();

	h = pict4svir_posiFast(mask, run_from, run_to);
	f->cd();
	h->Write();
	
	h = pict4svir_posiSlow(mask, run_from, run_to);
	f->cd();
	h->Write();

	h = pict4svir_sumFast(mask, run_from, run_to);
	f->cd();
	h->Write();
	
	h = pict4svir_sumSlow(mask, run_from, run_to);
	f->cd();
	h->Write();
	
	f->Close();
}
