#define OUT_PDF "mc_ibd.pdf"
void draw_mc_ibd(void)
{
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetLabelSize(0.06);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetLineWidth(2);
	
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR("Distance < 60 && DistanceZ > -40 && DistanceZ < 40");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");

	TFile *fp = new TFile("/space/danss_root2/mc_IBD_Gd_check_transcode_newScale.root");
	TTree *tp = (TTree *) fp->Get("DanssPair");
	if (!tp) return;
	TFile *fOut = new TFile("mc_ibd.root", "RECREATE");
	
	TH1D *hPositronEnergy = new TH1D("hPositronEnergy", "MC - recovered positron energy;E, MeV", 35, 1, 8);
	TH1D *hgtDiff = new TH1D("hgtDiff", "Neutron capture time;t_{capture}, us", 50, 0, 50);
	TH1D *hDistance = new TH1D("hDistance", "Distance between positron and neutron vertexes;R, cm", 25, 0, 100);
	TH1D *hDistanceZ = new TH1D("hDistanceZ", "Distance between positron and neutron vertexes, Z-projection;R_{Z}, cm", 100, -50, 50);
	TH1D *hX = new TH1D("hX", "Positron vertex X;X, cm", 25, 0, 100);
	TH1D *hY = new TH1D("hY", "Positron vertex Y;Y, cm", 25, 0, 100);
	TH1D *hZ = new TH1D("hZ", "Positron vertex Z;Z, cm", 100, 0, 100);
	TH1D *hNeutronHits = new TH1D("hNeutronHits", "SiPm hits in neutron capture;N_{hits}", 20, 0, 20);
	TH1D *hNeutronEnergy = new TH1D("hNeutronEnergy", "Energy in neutron capture;E_{N}, MeV", 50, 0, 10);
	TH1D *hAnnihilationHits = new TH1D("hAnnihilationHits", "SiPm hits beyond positron cluster;N_{hits}", 20, 0, 20);
	TH1D *hAnnihilationEnergy = new TH1D("hAnnihilationEnergy", "Energy beyond positron cluster;E_{#gamma}, MeV", 36, 0, 2);
	
	tp->Project("hPositronEnergy", "PositronEnergy", cX && cY && cZ && cR && c20 && cGamma && cPe);
	tp->Project("hgtDiff", "gtDiff", cX && cY && cZ && cR && cGamma && cPe);
	tp->Project("hDistance", "Distance", cX && cY && cZ && c20 && cGamma && cPe);
	tp->Project("hDistanceZ", "DistanceZ", cX && cY && cZ && c20 && cGamma && cPe);
	tp->Project("hX", "PositronX[0]+2", cY && cZ && cR && c20 && cGamma && cPe && "PositronX[0]>= 0");
	tp->Project("hY", "PositronX[1]+2", cX && cZ && cR && c20 && cGamma && cPe && "PositronX[1]>= 0");
	tp->Project("hZ", "PositronX[2]+0.5", cX && cY && cR && c20 && cGamma && cPe && "PositronX[2]>= 0");
	tp->Project("hNeutronHits", "NeutronHits", cX && cY && cZ && cR && c20 && cGamma && cPe);
	tp->Project("hNeutronEnergy", "NeutronEnergy", cX && cY && cZ && cR && c20 && cGamma && cPe);
	tp->Project("hAnnihilationHits", "AnnihilationGammas", cX && cY && cZ && cR && c20 && cPe);
	tp->Project("hAnnihilationEnergy", "AnnihilationEnergy", cX && cY && cZ && cR && c20 && cPe);
	
	hPositronEnergy->GetYaxis()->SetLabelSize(0.06);
	hgtDiff->GetYaxis()->SetLabelSize(0.06);
	hDistance->GetYaxis()->SetLabelSize(0.06);
	hDistanceZ->GetYaxis()->SetLabelSize(0.06);
	hX->GetYaxis()->SetLabelSize(0.06);
	hY->GetYaxis()->SetLabelSize(0.06);
	hZ->GetYaxis()->SetLabelSize(0.06);
	hNeutronHits->GetYaxis()->SetLabelSize(0.06);
	hNeutronEnergy->GetYaxis()->SetLabelSize(0.06);
	hAnnihilationHits->GetYaxis()->SetLabelSize(0.06);
	hAnnihilationEnergy->GetYaxis()->SetLabelSize(0.06);
	
	TCanvas *cv = new TCanvas("CV", "MC_IBD", 1200, 800);
	
	cv->SaveAs(OUT_PDF "[");
	
	cv->Divide(2, 1);
	cv->cd(1);
	hPositronEnergy->Draw();
	cv->cd(2);
	hgtDiff->Draw();
	cv->SaveAs(OUT_PDF);
	
	cv->Clear();
	cv->Divide(2, 1);
	cv->cd(1);
	hDistance->Draw();
	cv->cd(2);
	hDistanceZ->Draw();
	cv->SaveAs(OUT_PDF);
	
	cv->Clear();
	cv->Divide(2, 1);
	cv->cd(1);
	hX->Draw();
	cv->cd(2);
	hY->Draw();
	cv->SaveAs(OUT_PDF);
	
	cv->Clear();
	cv->Divide(2, 1);
	cv->cd(1);
	hZ->Draw();
	cv->SaveAs(OUT_PDF);
	
	cv->Clear();
	cv->Divide(2, 1);
	cv->cd(1);
	hNeutronHits->Draw();
	cv->cd(2);
	hNeutronEnergy->Draw();
	cv->SaveAs(OUT_PDF);
	
	cv->Clear();
	cv->Divide(2, 1);
	cv->cd(1);
	hAnnihilationHits->Draw();
	cv->cd(2);
	hAnnihilationEnergy->Draw();
	cv->SaveAs(OUT_PDF);
	cv->SaveAs(OUT_PDF "]");
	
	fOut->cd();
	hPositronEnergy->Write();
	hgtDiff->Write();
	hDistance->Write();
	hDistanceZ->Write();
	hX->Write();
	hY->Write();
	hZ->Write();
	hNeutronHits->Write();
	hNeutronEnergy->Write();
	hAnnihilationHits->Write();
	hAnnihilationEnergy->Write();
	fOut->Close();
}
