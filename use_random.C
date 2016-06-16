void use_random(char *fname, char *rname)
{
	TFile fs(fname);
	if (!fs.IsOpen()) {
		printf("Can not open file %s\n", fname);
		return;
	}
	TFile fr(rname);
	if (!fr.IsOpen()) {
		printf("Can not open file %s\n", rname);
		return;
	}
	TTree *ts = fs.Get("DanssPair");
	if (!ts) {
		printf("No tree DanssPair in file %s\n", fname);
		return;
	}
	TTree *tr = fr.Get("DanssPair");
	if (!tr) {
		printf("No tree DanssPair in file %s\n", rname);
		return;
	}
//		Standard cuts
	TCut cs("EventsBetween == 0 && gtFromPrevious > 50 && gtToNext > 100 && gtDiff > 0.6");
	gStyle->SetOptStat(10);
	gStyle->SetOptFit(1);

//		Distance
	TH1D hds("hds", "Positron to Neutron distance;cm", 40, 0, 160);
	TH1D hdr("hdr", "Positron to Neutron distance;cm", 40, 0, 160);
	TH1D hdc("hdc", "Positron to Neutron distance;cm", 40, 0, 160);
	
	hds.SetLineColor(kBlue);
	hdr.SetLineColor(kRed);
	hdc.SetLineColor(kGreen);

	ts->Project("hds", "Distance", cs);
	tr->Project("hdr", "Distance", cs);

	hds.Sumw2();
	hdr.Sumw2();
	hdc.Add(&hds, &hdr, 1.0, -1.0);

	TCanvas *cd = new TCanvas("CD", "Distance", 800, 800);
	hds.DrawCopy("hist");
	hdr.DrawCopy("hist,sames");
	hdc.DrawCopy("same");
	
	TLegend lg(0.6, 0.75, 0.8, 0.9);
	lg.AddEntry("hds", "Events", "l");
	lg.AddEntry("hdr", "Random background", "l");
	lg.AddEntry("hdc", "Signal", "l");
	lg.Draw();
	
	cd->Update();
//	cd->SaveAs("r_distance.png");

//		Time
	TH1D hts("hts", "Positron to Neutron time;us", 50, 0, 50);
	TH1D htr("htr", "Positron to Neutron time;us", 50, 0, 50);
	TH1D htc("htc", "Positron to Neutron time;us", 50, 0, 50);
	
	hts.SetLineColor(kBlue);
	htr.SetLineColor(kRed);
	htc.SetLineColor(kGreen);

	ts->Project("hts", "gtDiff", cs);
	tr->Project("htr", "gtDiff", cs);

	hts.Sumw2();
	htr.Sumw2();
	htc.Add(&hts, &htr, 1.0, -1.0);

	TF1 f2exp("f2exp", "expo(0) - expo(2) + pol0(4)", 5);

	TCanvas *ct = new TCanvas("CT", "gtDiff", 800, 800);
	htc.Fit("f2exp", "", "0", 1, 50);
	hts.DrawCopy("hist");
	htr.DrawCopy("hist,sames");
	htc.DrawCopy("sames");
	
	lg.Draw();
	
	ct->Update();
//	ct->SaveAs("r_gtDiff.png");

//		NeutronHits
	TH1D hhns("hhns", "Hits in Neutron capture", 20, 0, 20);
	TH1D hhnr("hhnr", "Hits in Neutron capture", 20, 0, 20);
	TH1D hhnc("hhnc", "Hits in Neutron capture", 20, 0, 20);
	
	hhns.SetLineColor(kBlue);
	hhnr.SetLineColor(kRed);
	hhnc.SetLineColor(kGreen);

	ts->Project("hhns", "NeutronHits", cs);
	tr->Project("hhnr", "NeutronHits", cs);

	hhns.Sumw2();
	hhnr.Sumw2();
	hhnc.Add(&hhns, &hhnr, 1.0, -1.0);

	TCanvas *chn=new TCanvas("CNH", "NeutronHits", 800, 800);
	hhns.DrawCopy("hist");
	hhnr.DrawCopy("hist,sames");
	hhnc.DrawCopy("same");
	
	lg.Draw();
	
	chn->Update();
//	chn->SaveAs("r_NeutronHits.png");

//		NeutronEnergy
	TH1D hhes("hhes", "Energy in Neutron capture", 50, 0, 10);
	TH1D hher("hher", "Energy in Neutron capture", 50, 0, 10);
	TH1D hhec("hhec", "Energy in Neutron capture", 50, 0, 10);
	
	hhes.SetLineColor(kBlue);
	hher.SetLineColor(kRed);
	hhec.SetLineColor(kGreen);

	ts->Project("hhes", "NeutronEnergy", cs);
	tr->Project("hher", "NeutronEnergy", cs);

	hhes.Sumw2();
	hher.Sumw2();
	hhec.Add(&hhes, &hher, 1.0, -1.0);

	TCanvas *che = new TCanvas("CHE", "NeutronEnergy", 800, 800);
	hhes.DrawCopy("hist");
	hher.DrawCopy("hist,sames");
	hhec.DrawCopy("same");
	
	lg.Draw();
	
	che->Update();
//	che->SaveAs("r_NeutronEnergy.png");

//		GammaHits
	TH1D hpns("hpns", "Hits beyond positron cluster", 10, 0, 10);
	TH1D hpnr("hpnr", "Hits beyond positron cluster", 10, 0, 10);
	TH1D hpnc("hpnc", "Hits beyond positron cluster", 10, 0, 10);
	
	hpns.SetLineColor(kBlue);
	hpnr.SetLineColor(kRed);
	hpnc.SetLineColor(kGreen);

	ts->Project("hpns", "AnnihilationGammas", cs);
	tr->Project("hpnr", "AnnihilationGammas", cs);

	hpns.Sumw2();
	hpnr.Sumw2();
	hpnc.Add(&hpns, &hpnr, 1.0, -1.0);

	TCanvas *cpn=new TCanvas("CPN", "AnnihilationGammas", 800, 800);
	hpns.DrawCopy("hist");
	hpnr.DrawCopy("hist,sames");
	hpnc.DrawCopy("same");
	
	lg.Draw();
	
	cpn->Update();
//	cpn->SaveAs("r_AnnihilationGammas.png");

//		GammaEnergy
	TH1D hpes("hpes", "Energy beyond positron cluster", 20, 0, 2);
	TH1D hper("hper", "Energy beyond positron cluster", 20, 0, 2);
	TH1D hpec("hpec", "Energy beyond positron cluster", 20, 0, 2);
	
	hpes.SetLineColor(kBlue);
	hper.SetLineColor(kRed);
	hpec.SetLineColor(kGreen);

	ts->Project("hpes", "AnnihilationEnergy", cs);
	tr->Project("hper", "AnnihilationEnergy", cs);

	hpes.Sumw2();
	hper.Sumw2();
	hpec.Add(&hpes, &hper, 1.0, -1.0);

	TCanvas *cpe= new TCanvas("CPE", "annihilationEnergy", 800, 800);
	hpes.DrawCopy("hist");
	hper.DrawCopy("hist,sames");
	hpec.DrawCopy("same");
	
	lg.Draw();
	
	cpe->Update();
//	cpe->SaveAs("r_AnnihilationEnergy.png");

//		Criteria
	hdc.Divide(&hds);
	htc.Divide(&hts);
	hhnc.Divide(&hhns);
	hhec.Divide(&hhes);
	hpnc.Divide(&hpns);
	hpec.Divide(&hpes);
	
	TFile fres("norm.root", "RECREATE");
	fres.cd();
	hdc.Write();
	htc.Write();
	hhnc.Write();
	hhec.Write();
	hpnc.Write();
	hpec.Write();
	fres.Close();		

	fs.Close();
	fr.Close();
}

