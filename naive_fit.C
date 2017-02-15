void naive_fit(TH1 *h)
{
	TF1 *oFun = new TF1("oFun", "[0] * (1 - [1] * sin(1.27*[2]*12.7/(x+1.804)) * sin(1.27*[2]*12.7/(x+1.804))) / (1 - [1] * sin(1.27*[2]*10.7/(x+1.804)) * sin(1.27*[2]*10.7/(x+1.804)))");
	oFun->SetParNames("Const", "sin^{2}(2#theta)", "#Delta m^{2}");
	oFun->SetRange(0.1, 10);
	oFun->SetParameters(0.71, 0, 2.5);
	oFun->FixParameter(0, 0.71);
	oFun->SetLineColor(kBlue);
	oFun->SetLineWidth(2);
	h->Fit(oFun, "", "", 1, 7);
//	TF1 *oFun1 = new TF1("oFun1", "[0] * (1 - [1] * 0.5*(cos(2.54*[2]*10.7/(x+1.804)) - cos(2.54*[2]*12.7/(x+1.804))))");
//	oFun1->SetParNames("Const", "sin^{2}(2#theta)", "#Delta m^{2}");
//	oFun1->SetParameters(0.71, 0, 0.5);
//	oFun1->FixParameter(0, 0.71);
//	oFun1->SetLineColor(kRed);
//	oFun1->SetLineWidth(2);
//	h->Fit(oFun1, "+", "sames", 1, 7);
}
