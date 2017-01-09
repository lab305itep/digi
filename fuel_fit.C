//	MC generated fuel positron spectra
//	x[0] - energy, par[0] - U235, par[1] - U238, par[2] - Pu239, par[3] - Pu241
double FuelFunction(double *x, double *par)
{
	const double fuel_coef[4][9] = {
		{24.4944, -655.916, 5292.02, -3672.99, 1135.48, -195.083, 19.4828, -1.06507, 0.0247332},	// U235
		{19.6536, -524.634, 3896.56, -2200.03, 503.498, -55.1004, 2.50978, 0.00797526, -0.00299389},	// U238
		{32.4051, -925.359, 7073.58, -5452.01, 1884.54, -361.221, 39.9934, -2.40466, 0.0609398},	// Pu239
		{33.7827, -851.737, 5946.52, -4116.52, 1251.52, -207.387, 19.6074, -0.999286, 0.0213963}};	// Pu241
	double E;
	int i, j;
	double S, V;

	E = x[0];
 
	S = 0;
	for (j=0; j<4; j++) {
		V = par[j];
		for (i=0; i<9; i++) {
			S += V * fuel_coef[j][i];
			V *= E;
		}
	}
	return S;
}

//	Ratio for different U235/Pu239 mixture to the reactor start
//	x[0] - energy, par[0] - Pu239 amount, Pu241 and U238 - constant, U235 adding to unit
//	
double FuelRatio(double *x, double *par)
{
	const double U235_0  = 0.69;
	const double U238    = 0.07;
	const double PU239_0 = 0.20;
	const double PU241   = 0.04;
	
	double E = x[0];
	double PU239 = par[0];
	double U235 = 1 - PU239 - U238 - PU241;
	
	double ppar[4];
	double ref;
	double now;
	
	ppar[0] = U235;
	ppar[1] = U238;
	ppar[2] = PU239;
	ppar[3] = PU241;
	now = FuelFunction(x, ppar);
	
	ppar[0] = U235_0;
	ppar[1] = U238;
	ppar[2] = PU239_0;
	ppar[3] = PU241;
	ref = FuelFunction(x, ppar);
	
	return now/ref;
}
