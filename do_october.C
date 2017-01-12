#define N_OCTOBER 7
void do_october(void)
{
	const int t_period[N_OCTOBER][2] = {
		{1475218800, 1475398300},
		{1475398300, 1475577800},
		{1475578800, 1475851400},
		{1476095100, 1476260900},
		{1476261900, 1476457100},
		{1476458100, 1476696500},
		{1476697500, 1476875600},
		{1476876600, 1477052400},
		{1477053000, 
	};
	
	int i;
	TH1D *h[N_OCTOBER];
	char str[64];
	double val, err;
	
	for (i=0; i<N_OCTOBER; i++) {
		sprintf(str, "PositronEnergy_p%d", i+1);
		h[i] = spectr3("root/5469_6694.root", "root/5469_6694_r.root", str, t_period[i][0], t_period[i][1], 0.025);
		val = h[i]->IntegralAndError(1, 35, err);
		printf("Position %d: Freq = %7.1f +- %4.1f mHz\n", i+1, val, err);
	}
}
