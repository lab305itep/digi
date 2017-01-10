void do_april(void)
{
	const unsigned int t_period[10][2] = {
		{1461229920, 1461321060},
		{1461321840, 1461403200},
		{1461403800, 1461488040},
		{1461488340, 1461587040},
		{1461587220, 1461664680},
		{1461665040, 1461748140},
		{1461748500, 1461853800},
		{1461854160, 1461920400},
		{1461920940, 1462006620},
		{1462007220, 1462076460}
	};

	int i;
	TH1D *h[10];
	char str[64];
	double err, val;
	
	for (i=0; i<10; i++) {
		sprintf(str, "PositronEnergy_p%d", i+1);
		h[i] = spectr3("root/2306_2834.root", "root/2306_2834_r.root", str, t_period[i][0], t_period[i][1]);
		val = h[i]->IntegralAndError(1, 35, err);
		printf("Position %d: Freq = %7.1f +- %4.1f mHz\n", i+1, val, err);
	}
}
