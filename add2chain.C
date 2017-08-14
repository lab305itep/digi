void Add2Chain(TChain *ch, int run_from, int run_to, int mask, const char *f_template)
{
	const char stat_file_name[] = "stat_all.txt";
	FILE *f;
	char str[1024];
	char fname[1024];
	char *ptr;
	int i, num, cond, N;

	f = fopen(stat_file_name, "rt");
	if (!f) {
		printf("Stat file %s not found: %m\n", stat_file_name);
		return -1;
	}

	N = 0;
	for(i=0;;i++) {
		ptr = fgets(str, sizeof(str), f);
		if (!ptr || feof(f)) break;
		ptr = strtok(str, " \t");
		if (!isdigit(ptr[0])) continue;
		num = strtol(ptr, NULL, 10);
		ptr = strtok(NULL, " \t");
		if (!ptr) {
			printf("Strange record at line %d\n", i);
			continue;
		}
		cond = strtol(ptr, NULL, 10);
		if (num < run_from) continue;
		if (num > run_to) break;
		if (cond <= 0) continue;
		if (mask & (1 << (cond - 1))) {
			sprintf(fname, f_template, num);
			ch->AddFile(fname);
			N++;
		}
	}
	fclose(f);
	printf("%d runs found.\n", N);
}
