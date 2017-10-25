#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define DEVIATION	(6*50)

int main(int argc, char **argv)
{
	const int PosVal[3] = {0, 6*990, 6*1980};
	char str[2049];
	int i, j, k, irc;
	int iFirst, iLast;
	FILE *f;
	struct {
		int len;
		int num;
		int ip;
		int type;
		int utime;
		int pos[4];
	} posrec;
	int p[4];
	
	if (argc < 4) {
		printf("Usage: %s data_directory first_num last_num\n", argv[0]);
		return 100;
	}
	
	iFirst = strtol(argv[2], NULL, 10);
	iLast = strtol(argv[3], NULL, 10);
	
	for(i = iFirst; i <= iLast; i++) {
		sprintf(str, "bzcat %s/danss_data_%6.6d.data.bz2 2>/dev/null", argv[1], i);
		printf("%6.6d ", i);
		f = popen(str, "r");
		if (!f) {
			printf("  -1\n");
			continue;
		}
		irc = fread(&posrec, sizeof(posrec), 1, f);
		if (irc != 1) {
			printf("   0\n");
			pclose(f);
			continue;
		}
		if (posrec.len != sizeof(posrec) || posrec.num || posrec.ip != 0x7F000001 || posrec.type != 0x2000000) {
			printf("   1\n");
			pclose(f);
			continue;
		}
		memset(p, 0, sizeof(p));
		for(j=0; j<4; j++) {
			for (k=0; k<3; k++) if (abs(posrec.pos[j] - PosVal[k]) < DEVIATION) p[j] = k+2;
			printf(" %5d[%1d] ", posrec.pos[j], p[j]);
		}
		if (!p[0] || p[1] != p[0] || p[2] != p[0] || p[3] != p[0]) {
			printf("   1\n");
			pclose(f);
			continue;
		}
		printf("   %d\n", p[0]);
		pclose(f);
	}
	
	return 0;
}
