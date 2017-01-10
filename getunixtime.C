#include <time.h>

time_t getunixtime(const char *str)
{
	struct tm t;
	
	memset(&t, 0, sizeof(t));
	strptime(str, "%d.%m.%Y %H:%M", &t);
	return mktime(&t);
}
