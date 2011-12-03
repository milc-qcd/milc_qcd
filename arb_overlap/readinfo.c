/************* readinfo.c *****************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

#define MAX_REC_LEN 128
int read_gauge_info_i(FILE* fp, char *key, int* value)
{
    char *gleich = "=";
    char line[MAX_REC_LEN];
    char *pos;
    int success=0;

    
    
    rewind(fp);
    while (!feof(fp))
    {
	if (fgets(line, MAX_REC_LEN, fp))
	{
	    if (strstr(line,key) && (pos = strstr(line, gleich)))
	    {
		pos++;
		sscanf(pos,"%i",value);
		node0_printf("FOUND %s to be %i \n",key, *value);
		success=1;
		fflush(stdout);
	    }
	} else {
	    break;
	}
    }
    return success;
}


int read_gauge_info_d(FILE* fp, char *key, double* value)
{
    char *gleich = "=";
    char line[MAX_REC_LEN];
    char *pos;
    int success=0;
    Real val;

    
    
    rewind(fp);
    while (!feof(fp))
    {
	if (fgets(line, MAX_REC_LEN, fp))
	{
	    if (strstr(line,key) && (pos = strstr(line, gleich)))
	    {
		pos++;
		sscanf(pos,"%e",&val);
		*value=val;
		node0_printf("FOUND %s to be %f \n",key, *value);
		success=1;
		fflush(stdout);
	    }
	} else {
	    break;
	}
    }
    return success;
}
