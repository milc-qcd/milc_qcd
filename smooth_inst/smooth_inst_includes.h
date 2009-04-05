#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/dirs.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

/* hex 00010203 for Real 00010223 for double */
#define TOPO_VERSION_NUMBER (66051+8*(sizeof(Real)-4))

int setup(void);
int readin(const int prompt);
int initial_set(void);
void make_loop_table2(void);
void smooth(void);
void instanton_density(void);
void load_fatlinks(void);
void save_topo(char *filenam);
void char_num(int dig[13], int chr[2], int *ch, const int length);
void path(int *dir,int *sign, const int length);
