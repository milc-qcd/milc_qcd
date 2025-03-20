/* The Plaquette gauge action */

#ifdef GAUGE_ACTION_PART1
#undef MAX_LENGTH
#undef MAX_NUM
#define NREPS 1
#define NLOOP 1
#define MAX_LENGTH 4
#define MAX_NUM 6
#endif

#ifdef GAUGE_ACTION_PART2
    static int loop_ind[NLOOP][MAX_LENGTH] = {
    { XUP, YUP, XDOWN, YDOWN }
    };
    static int loop_length_in[NLOOP] = {4};

    for(j=0;j<NLOOP;j++){
	loop_num[j] = 0;
	loop_length[j] = loop_length_in[j];
	for(i=0;i<NREPS;i++){
	    loop_coeff[j][i] = 0.0;
	}
    }

    /* set up the loop coefficients */
    loop_coeff[0][0]= 1.0;
    strcpy(gauge_action_description,"\"Single plaquette action\"");
    node0_printf("Single plaquette action\n");
#endif
