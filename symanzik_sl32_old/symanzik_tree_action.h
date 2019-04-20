/* Tree level Symanzik improved gauge action */
#ifdef GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#define NREPS 1
#define NLOOP 2
#define MAX_LENGTH 6
#define MAX_NUM 12
#endif

#ifdef GAUGE_ACTION_PART2
    static int loop_ind[NLOOP][MAX_LENGTH] = {
    { XUP, YUP, XDOWN, YDOWN, NODIR, NODIR },
    { XUP, XUP, YUP, XDOWN, XDOWN , YDOWN},
    };
    static int loop_length_in[NLOOP] = {4,6};

    for(j=0;j<NLOOP;j++){
	loop_num[j] = 0;
	loop_length[j] = loop_length_in[j];
	for(i=0;i<NREPS;i++){
	    loop_coeff[j][i] = 0.0;
	}
    }

    /* set up the loop coefficients */
    /* Loop coefficients from Urs */
    loop_coeff[0][0]= 1.0;
    loop_coeff[1][0]=  -1.00/(20.0*u0*u0);
    strcpy(gauge_action_description,"\"Symanzik 1x1 + 1x2 action\"");
    node0_printf("Symanzik 1x1 + 1x2 action\n");
#endif
