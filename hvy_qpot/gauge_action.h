/* Loop tables for glueball operators */
#ifdef GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#define NREPS 1
#define NLOOP 4
#define MAX_LENGTH 6
#define MAX_NUM 16
#define SPATIAL_ONLY
#endif

#ifdef GAUGE_ACTION_PART2
    static int loop_ind[NLOOP][MAX_LENGTH] = {
    { XUP, YUP, XDOWN, YDOWN, NODIR, NODIR },  /* plaquette */      
    { XUP, XUP, YUP, XDOWN, XDOWN , YDOWN},    /* 1x2 loop */       
    { XUP, YUP, ZUP, YDOWN, XDOWN , ZDOWN},    /* chair */          
    { XUP, YUP, ZUP, XDOWN, YDOWN , ZDOWN},    /* equator of cube */
    };
    static int loop_length_in[NLOOP] = {4,6,6,6};

    for(j=0;j<NLOOP;j++){
	loop_num[j] = 0;
	loop_length[j] = loop_length_in[j];
	for(i=0;i<NREPS;i++){
	    loop_coeff[j][i] = 0.0;
	}
    }

    /* Loop coefficients from Urs */
    loop_coeff[0][0]= 1.0;
    loop_coeff[1][0]= 1.0;
    loop_coeff[2][0]= 1.0;
    loop_coeff[3][0]= 1.0;
    strcpy(gauge_action_description,"\"Glueball 1x1, 1x2, chair, oval\"");
    node0_printf("Glueball ops 1x1, 1x2, chair, oval\n");

#endif
