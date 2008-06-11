/*************  dslash_w_space.c *******************************/
/* MIMD version 7 */

#include "generic_wilson_includes.h"

/* Temporary work space for dslash_w_field and dslash_w_field_special */ 
static half_wilson_vector *htmp[8] ;

/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

void malloc_dslash_w_3D_temps(){
  int j;

  if(!temp_not_allocated)return;
  for( j = 0; j < 8; j++ ){
    htmp[j] =(half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(htmp[j] == NULL){
      printf("node %d can't malloc htmp[%d]\n",this_node,j);
      terminate(1);
    }
  }
  temp_not_allocated = 0 ;
}

void cleanup_dslash_w_3D_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<8;i++) {
      free(htmp[i]) ; 
    }
  temp_not_allocated=1 ;
}

void grow_add_three_wvecs2( wilson_vector *a, half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3, int sign, int sum ){
  int i;
  if(sum==0)
    {
      /* wp_grow( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESMINUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESMINUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESPLUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESPLUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
    }
  else
    {
      /*wp_grow_add( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TMI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TMI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TPI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TPI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
    }
  
  /* wp_grow_add( b2,a,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUM( a->d[2].c[i], b2->h[1].c[i]);
	CSUB( a->d[3].c[i], b2->h[0].c[i], a->d[3].c[i] );
      }
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUB( a->d[2].c[i], b2->h[1].c[i], a->d[2].c[i] );
	CSUM( a->d[3].c[i], b2->h[0].c[i]);
      }
    }
  
  /* wp_grow_add( b3,a,ZUP,sign); */
  
  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TMI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TPI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  else
    {
      /* case ZDOWN:*/
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TPI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TMI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  

}

void wp_shrink_3dir2( wilson_vector *a,  half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3, int sign ){
  register int i; /*color*/
  
/* wp_shrink( a,b1,XUP,sign); */

  if(sign==PLUS)
    {
    /* case XUP: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].real;
      }
    }
  else
    {
      /* case XDOWN: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].real;
      }
    }
  
  
  /*    wp_shrink( a,b2,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].imag;
      }
      
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].imag;
      }
    }
  
  /*    wp_shrink( a,b3,ZUP,sign); */

  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real - a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag + a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real + a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag - a->d[3].c[i].real;
      }
    }
  else
    {
      /* case ZDOWN: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real + a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag - a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real - a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag + a->d[3].c[i].real;
      }
      
    }

}

/*
	dslash_w_3D(src, dest, isign, l_parity);
Compute SUM_dirs_xyz ( 
    ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
)

*/

void dslash_w_3D_site( field_offset src, field_offset dest, int isign, 
		       int parity) {
half_wilson_vector hwvx,hwvy,hwvz;

register int i;
register site *s;
register int dir,otherparity=0;
msg_tag *tag[8];

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    malloc_dslash_w_3D_temps();

    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
      wp_shrink_3dir2( (wilson_vector *)F_PT(s,src), htmp[XUP]+i,
	    htmp[YUP]+i, htmp[ZUP]+i, isign);
      //if(i==259)printf("i = %d, (1+gx)G = %e\n", i, s->htmp[XUP].h[1].c[0].real);
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	tag[dir]=start_gather_field( htmp[dir], sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
      wp_shrink_3dir2( (wilson_vector *)F_PT(s,src),
		       &hwvx, &hwvy, &hwvz,-isign);
	mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, htmp[XDOWN]+i);
	mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, htmp[YDOWN]+i);
	mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, htmp[ZDOWN]+i);
	//if(i==258)printf("i = %d, U+(1-gx)G = %e\n", i, s->htmp[XDOWN].h[1].c[0].real);

    }

    for( dir=XUP; dir <= ZUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather_field(htmp[OPP_DIR(dir)], 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= ZUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	mult_su3_mat_hwvec( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	//if(i==2)printf("i = %d, received (should be i=259)(1+gx)G = %e\n", i, (*(half_wilson_vector * )(gen_pt[XUP][i])).h[1].c[0].real);

	/* printf("(%d %d %d %d) up[YUP] %e %e\n", s->x,s->y,s->z,s->t, s->link[YUP].e[0][0].real,s->link[YUP].e[0][0].imag); */
	grow_add_three_wvecs2( (wilson_vector *)F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= ZUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
      //if(i==2)printf("i = %d, received (should be i=258)(1-gx)G = %e\n", i, (*(half_wilson_vector * )(gen_pt[XDOWN][i])).h[1].c[0].real);
	grow_add_three_wvecs2( (wilson_vector *)F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w_space() ) */

void dslash_w_3D_field( wilson_vector *src, wilson_vector *dest, 
			int isign, int parity) {
half_wilson_vector hwvx,hwvy,hwvz;

register int i;
register site *s;
register int dir,otherparity=0;
msg_tag *tag[8];

  /* The calling program must clean up the temps! */
  malloc_dslash_w_3D_temps();

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
      wp_shrink_3dir2( src+i, htmp[XUP]+i, htmp[YUP]+i, htmp[ZUP]+i, isign);
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	tag[dir]=start_gather_field( htmp[dir], sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
      wp_shrink_3dir2( src+i, &hwvx, &hwvy, &hwvz,-isign);
	mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, &(htmp[XDOWN][i]));
	mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, &(htmp[YDOWN][i]));
	mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, &(htmp[ZDOWN][i]));
    }

    for( dir=XUP; dir <= ZUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather_field(htmp[OPP_DIR(dir)], 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= ZUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	mult_su3_mat_hwvec( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	grow_add_three_wvecs2( dest+i,
	    &hwvx, &hwvy, &hwvz, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= ZUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
      grow_add_three_wvecs2( dest+i,
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= ZUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w_space() ) */

