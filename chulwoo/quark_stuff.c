#ifdef  ASQ_OPTIMIZED_FORCE
void u_shift_fermion_special(su3_vector *src, su3_vector *dest, int dir, 
 msg_tag **tag, int start) {
  su3_vector *tmpvec ; 
  register site *s ;
  register int i ;
  int *nei;
  
  if(GOES_FORWARDS(dir)) /* forward shift */
    {
      nei = neighbor(dir);
#ifndef PPC440QCDOC
  u_shift_fermion_plus_called++;
#endif
      if(start)
      *tag = start_gather_field(src, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0]);
      else
           restart_gather_field_special(src, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0], *tag);
#if 0
      FORALLSITES(i,s)
      if(nei[i]!=NOWHERE){
	mult_su3_mat_vec(&(s->link[dir]),(su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      }
#else
     all_mult_su3_mat_vec(dir,gen_pt[0],dest,nei,1, sites_on_node, lattice);
#endif
      wait_gather(*tag);
#if 0
      FORALLSITES(i,s)
      if(nei[i]==NOWHERE){
	mult_su3_mat_vec(&(s->link[dir]),(su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      }
#else
     all_mult_su3_mat_vec(dir,gen_pt[0],dest,nei,0, sites_on_node, lattice);
#endif
    }
  else /* backward shift */
    {
      nei = neighbor(OPP_DIR(dir));

      tmpvec = (su3_vector *)emalloc( sites_on_node*sizeof(su3_vector) );
/*  parallel transport for non-local points */

#if 0
      FORALLSITES(i,s)
      if(nei[i]==NOWHERE){
	mult_adj_su3_mat_vec(&(s->link[OPP_DIR(dir)]),&(src[i]), &tmpvec[i]);
      }
#else
     all_mult_adj_su3_mat_vec(OPP_DIR(dir),src,tmpvec,nei,0,sites_on_node,
lattice);
#endif
      if(start)
      *tag = start_gather_field(tmpvec, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0]);
      else
           restart_gather_field_special(tmpvec, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0], *tag);
/*  parallel transport for local points */
#if 0
      FORALLSITES(i,s)
      if(nei[i]!=NOWHERE){
	mult_adj_su3_mat_vec(&(s->link[OPP_DIR(dir)]),&(src[i]), &tmpvec[i]);
      } 
#else
     all_mult_adj_su3_mat_vec(OPP_DIR(dir),src,tmpvec,nei,1,sites_on_node,lattice);
#endif
      wait_gather(*tag);
      /* copy the gen_pt to the dest */
      FORALLSITES(i,s)
	dest[i] = *(su3_vector *)gen_pt[0][i];
      efree(tmpvec) ;
    }

}

void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,Real coeff) {
  register site *s ;
  register int i ,ii,jj;  
  register Real tmp_coeff ,temp0,temp1,temp2,temp;
  register complex tcmplx0,tcmplx1;

  if(GOES_BACKWARDS(dir))
    {
      dir = OPP_DIR(dir) ; 
      coeff = -coeff ;
    }
  FORALLSITES(i,s){
    if(s->parity==ODD) 
      tmp_coeff = -coeff ;
    else
      tmp_coeff = coeff ;
#if 0
    uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
    su3_projector(&(back[i]), &(forw[i]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff, &tmat2 );
    make_anti_hermitian( &tmat2, &(s->mom[dir]) ); 
#else

    temp0 = back[i].c[0].imag * forw[i].c[0].real -
            back[i].c[0].real * forw[i].c[0].imag;
    temp1 = back[i].c[1].imag * forw[i].c[1].real -
            back[i].c[1].real * forw[i].c[1].imag;
    temp2 = back[i].c[2].imag * forw[i].c[2].real -
            back[i].c[2].real * forw[i].c[2].imag;
    temp = (temp0+temp1+temp2)*0.3333333333;
   
    s->mom[dir].m00im += tmp_coeff*(temp0-temp);
    s->mom[dir].m11im += tmp_coeff*(temp1-temp);
    s->mom[dir].m22im += tmp_coeff*(temp2-temp);
    
    CMUL_J(back[i].c[0],forw[i].c[1],tcmplx0);
    CMUL_J(back[i].c[1],forw[i].c[0],tcmplx1);
    s->mom[dir].m01.real += tmp_coeff*(tcmplx0.real -tcmplx1.real)*0.5;
    s->mom[dir].m01.imag += tmp_coeff*(tcmplx0.imag +tcmplx1.imag)*0.5;
    CMUL_J(back[i].c[0],forw[i].c[2],tcmplx0);
    CMUL_J(back[i].c[2],forw[i].c[0],tcmplx1);
    s->mom[dir].m02.real += tmp_coeff*(tcmplx0.real -tcmplx1.real)*0.5;
    s->mom[dir].m02.imag += tmp_coeff*(tcmplx0.imag +tcmplx1.imag)*0.5;
    CMUL_J(back[i].c[1],forw[i].c[2],tcmplx0);
    CMUL_J(back[i].c[2],forw[i].c[1],tcmplx1);
    s->mom[dir].m12.real += tmp_coeff*(tcmplx0.real -tcmplx1.real)*0.5;
    s->mom[dir].m12.imag += tmp_coeff*(tcmplx0.imag +tcmplx1.imag)*0.5;
#endif

}
