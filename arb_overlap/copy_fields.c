/**************** copy_fields.c *************************/
/* MIMD version 7 */
#include "arb_ov_includes.h"

///* copy a gauge field - an array of four su3_matrices */
//void gauge_field_copy(field_offset src,field_offset dest) {
//register int i,dir,src2,dest2;
//register site *s;
//    FORALLSITES(i,s){
//	src2=src; dest2=dest;
//        for(dir=XUP;dir<=TUP; dir++){
//	    su3mat_copy( (su3_matrix *)F_PT(s,src2),
//		(su3_matrix *)F_PT(s,dest2) );
//	    src2 += sizeof(su3_matrix);
//	    dest2 += sizeof(su3_matrix);
//	}
//    }
//}

/* exchange gauge fields - an array of four su3_matrices */
void gauge_field_change(field_offset src,field_offset dest) {
register int i,dir,src2,dest2;
register site *s;
su3_matrix tmp;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    tmp=* (su3_matrix *)F_PT(s,dest2);
	    su3mat_copy( (su3_matrix *)F_PT(s,src2),
		(su3_matrix *)F_PT(s,dest2) );
	    su3mat_copy( &tmp,
		(su3_matrix *)F_PT(s,src2) );
	    
	    src2 += sizeof(su3_matrix);
	    dest2 += sizeof(su3_matrix);
	}
    }
}

/* copy a momenta - an array of four ah_matrices */
void mom_copy(field_offset src,field_offset dest) {
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    *(( anti_hermitmat *)F_PT(s,dest2))=* (( anti_hermitmat *)F_PT(s,src2)) ;
	    src2 += sizeof( anti_hermitmat);
	    dest2 += sizeof( anti_hermitmat);
	}
    }
}
