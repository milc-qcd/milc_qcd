void set_I_from_field(QDP_Int *dest, field_offset src) { int i;
 site *s;
 QLA_Int *temp;
 temp = QDP_expose_I (dest);
 for(i=0,s=lattice;
     i<sites_on_node;
     i++,s++) {
   memcpy((void *)&temp[i], ((char *)( s ) + (src)), sizeof(QLA_Int));
 }
 QDP_reset_I (dest);
}
;
void set_R_from_field(QDP_D_Real *dest, field_offset src) {
  int i;
  site *s;
  QLA_Real *temp;
  temp = QDP_D_expose_R( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], ((char *)( s ) + (src)), sizeof(QLA_Real));
  }
  QDP_D_reset_R( dest );
}
;
void set_V_from_field(QDP_D3_ColorVector *dest, field_offset src) {
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_D3_expose_V( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], ((char *)( s ) + (src)), sizeof(QLA_ColorVector));
  }
  QDP_D3_reset_V( dest );
}
;
void set_M_from_field(QDP_D3_ColorMatrix *dest, field_offset src) {
  int i;
  site *s;
  QLA_ColorMatrix *temp;
  temp = QDP_D3_expose_M( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], ((char *)( s ) + (src)), sizeof(QLA_ColorMatrix));
  }
  QDP_D3_reset_M( dest );
}
;
;


void set_field_from_I(field_offset dest, QDP_Int *src) {
  int i;
  site *s;
  QLA_Int *temp;
  temp = QDP_expose_I (src);
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy(((char *)( s ) + (dest)), (void *)&temp[i], sizeof(QLA_Int));
  }
  QDP_reset_I (src);
}
;
void set_field_from_R(field_offset dest, QDP_D_Real *src) {
  int i;
  site *s;
  QLA_Real *temp;
  temp = QDP_D_expose_R( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy(((char *)( s ) + (dest)), (void *)&temp[i], sizeof(QLA_Real));
  }
  QDP_D_reset_R( src );
}
;
void set_field_from_V(field_offset dest, QDP_D3_ColorVector *src) {
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_D3_expose_V( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy(((char *)( s ) + (dest)), (void *)&temp[i], sizeof(QLA_ColorVector));
  }
  QDP_D3_reset_V( src );
}
;
void set_field_from_M(field_offset dest, QDP_D3_ColorMatrix *src) {
  int i;
  site *s;
  QLA_ColorMatrix *temp;
  temp = QDP_D3_expose_M( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy(((char *)( s ) + (dest)), (void *)&temp[i], sizeof(QLA_ColorMatrix));
  }
  QDP_D3_reset_M( src );
}
;
;


void set_I_from_temp(QDP_Int *dest, int *src) {
  int i;
  site *s;
  QLA_Int *temp;
  temp = QDP_expose_I (dest);
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_Int));
  }
  QDP_reset_I (dest);
}
;
void set_R_from_temp(QDP_D_Real *dest, double *src) {
  int i;
  site *s;
  QLA_Real *temp;
  temp = QDP_D_expose_R( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_Real));
  }
  QDP_D_reset_R( dest );
}
;
void set_V_from_temp(QDP_D3_ColorVector *dest, su3_vector *src) {
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_D3_expose_V( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_ColorVector));
  }
  QDP_D3_reset_V( dest );
}
;
void set_M_from_temp(QDP_D3_ColorMatrix *dest, su3_matrix *src) {
  int i;
  site *s;
  QLA_ColorMatrix *temp;
  temp = QDP_D3_expose_M( dest );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_ColorMatrix));
  }
  QDP_D3_reset_M( dest );
}
;
;


void set_temp_from_I(int *dest, QDP_Int *src) {
  int i;
  site *s;
  QLA_Int *temp;
  temp = QDP_expose_I (src);
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_Int));
  }
  QDP_reset_I (src);
}
;
void set_temp_from_R(double *dest, QDP_D_Real *src) {
  int i;
  site *s;
  QLA_Real *temp;
  temp = QDP_D_expose_R( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_Real));
  }
  QDP_D_reset_R( src );
}
;
void set_temp_from_V(su3_vector *dest, QDP_D3_ColorVector *src) {
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_D3_expose_V( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_ColorVector));
  }
  QDP_D3_reset_V( src );
}
;
void set_temp_from_M(su3_matrix *dest, QDP_D3_ColorMatrix *src) {
  int i;
  site *s;
  QLA_ColorMatrix *temp;
  temp = QDP_D3_expose_M( src );
  for(i=0,s=lattice;
      i<sites_on_node;
      i++,s++) {
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_ColorMatrix));
  }
  QDP_D3_reset_M( src );
}
;
;


void set4_I_from_temp(QDP_Int *dest[], int *src) {
  int i, dir;
  site *s;
  QLA_Int *temp;
  for(dir=0;
      dir<4;
      dir++) {
    temp = QDP_expose_I (dest[dir]);
    for(i=0,s=lattice;
	i<sites_on_node;
	i++,s++) {
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_Int));
    }
    QDP_reset_I (dest[dir]);
  }
}
;
void set4_R_from_temp(QDP_D_Real *dest[], double *src) {
  int i, dir;
  site *s;
  QLA_Real *temp;
  for(dir=0;
      dir<4;
      dir++) {
    temp = QDP_D_expose_R( dest[dir] );
    for(i=0,s=lattice;
	i<sites_on_node;
	i++,s++) {
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_Real));
    }
    QDP_D_reset_R( dest[dir] );
  }
}
;
void set4_V_from_temp(QDP_D3_ColorVector *dest[], su3_vector *src) {
  int i, dir;
  site *s;
  QLA_ColorVector *temp;
  for(dir=0;
      dir<4;
      dir++) {
    temp = QDP_D3_expose_V( dest[dir] );
    for(i=0,s=lattice;
	i<sites_on_node;
	i++,s++) {
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_ColorVector));
    }
    QDP_D3_reset_V( dest[dir] );
  }
}
;
void set4_M_from_temp(QDP_D3_ColorMatrix *dest[], su3_matrix *src) {
  int i, dir;
  site *s;
  QLA_ColorMatrix *temp;
  for(dir=0;
      dir<4;
      dir++) {
    temp = QDP_D3_expose_M( dest[dir] );
    for(i=0,s=lattice;
	i<sites_on_node;
	i++,s++) {
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_ColorMatrix));
    }
    QDP_D3_reset_M( dest[dir] );
  }
}
;
;

