#ifndef _MGRID_D3_GENERIC_H
#define _MGRID_D3_GENERIC_H

#define GRID_Real                       GRID_D_Real
#define GRID_asqtad_create_L_from_MILC	GRID_D3_asqtad_create_L_from_MILC
#define GRID_asqtad_destroy_L		GRID_D3_asqtad_destroy_L
#define GRID_asqtad_invert_multi        GRID_D3_asqtad_invert_multi
#define GRID_asqtad_invert_block        GRID_D3_asqtad_invert_block
#define GRID_asqtad_invert              GRID_D3_asqtad_invert
#define GRID_create_V			GRID_D3_create_V
#define GRID_create_nV		        GRID_D3_create_nV
#define GRID_destroy_V		        GRID_D3_destroy_V
#define GRID_destroy_nV		        GRID_D3_destroy_nV
#define GRID_create_V_from_vec		GRID_D3_create_V_from_vec
#define GRID_create_nV_from_vecs	GRID_D3_create_nV_from_vecs
#define GRID_extract_V_to_vec		GRID_D3_extract_V_to_vec
#define GRID_extract_nV_to_vecs	        GRID_D3_extract_nV_to_vecs
#define GRID_ColorVector                GRID_D3_ColorVector
#define GRID_ColorVectorBlock           GRID_D3_ColorVectorBlock
#define GRID_FermionLinksAsqtad         GRID_D3_FermionLinksAsqtad

#define GRID_ColorVectorArray GRID_D3_ColorVectorArray
#define GRID_implicitly_restarted_lanczos GRID_D3_implicitly_restarted_lanczos
#define GRID_create_V_array GRID_D3_create_V_array
#define GRID_destroy_V_array GRID_D3_destroy_V_array
#define GRID_create_V_array_from_vec_array GRID_D3_create_V_array_from_vec_array
#define GRID_extract_V_array_to_vec_array GRID_D3_extract_V_array_to_vec_array

#endif /* _MGRID_D3_GENERIC_H */
