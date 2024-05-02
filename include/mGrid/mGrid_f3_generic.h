#ifndef _MGRID_F3_GENERIC_H
#define _MGRID_F3_GENERIC_H

#define GRID_Real                       GRID_F_Real
#define GRID_asqtad_create_L_from_MILC	GRID_F3_asqtad_create_L_from_MILC
#define GRID_extract_MILC_from_L        GRID_F3_extract_MILC_from_L     
#define GRID_asqtad_destroy_L		GRID_F3_asqtad_destroy_L
#define GRID_asqtad_invert_multi        GRID_F3_asqtad_invert_multi
#define GRID_asqtad_invert_mixed        GRID_F3_asqtad_invert_mixed
#define GRID_asqtad_invert_block        GRID_F3_asqtad_invert_block
#define GRID_asqtad_invert_mixed_block  GRID_F3_asqtad_invert_mixed_block
#define GRID_asqtad_invert              GRID_F3_asqtad_invert
#define GRID_create_M_from_mat4         GRID_F3_create_M_from_mat4
#define GRID_create_V			GRID_F3_create_V
#define GRID_create_nV			GRID_F3_create_nV
#define GRID_destroy_V			GRID_F3_destroy_V
#define GRID_destroy_nV			GRID_F3_destroy_nV
#define GRID_create_V_from_vec		GRID_F3_create_V_from_vec
#define GRID_create_nV_from_vecs	GRID_F3_create_nV_from_vecs
#define GRID_extract_V_to_vec		GRID_F3_extract_V_to_vec
#define GRID_extract_nV_to_vecs		GRID_F3_extract_nV_to_vecs
#define GRID_hisq_aux_links             GRID_F3_hisq_aux_links
#define GRID_hisq_links                 GRID_F3_hisq_links
#define GRID_ColorMatrix                GRID_F3_ColorMatrix
#define GRID_ColorVector                GRID_F3_ColorVector
#define GRID_ColorVectorBlock           GRID_F3_ColorVectorBlock
#define GRID_FermionLinksAsqtad         GRID_F3_FermionLinksAsqtad

#define GRID_ColorVectorArray GRID_F3_ColorVectorArray
#define GRID_implicitly_restarted_lanczos GRID_F3_implicitly_restarted_lanczos
#define GRID_create_V_array GRID_F3_create_V_array
#define GRID_destroy_V_array GRID_F3_destroy_V_array
#define GRID_create_V_array_from_vec_array GRID_F3_create_V_array_from_vec_array
#define GRID_extract_V_array_to_vec_array GRID_F3_extract_V_array_to_vec_array

#endif /* _MGRID_F3_GENERIC_H */
