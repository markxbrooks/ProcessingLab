/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_transform.h:  Headerfile for ind_transform.f */

#undef RS_2_XYZ_DEBUG
#undef CUT_XYZ_DEBUG

void xyz_2_rs_long(struct vector *, int, detector_setting, 
	      float, float , float ,struct rs_coord *, float *, int *, struct vector *);
void xyz_2_rs(struct vector *, int, detector_setting, 
	      float, float , float ,struct rs_coord *);
void rs_2_xyz(int, detector_setting,
	     struct rs_coord *, struct vector *, int called_from_mosflm);
int cut_xyz(int, struct vector *, float, float, struct vector *);
int xyz_2_rs_single_use_mos(struct vector, detector_setting,
	 float, float *, float *, float *,struct rs_coord *, float *, int *, struct vector *);
int xyz_2_rs_single(struct vector, detector_setting,
	 float, float *, float *, float *,struct rs_coord *, float *, int *, struct vector *);
