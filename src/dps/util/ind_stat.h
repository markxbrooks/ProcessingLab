/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_stat.h:   */

#ifndef IND_STAT_H

void calc_dr_ds(struct rs_coord *, struct rs_coord *, int , 
	   float *, float *, float *);

int num_ok_hkl(struct hkl_fract *, int , struct hkl_int *, float);


#define IND_STAT_H
#endif
