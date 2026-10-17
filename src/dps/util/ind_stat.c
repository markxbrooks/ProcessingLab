/*=======================================================================
 * All files in the distribution of the DPS system are Copyright
 * 1996 by the Computational Biology group in the Department of Biological
 * Sciences at Purdue University.  All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that this entire copyright notice is duplicated in all such
 * copies, and that any documentation, announcements, and other materials
 * related to such distribution and use acknowledge that the software was
 * developed by the Computational Biology group in the Department of
 * Biological Sciences at Purdue University, W. Lafayette, IN by Ingo
 * Steller and Michael G. Rossmann. No charge may be made for copies,
 * derivations, or distributions of this material without the express
 * written consent of the copyright holder.  Neither the name of the
 * University nor the names of the authors may be used to endorse or
 * promote products derived from this material without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE.
 *======================================================================*/

/*=====================================================================*
 *                                                                     *
 *                         Data Processing Suit                        *
 *                                                                     *
 *                            Utility-Library                          *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                            File: ind_stat.c                         *
 *                                                                     *
 *=====================================================================*/

#include <math.h>
#include <stdio.h>
#include "util.h"

/* calc_dr_dr(): On two given experimental and predicted r/s lists	*
 *		 calculate the mean differences dr, ds.			*/

void calc_dr_ds(struct rs_coord *rs, struct rs_coord *rs_pred, int number, 
float *dr, float *ds, float *sum_sqr)
{
	int i;	/* counter */
	int div;	/* Number of summands */
	float help_dr, help_ds;

	*dr = 0.0;
	*ds = 0.0;
	*sum_sqr=0.0;
	div = 0;
	/* Go through the lists */
	for(i=0; i < number; i++) {
		/* If predict r and s values are 0.0, then there was
                   a problem in calcualting them. See also xyz_2_rs() 	*/
		if ((rs_pred[i].r != 0.0) || (rs_pred[i].s != 0.0)) {
			help_dr=rs[i].r-rs_pred[i].r;
			help_ds=rs[i].s-rs_pred[i].s;
			*dr = *dr + fabs(help_dr);
			*ds = *ds + fabs(help_ds);
			*sum_sqr=*sum_sqr+help_dr*help_dr+help_ds*help_ds;
			div++;
		}
	}
	*dr = *dr / div;
	*ds = *ds / div;
}

/* num_hkl_ok(): Number of hkl with fractional part less then a given
                 value. A list of integer and no integer hkl must be given.
		 See also xyz_2_hkl()					*/

int num_ok_hkl(struct hkl_fract hkl[], int number_hkl, 
struct hkl_int hkl_calc[], float max_diff)
{
	int i;		/* counter */
	int ok = 0;	/* counter for ok */

	for(i=0;i<number_hkl;++i) {
		/* Calculate difference in hkl do decide on good or bad */
		if ( ( fabs(hkl[i].h-hkl_calc[i].h) < max_diff) &&
		    ( fabs(hkl[i].k-hkl_calc[i].k) < max_diff) &&
		    ( fabs(hkl[i].l-hkl_calc[i].l) < max_diff) ) {
			++ok;
		}
	}
	return(ok);
}

void bin_ok_hkl(struct hkl_fract hkl[], int number_hkl, 
struct hkl_int hkl_calc[], float max_diff)
{
	int i;		/* counter */
	int ok = 0;	/* counter for ok */
	int bin[3][20] = {};


	for(i=0;i<number_hkl;++i) {
	  int ih = (int)(fabs(hkl[i].h-hkl_calc[i].h) * 20);
	  int ik = (int)(fabs(hkl[i].k-hkl_calc[i].k) * 20);
	  int il = (int)(fabs(hkl[i].l-hkl_calc[i].l) * 20);
	  bin[0][ih]++;
	  bin[1][ik]++;
	  bin[2][il]++;
	}

	printf("HISTOGRAM OF DEVIATION FROM INTEGER\n");
	printf("          HD  KD  LD\n");
	for(i=0;i<=10;i++) {
	  printf("BIN %4.2f %3d %3d %3d\n", 0.05 * i, bin[0][i], bin[1][i], bin[2][i]);
	}
}
