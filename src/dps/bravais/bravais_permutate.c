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
 *                              dps_bravais                            *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                       File: bravais_permutate.c                     *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "bravais.h"
#include "bravais_permutate.h"

/*
 *	The orignal test constant.
 *
 *	#define	EPSILON		(0.00000001)
 */
#define	EPSILON		(0.1)		/* the new test */

/* #define	VEC_OK		(0.10)   Done in bravais_permutate.h */

/* gen_permutations(): Generate new OM from the vectors of the given OM. 
 *		       Check if the basis vector are not more the VEC_OK
 *		       longer than the original. Return the number of 
 *		       generated OM's. (Kabsch, J. Appl. Cryst. 1993)
 */

int gen_permutations(struct matrix_3x3 *om, struct matrix_3x3 *oms)

{
	float len_max[3], len_min[3];	    /* Variables for storing the maximum allowed
						       vector length */
	float volume;			    /* Variable for storing the volume of the
						       original cell. */
	int new_count;	    /* Counter for the generated om's */
	int i0, i1, i2;	    /* Counter */
	int j0, j1, j2;	    /* Counter */
	int k0, k1, k2;	    /* Counter */
	float help;
	float xx;

	len_max[0]=vec_lenf(om->r[0])*(1+VEC_OK);
	len_max[1]=vec_lenf(om->r[1])*(1+VEC_OK);
	len_max[2]=vec_lenf(om->r[2])*(1+VEC_OK);
	len_min[0]=vec_lenf(om->r[0])*(1-VEC_OK);
	len_min[1]=vec_lenf(om->r[1])*(1-VEC_OK);
	len_min[2]=vec_lenf(om->r[2])*(1-VEC_OK);
	volume = vec_volumef(om->r[0], om->r[1], om->r[2]);
	/*    printf(" len_min: %f %f %f len_max: %f %f %f volume: %f\n",
		len_min[0],len_min[1],len_min[2],len_max[0],len_max[1],len_max[2],volume);*/
	new_count=0;
	for(i0=-1;i0<=1;i0++){
		for(i1=-1;i1<=1;i1++){
			for(i2=-1;i2<=1;i2++) {
				oms[new_count].r[0].x=i0*(om->r[0].x)+i1*(om->r[1].x)+i2*(om->r[2].x);
				oms[new_count].r[0].y=i0*(om->r[0].y)+i1*(om->r[1].y)+i2*(om->r[2].y);
				oms[new_count].r[0].z=i0*(om->r[0].z)+i1*(om->r[1].z)+i2*(om->r[2].z);
				help=vec_lenf(oms[new_count].r[0]);
				if((help < len_max[0]) && (help > len_min[0])) {
					for(j0=-1;j0<=1;j0++){
						for(j1=-1;j1<=1;j1++){
							for(j2=-1;j2<=1;j2++) {
								oms[new_count].r[1].x=j0*(om->r[0].x)+j1*(om->r[1].x)+j2*
								    (om->r[2].x);
								oms[new_count].r[1].y=j0*(om->r[0].y)+j1*(om->r[1].y)+j2*
								    (om->r[2].y);
								oms[new_count].r[1].z=j0*(om->r[0].z)+j1*(om->r[1].z)+j2*
								    (om->r[2].z);
								help=vec_lenf(oms[new_count].r[1]);
								if((help < len_max[1]) && (help > len_min[1])) {
/* code folded from here */
	for(k0=-1;k0<=1;k0++){
		for(k1=-1;k1<=1;k1++){
			for(k2=-1;k2<=1;k2++) {
				oms[new_count].r[2].x=k0*(om->r[0].x)+k1*(om->r[1].x)+k2*(om->r[2].x);
				oms[new_count].r[2].y=k0*(om->r[0].y)+k1*(om->r[1].y)+k2*(om->r[2].y);
				oms[new_count].r[2].z=k0*(om->r[0].z)+k1*(om->r[1].z)+k2*(om->r[2].z);
				help=vec_lenf(oms[new_count].r[2]);
				if((help < len_max[2]) && (help > len_min[2]) && (fabs(volume - vec_volumef(oms[new_count].r[0], oms[new_count].r[1], oms[new_count].r[2])) < EPSILON)) {
					new_count++;
					oms[new_count]=oms[new_count-1];
				}
			}
		}
	}
/* unfolding */
								}
							}
						}
					}
				}
			}
		}
	}
	return(new_count);
}
