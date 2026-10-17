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
 *                              dps_index                              *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                        File: ind_fortran.c                          *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "index_def.h"
#include "ind_gen_direct.h"

/* index_fortran_(): This routine is the call for the indexing from a	*
 *		    FORTRAN program.					*/

void index_fortran_( int *numvect, int *numvectref, int *maxcyc, int *inicyc,
                     float *convrg, int *numormat, float *stepscal,
                     int *numrlp, int *numsptdps, int *npkfftmin, int *nmultfft,
                     int *idumpth, float *theta1d, float *theta2d,
                     float *sinthmin, float *volmin, float *hkldev, float *hklmin,
                     int dumpfft[], int *dumpvec, int *dumpref, int *dumpgen,
                     int *dumpch,
                     int *number_xyz, int*nlattice, float *x,  float *y, float *z, 
		     int *usespot, float *d_max,
                     float *d_min, int *index_mode, float *a, float *wave, int *ierrflg, float *tarcell, int *centering)

{
  int i,icount,dpserrflg;                          /* Counter */
  int numvectmax = *numvect + (int)sqrt(*d_max);
  int numvectrefmax = *numvectref + (int)sqrt(*d_max);
	struct vector xyz[N_XYZ];       /* Array for reciprocal space coordinates 
	                                       in a selected resolution range*/
	struct matrix_3x3 om;           /* Generated orientation matrices */
	float origcell[6];
	for (i = 0; i < 6; i++) origcell[i] = tarcell[i];
	if (*centering < 0) *centering = 0;
	if (*centering > 6) *centering = 0;

	/* Make a list of type vector from the coordinates x,y,z. */
	icount = 0;
	for(i=0;i<(*number_xyz);i++) {
	  if(usespot[i]==*nlattice)
	    {
	      /*	      if(icount<10)
		printf("%d %d, xyz[%d] = %f %f %f\n",usespot[i],*nlattice,i,x[i],y[i],z[i]);
	      */
	      xyz[icount].x=x[i];
	      xyz[icount].y=y[i];
	      xyz[icount++].z=z[i];
	    }
	}
	if (tarcell[0] != 0) {
	  printf("Target cell: %f %f %f %f %f %f ", 
		 tarcell[0], tarcell[1], tarcell[2],
		 tarcell[3], tarcell[4], tarcell[5]);
	  if (*centering != 0) printf("with %c centering", "PABCFIH"[*centering]);
	  printf("\n");

	  struct matrix_3x3 orig_cell, new_cell;
	  abc_to_cell(tarcell, &orig_cell);
	  new_cell = orig_cell;

	  if (*centering == 1) { // A centering
	    // C = (b + c) / 2
	    new_cell.r[2].x = (orig_cell.r[1].x + orig_cell.r[2].x) / 2;
	    new_cell.r[2].y = (orig_cell.r[1].y + orig_cell.r[2].y) / 2;
	    new_cell.r[2].z = (orig_cell.r[1].z + orig_cell.r[2].z) / 2;
	  } else if (*centering == 2) { // B centering
	    // A = (a + c) / 2
	    new_cell.r[0].x = (orig_cell.r[0].x + orig_cell.r[2].x) / 2;
	    new_cell.r[0].y = (orig_cell.r[0].y + orig_cell.r[2].y) / 2;
	    new_cell.r[0].z = (orig_cell.r[0].z + orig_cell.r[2].z) / 2;
	  } else if (*centering == 3) { // C centering
	    // B = (a + b) / 2
	    new_cell.r[1].x = (orig_cell.r[0].x + orig_cell.r[1].x) / 2;
	    new_cell.r[1].y = (orig_cell.r[0].y + orig_cell.r[1].y) / 2;
	    new_cell.r[1].z = (orig_cell.r[0].z + orig_cell.r[1].z) / 2;
	  } else if (*centering == 4) { // F centering
	    // B = (a + b) / 2
	    new_cell.r[1].x = (orig_cell.r[0].x + orig_cell.r[1].x) / 2;
	    new_cell.r[1].y = (orig_cell.r[0].y + orig_cell.r[1].y) / 2;
	    new_cell.r[1].z = (orig_cell.r[0].z + orig_cell.r[1].z) / 2;
	    // C = (a + c) / 2
	    new_cell.r[2].x = (orig_cell.r[0].x + orig_cell.r[2].x) / 2;
	    new_cell.r[2].y = (orig_cell.r[0].y + orig_cell.r[2].y) / 2;
	    new_cell.r[2].z = (orig_cell.r[0].z + orig_cell.r[2].z) / 2;
	  } else if (*centering == 5) { // I centering
	    // C = (a + b + c) / 2
	    new_cell.r[2].x = (orig_cell.r[0].x + orig_cell.r[1].x + orig_cell.r[2].x) / 2;
	    new_cell.r[2].y = (orig_cell.r[0].y + orig_cell.r[1].y + orig_cell.r[2].y) / 2;
	    new_cell.r[2].z = (orig_cell.r[0].z + orig_cell.r[1].z + orig_cell.r[2].z) / 2;
	  } else if (*centering == 6) { // H centering
	    // A = (2a + b + c) / 3
	    new_cell.r[0].x = (2 * orig_cell.r[0].x + orig_cell.r[1].x + orig_cell.r[2].x) / 3;
	    new_cell.r[0].y = (2 * orig_cell.r[0].y + orig_cell.r[1].y + orig_cell.r[2].y) / 3;
	    new_cell.r[0].z = (2 * orig_cell.r[0].z + orig_cell.r[1].z + orig_cell.r[2].z) / 3;
	    // B = (-a + b + c) / 3
	    new_cell.r[1].x = (-orig_cell.r[0].x + orig_cell.r[1].x + orig_cell.r[2].x) / 3;
	    new_cell.r[1].y = (-orig_cell.r[0].y + orig_cell.r[1].y + orig_cell.r[2].y) / 3;
	    new_cell.r[1].z = (-orig_cell.r[0].z + orig_cell.r[1].z + orig_cell.r[2].z) / 3;
	    // C = (-a - 2b + c) / 3
	    new_cell.r[2].x = (-orig_cell.r[0].x - 2 * orig_cell.r[1].x + orig_cell.r[2].x) / 3;
	    new_cell.r[2].y = (-orig_cell.r[0].y - 2 * orig_cell.r[1].y + orig_cell.r[2].y) / 3;
	    new_cell.r[2].z = (-orig_cell.r[0].z - 2 * orig_cell.r[1].z + orig_cell.r[2].z) / 3;
	  }
	  if (*centering != 0) {
	    cell_to_abc(&new_cell, tarcell);
	    printf("%f %f %f %f %f %f in primitive setting.\n", 
		   tarcell[0], tarcell[1], tarcell[2], tarcell[3], tarcell[4], tarcell[5]);	    
	  }
	}
	printf("Trying to index lattice %i with %i reflections\n",
	       *nlattice,icount);


	if (*ierrflg > 0)
	  print_vector_list( icount, xyz);
	dpserrflg = 1;
	while ( dpserrflg != 0 && *numvect <= numvectmax && *numvectref <= numvectrefmax ) {
	  dpserrflg = do_index(*numvect, *numvectref, *maxcyc, *inicyc,
			       *convrg, *numormat,  *stepscal,
			       *numrlp, *numsptdps, *npkfftmin, *nmultfft,
			       *idumpth, *theta1d, *theta2d,
			       *sinthmin,  *volmin,  *hkldev,  *hklmin,
			       dumpfft, *dumpvec, *dumpref, *dumpgen, *dumpch,
			       icount, xyz, *d_max, *d_min, tarcell, *index_mode, &om);
	  if (origcell[0] != 0 ) break;
	  if ( dpserrflg != 0 ) {
	    *numvect += 10;
	    *numvectref += 10;
	  }
	}
	*ierrflg = dpserrflg;
	
	/*
	  HRP 29.03.2010 change the matrix from DPS to Mosflm here to save confusion
	  later 
	*/
	a[0] = *wave * om.r[0].z;
	a[1] = *wave * om.r[0].x;
	a[2] = *wave * om.r[0].y;
	a[3] = *wave * om.r[1].z;
	a[4] = *wave * om.r[1].x;
	a[5] = *wave * om.r[1].y;
	a[6] = *wave * om.r[2].z;
	a[7] = *wave * om.r[2].x;
	a[8] = *wave * om.r[2].y;

	// restore the tarcell destroyed in ind_choose.c
	for (i = 0; i < 6; i++) tarcell[i] = origcell[i];
}


