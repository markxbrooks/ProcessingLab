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
 *                              File: index.c                          *
 *                                                                     *
 *=====================================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "index.h"
#include "ind_gen_direct.h"
/* do_index(): This routine is the main routine for indexing. On a given*
 *	    list of xyz-coordinates and the maximum expected length	*
 *	    of the cell edges,  it tries to index the reflections. The	*
 *	    output is a 3 by 3 matrix with the orientation matrix.	*
 *	    At the moment is it only possible to run index in automatic *
 *	    mode. In this mode index will select the orientation matrix *
 *	    by itself.	Later other modes with certain influence of the *
 *	    will be included.						*/

int do_index(int numvect,
             int numvectref,
             int maxcyc,
             int inicyc,
             float convrg,
             int numormat,
             float stepscal,
             int numrlp,
             int numsptdps,
             int npkfftmin,
             int nmultfft,
             int idumpth,
             float theta1d,
             float theta2d,
             float sinthmin,
             float volmin,
             float hkldev,
             float hklmin,
             int dumpfft[],
             int dumpvec,
             int dumpref,
             int dumpgen,
             int dumpch,
             int number_xyz, 
	     struct vector xyz[], 
	     float d_max, 
	     float d_min, 
	     float tarcell[],
	     Index_Mode index_mode, 
	     struct matrix_3x3 *om_out)

{
  int i,j,dpserrflg,tarcellflg;				/* Counter */
  int number_dir;			/* Number of direct lattice vectors (t_vectors)*/
  int number_cut;
  struct tvector t_vector[N_TVECTOR];	/* Array for the t-vectors. */
  int num_om = 0;			/* Number of generated orientation matrices (num_om < NUM_OM_MAX) */
	int selected_om;			/* Finally selected orientation matrix */
	struct cell_struct cell, cell_rec;	/* Direct and reciprocal cell generated from om */
	struct matrix_3x3 om[MAX_OM_COUNT];	/* Orientation matrix array*/
	float cut_crit = CUT_CRIT;		/* Criterium for cutting t-vector list */
	if(tarcell[0] == 0.0){
	  tarcellflg = 0;
	    } else {
	  tarcellflg = 1;
	}

        /* printf("values of dumpfft\t%d %d %d \n",dumpfft[0],dumpfft[1],dumpfft[2]); */

	/* First generate a list of possible direct space vectors. This is
       the main part of the indexing program				*/
	/* for debugging...*/

	number_dir = gen_direct_lattice_vectors(numvect, numvectref, 
                     maxcyc, inicyc, convrg, numormat,  stepscal,
                     numrlp, numsptdps, npkfftmin, nmultfft,
		     idumpth, theta1d, theta2d,
		      dumpfft, dumpvec, dumpref, dumpgen,
                     d_max, d_min, xyz, &number_xyz, &t_vector[0]);

	/* Debug output of vector list */
        if (dumpvec != 0) {
          printf("t_vectors after gen_direct_lattice_vectors\n");
	  print_tvector_list(number_dir,&t_vector[0]);
	}

	sort_tvectors_on_n_and_l(number_dir, &t_vector[0]);

	/* Debug output of sorted vector list */
        if (dumpvec != 0) {
	  printf("t_vectors after sort on n and l\n");
	  print_tvector_list(number_dir,&t_vector[0]);
	}

	/* Merge t_vectors with same length, same or opposite direction and
       similar number of fitting reflections. For more detail look at the
       comments in ind_sort_merge.c					*/

	number_dir = merge_tvectors_same_n(number_dir, &t_vector[0]);


	/* Debug output of merged vector list */
        if (dumpvec != 0) {
	  printf("t_vectors after merge\n");
	  print_tvector_list(number_dir,&t_vector[0]);
	}


	/* Now finally sort the t_vector list in respect to the length of 
       the vectors (Shortest first!). This will generate the smallest unit
       cell as the first generated orientation matrix.			*/

       sort_tvectors_on_l(number_dir, &t_vector[0]);


       /* Debug output of vector list after sort */
       if (dumpvec != 0) {
	  printf("t_vectors after sort on length\n");
	  print_tvector_list(number_dir,&t_vector[0]);
	  printf("\n\n");
	  printf("Angle between vectors\n");
	  printf("    ");
	  for (i = 0; i < number_dir; i++) printf("%3d ", i);
	  printf("\n");
	  for (i = 0; i <= number_dir; i++) {
	    printf("%3d ", i);
	    for (j = 0; j < number_dir; j++) {
	      if (j >= i) printf("  | ");
	      else printf("%3.0f ", vec_angle(t_vector[i].vec, t_vector[j].vec) / M_PI * 180.0);
	    }
	    printf("\n");
	  }
	  printf("\n\n");

	}


	/* Choose non linear combinations of direct space vectors to generate unit
	   cells, make the reciprocal unit cell (orientation matrix) out of this
	   and store them in an array. Generate MAX_OM_COUNT or as much as possible
	   orientation matrices.    */

	num_om = choose_non_lin(dumpch, sinthmin, volmin, hkldev, hklmin,
				number_xyz, xyz, number_dir, &t_vector[0], MAX_OM_COUNT, tarcell, &om[0]);

	dpserrflg = 0;
	if (num_om < 1) {
	    /* No matrix found, exit */
	    printf("DIRECT SPACE VECTORS DID NOT YIELD AN ORIENTATION MATRIX\n\n");
            printf("Try indexing with a smaller or larger maximum cell edge or try REFIX!\n");
            printf("A useful idea is to use the MEASURE CELL option before indexing\n");
	    dpserrflg = 1;
	    return (dpserrflg);
       	}


	/* Finally, select the 'best' of these as the result. More on how the
       best is chosen in ind_sort_merge.c				*/

	selected_om = choose_best_om(dumpch, hkldev, num_om,  om,  number_xyz,  xyz,tarcellflg);
	if (dumpch != 0) printf("selected %d\n",selected_om);

	/* Sort, copy and return the final orientation matrix to the calling routine  */

	*om_out = om[selected_om];
	/* hrp 01.09.2004
	   I don't think this is necessary as the reduce_cell stuff will do
	   all the permutations anyway, and test against type I and type II 
	   cells
	   sort_cell_vectors(om_out);
	*/
	reduce_cell(om_out);
	return (dpserrflg);
}


