/* locally modified version of bravais.c to return the 44 orientation matrices
   (i.e. one for each characteristic lattice), the calculated FOM, and the 
   Bravais Lattice rather than write out to a file. 
   This is a bit of a pain as I don't really want to start having to hack 
   Ingo's code and creating version problems, but that's life.
   
                                                            Harry, August 1998
*/

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
 *                            File: bravais.c                          *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/util.h"
#include "bravais.h"

void bravais(struct matrix_3x3 *om, struct matrix_3x3 cl_om[45], 
float crit[45],char sys_id[45][3], int output_type)

{
	int i, j;				/* Counter */
	struct matrix_3x3 dir_om;		/* Direct cell form om */
	struct matrix_3x3 perm_om[1000];	/* Permutated OM's */
	int num_om;				/* Number of permutated OM's */
	int selected;			/* Number of selected om */
	float min;				/* Minimal value found for criterium */
	float help;
	struct cell_struct cell;
	struct matrix_3x3 trans_om, trans_om_rec; /* Transformed om direct and rec*/

	/* Mapping the system string to the lattice character (1-44) */
	char system[45][3] = {
		'0','0','\0',
		'c','F','\0',
		'h','R','\0',
		'c','P','\0',
		'h','R','\0',
		'c','I','\0',
		't','I','\0',
		't','I','\0',
		'o','I','\0',
		'h','R','\0',
		'm','C','\0',
		't','P','\0',
		'h','P','\0',
		'o','C','\0',
		'm','C','\0',
		't','I','\0',
		'o','F','\0',
		'm','C','\0',
		't','I','\0',
		'o','I','\0',
		'm','C','\0',
		't','P','\0',
		'h','P','\0',
		'o','C','\0',
		'h','R','\0',
		'm','C','\0',
		'o','F','\0',
		'm','C','\0',
		'm','C','\0',
		'm','C','\0',
		'm','C','\0',
		'a','P','\0',
		'o','P','\0',
		'm','P','\0',
		'm','P','\0',
		'm','P','\0',
		'o','C','\0',
		'm','C','\0',
		'o','C','\0',
		'm','C','\0',
		'o','C','\0',
		'm','C','\0',
		'o','I','\0',
		'm','I','\0',
		'a','P','\0'	};
	for (i=0;i<=44;i++)
	{
		for (j=0;j<=2;j++)
		{
			sys_id[i][j] = system[i][j];
		}
	}
	/* Also the cell should have been already reduced, we reduce it again */
	/* print_3x3(*om);*/
	reduce_cell(om);

	/* Output the input values
    if (output_type == 0) {
      printf("Original Orientation Matrix:\n");
      print_3x3(*om);
    }
     */
	/* Make sure the order of the cell vectors is correct (short to long) */

	om_from_vectors(om->r[0], om->r[1], om->r[2], &dir_om);
	reccell_from_om(dir_om, &cell);
	if (cell.a > cell.b || cell.b > cell.c) {
		resort_cell_vectors(om,&cell);
		om_from_vectors(om->r[0], om->r[1], om->r[2], &dir_om);
		reccell_from_om(dir_om, &cell);
		if (cell.a > cell.b || cell.b > cell.c)
			printf ("Warning, can't get cell conforming with convention!\n");
	}
	if (output_type == 0)
		print_cell(cell);
	/* generate all possible permutations of the vectors */
	num_om = gen_permutations(&dir_om, &(perm_om[0]));
	if (output_type == 0) {
		printf(" number of om generated: %d\n",num_om);
		printf("Final Results:\n");
	}

	/* generate all possible permutations of the vectors */
	num_om = gen_permutations(&dir_om, &(perm_om[0]));
	/*
    printf(" number of om generated: %d\n",num_om);
    printf("Final Results:\n");
    */
	/* Loop through the 44 lattice characters. */
	for(i=1;i<=44;i++) {
		min=100000000000.0;

		/* For each character get the permutation with the smallest criterion 
	   and print this cell with the criterion */

		for(j=0;j<num_om;j++) {

			/* calculate the criterion */
			help=get_goof(i, &(perm_om[j]));
			if(help < min) {
				min=help;
				selected=j;
			}
		}

		/* Print the output for this character */
		if (output_type == 0)
			printf("Lattice: %2d System: %s crit: %10.1f \n", i,
			    system[i], min);
		transform_om(i, &perm_om[selected], &trans_om);
		om_from_vectors((trans_om.r[0]), (trans_om.r[1]),
		    (trans_om.r[2]),
		    &trans_om_rec);
		if (output_type == 0) {
			print_3x3(trans_om_rec);
			print_cell_from_om(trans_om_rec);
		}
		else {
			for (j=0;j<=2;j++) {
				sys_id[i][j] = system[i][j];
			}
			cl_om[i] = trans_om_rec;
			crit[i] = min;
		}

	}
}
