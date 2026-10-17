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
 *                         File: dps_bravais.c                         *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "bravais.h"
#include "dps_bravais_read.h"

main(int argc, char *argv[])

{
	int i;			    /* Counter */
	char *om_filename;	    	    /* Output data...OM */
	struct matrix_3x3 om;	    /* Read OM */
	struct matrix_3x3 om_compare;   /* OM for comparison reasons */
	struct matrix_3x3 cl_om[45];
	float crit[45];
	char sys_id[45][3];
	int output_type;
	BravaisMode mode;		    /* Mode for Bravais */
	bravais_control control;	    /* Control for BRavais */

	/* Check arguments */

	/* No arguments, print help message and exit */
	if (argc > 1) {
		printf("agc %d\n",argc);
		print_help();
		exit(1);
	}

	/* Read Control Input */
	if (!read_control_do_read(stdin, &control)) {
		printf("Error Reading Control Input File.\n");
		exit(1);
	}

	/* Check if om filname was given */
	if (control.om_filename == NULL) {
		printf("No orientation matrix ouput file specified...exit.\n");
		print_help();
		exit(1);
	}

	/* Read om from file */
	read_om(control.om_filename, &om);

	if(control.mode == BRAVAIS_COMPARE) {
		/*	read_om((control.om_filename), &om_compare);
	bravais_compare(&om_compare, &om);
*/
		printf("This feature has not been implmented....");
	}
	else {
		/* Run bravais */
		bravais(&om,cl_om,crit,sys_id,output_type);
	}




}
