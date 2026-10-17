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
 *                          File: dps_index.c                          *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "index.h"
#include "dps_index_read.h"

main(int argc, char *argv[])

{
	int i,called_from_mosflm;	    /* Counter */
	char *input_filename, 	    /* Filename for Control Input, */
	*detector_filename,	    /* Detector geometry info (Input) */
	*om_filename;	    	    /* Output data...OM */
	index_control control;	    /* Control parameter for the index run */
	detector_setting det_set;	    /* Detector setting  (see util_def.h) */
	int number_rs;		    /* Number of reflections in Inputfile */
	struct rs_coord *rs;	    /* Array for file coordinates */
	int number_xyz;		    /* Number of xyz points in resmin/resmiax */
	struct vector xyz_all[N_RS];    /* Array for reciprocal space coordinates */
	struct vector xyz[N_XYZ];	    /* Array for reciprocal space coordinates 
					       in a selected resolution range*/
	struct matrix_3x3 om;	    /* Generated orientation matrices */
	/* Check arguments */

	/* No arguments, print help message and exit */
	if (argc > 1) {
		printf("agc %d\n",argc);
		print_help();
		exit(1);
	}

	/* Read Control Input */
	control.om_filename = NULL;
	if (!read_control_do_read(stdin, &control, &det_set)) {
		printf("Error Reading Control Input File from stdin.\n");
		exit(1);
	}


	/* Read Reflections from file and convert all to reciprocal space 
       coordinates. Than restrict the list to a given resolution range.	    */

	/* Read reflections */
	rs = NULL;
	number_rs = read_multiple_rs(&control, &rs);

	/* Calculate reciprocal space coordinates */
	rs_2_xyz(number_rs, det_set, rs, &xyz_all[0],called_from_mosflm);

	/* Restrict resolution */
	number_xyz = cut_xyz(number_rs, xyz_all, control.resmin, control.resmax, &xyz[0]);


	/* Alternatively read reciprocal space coordinates directly */
#ifdef UNDEF
	/* Read reciprocal space coordinates */
	number_xyz=readxyz(peaks_filename,&xyz[0]);

	/* Should be deleted!!! */
	number_rs=number_xyz;

#endif

#ifdef DPS_INDEX_DEBUG
	print_vector_list(number_xyz, xyz);
	printf("Reflections in input: %d xyz coordinates in resmin/resmax: %d\n", number_rs,number_xyz);
#endif

	/* Ok, now let's index this stuff... */
	do_index(number_xyz, xyz, control.d_max, control.mode, &om);

	/* Show it to us.... */
	printf("Final Results:\n");
	print_3x3(om);

	/* Generate cell parameters and print them out */
	print_cell_from_om(om);

	fflush(stdout);
	if(control.om_filename != NULL) {
		write_om(control.om_filename, om);
	}

}
