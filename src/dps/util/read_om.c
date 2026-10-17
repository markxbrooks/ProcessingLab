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
 *                            File: read_om.c                          *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util_def.h"
#include "read_om.h"

/* read_om(): reads the orientation matrix from a file */

Boolean read_om(char *filename, struct matrix_3x3 *om)

{
	FILE *input;
	char line[100];

	input = fopen(filename, "r");
	if (input == NULL) {
		printf("Can't open %s!! Exit/n", filename);
		return(FALSE);
	};

	fgets(line, sizeof(line), input);
	sscanf(line,"%f %f %f\n",&(om->r[0].x),&(om->r[1].x),&(om->r[2].x));
	fgets(line, sizeof(line), input);
	sscanf(line,"%f %f %f\n",&(om->r[0].y),&(om->r[1].y),&(om->r[2].y));
	fgets(line, sizeof(line), input);
	sscanf(line,"%f %f %f\n",&(om->r[0].z),&(om->r[1].z),&(om->r[2].z));
	fclose(input);
	return(TRUE);
}

/* write_om(): writes the orientation matrix to a file */

Boolean write_om(char *filename, struct matrix_3x3 om)

{
	FILE *output;
	char line[100];

	output = fopen(filename, "w");
	if (output == NULL) {
		printf("Can't open %s!! Exit/n", filename);
		return(FALSE);
	};

	fprintf(output,"% 10.6f % 10.6f % 10.6f\n",om.r[0].x,om.r[1].x,om.r[2].x);
	fprintf(output,"% 10.6f % 10.6f % 10.6f\n",om.r[0].y,om.r[1].y,om.r[2].y);
	fprintf(output,"% 10.6f % 10.6f % 10.6f\n",om.r[0].z,om.r[1].z,om.r[2].z);
	fclose(output);
	return(TRUE);
}

