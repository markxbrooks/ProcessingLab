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
 *                           File: ind_help.c                          *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <math.h>
#include "util.h"
#include "index.h"

/* print_help(): Print command line options			 	*/

void print_help(void)
{
	printf("Usage: index < control_input \n");

}

/* print_tvector_list(): Print the list of generated direct space vectors. */

void print_tvector_list(int num,  struct tvector t[])

{
	int i;		/* Counter */

	printf("List of generated direct space vectors:\n");
	printf("No.    theta1   theta2    x        y        z       len       fit\n");
	for(i=0;i<=num;i++) {
		printf("%3d: %8.4f %8.4f % 8.3f % 8.3f % 8.3f %8.2f %7.2f\n", i, t[i].theta1,t[i].theta2,t[i].vec.x, t[i].vec.y, t[i].vec.z, t[i].len, t[i].fit);
	}
}

