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
 *                         File: dps_frameinfo.c                       *
 *                                                                     *
 *=====================================================================*/

/* dps_frameinfo is a program to identify the type of a image file.
   It is primarily a program to check the read_frame() function and
   should serve as an example how to use this function. It can also
   serve to check, if a certain image file format can be detected by
   the DPS programs.				

   Input: dps_frameinfo frame.file					*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"

main(int argc, char *argv[])

{
	FrameInfo frameinfo;
	od_type *frame, *frame2, *frame3, *frame4;
	int i;
	unsigned long int sum = 0;
	unsigned char *test;
	int j;

	if(argc == 1) {
		printf("Usage:  dps_frameinfo frame.file\n\n");
		exit(1);
	}

	frame = NULL;
	if(read_frame(argv[1], &frameinfo, &frame) == 0) {

		printf("maxod: %d minod %d\n",frameinfo.maxod, frameinfo.minod);

		/* Print out the first 200 values of the file */
		j=0;
		printf("0x%2x ",0);
		for(i=0;i<200;i++) {
			j++;
			printf("%5d ",frame[i]);
			if(j==16) {
				j=0;
				printf("\n0x%2x ",i+1);
			}
		}
		printf("\n");
	}
}
