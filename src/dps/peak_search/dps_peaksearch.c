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
 *                         File: dps_peaksearch.c                      *
 *                                                                     *
 *=====================================================================*/

/* dps_peaksearch

   Input: dps_peakssearch frame.file					*/

/* This program does a peak search on a given image and returns a list of 
   R, S coordinates in pixel. It uses the read_frame routine from the util
   library and a modified algorithm of Sangsoo Kim.			*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"

main(int argc, char *argv[])

{

	int stepx = 3;	/* Initial stepsize for scanning through the image */
	int stepy = 3;	/* Initial stepsize for scanning through the image */
	float back_count = 2*(2*stepx+1)+4*stepy-2; /* Number of background pixels */
	float spot_count = (2*stepx-1)*(2*stepy-1); /* Number of spot pixels (see below) */
	FrameInfo frameinfo; /* frameinfo structure (return value of read_frame) */
	od_type *frame; /* The frame... */
	int i, j; /* counter */
	int y, dx; /* more counter */
	od_type value; /* variable to store an actual od_value */
	int k, l; /* many more counter... */
	int x_max, y_max;
	od_type maxval;
	int bma_x, bmi_x, bma_y, bmi_y;
	int back, spot;
	FILE *out;
	static float noise_thresh = 1.2;
	double A, B, I, sigmaI;

	if(argc < 3) {
		printf("Usage:  dps_peaksearch frame.file peaks.file\n\n");
		exit(1);
	}

	/* Open the output file */
	out = fopen(argv[2],"w");
	if (out == NULL) {
		printf("Could not open output file [%s]!\n", argv[2]);
	}
	if(argc == 4)
	{
		noise_thresh = atof(argv[3]);
		fprintf(stdout,"dps_peak_search: noise_thresh set to %8.2f\n",noise_thresh);
	}
	frame = NULL;

	/* Read the frame */
	if(0 == read_frame(argv[1], &frameinfo, &frame)) {

		/* Write header line to peaks file */
		fprintf(out,"DPS-PF A1.0\n");

		/* The next two loops go over the whole frame with stepsize 
	     * stepx and stepy. */
		for(i=2*stepy;i<frameinfo.leny-2*stepy;i=i+stepy) {
			for(j=2*stepx;j<frameinfo.lenx-2*stepx;j=j+stepx) {

				/* y hold the index of the pixel at i,j */
				y = i*frameinfo.lenx+j;

				/* dx is the difference in index between i,j and
			 * i, j+stepx */
				dx = stepy*frameinfo.lenx;

				value=frame[y];

				/* Check if we have a maximum at i,j */
				if ((value > frame[y+stepy]) && 	
				    (value > frame[y-stepy]) &&
				    (value > frame[y+dx])    &&
				    (value > frame[y-dx])    &&
				    (value > frame[y+stepy+dx]) && 	
				    (value > frame[y-stepy+dx]) &&
				    (value > frame[y+stepy-dx]) &&
				    (value > frame[y-stepy-dx])        ) {

					/* If we have a maximum, try to find the maximum
			     * in a box around i, j with size 2stepy * 2stepx
			     * and stepsize 1. */
					maxval=frame[y];

					/* x_max, y_max will hold the final maximum */
					x_max=j;
					y_max=i;
					for(k=i-stepy;k<=i+stepy;k++) {
						for(l=j-stepx;l<=j+stepx;l++) {

							/* Same as above only with stepsize 1 */
							y = k*frameinfo.lenx+l;
							dx = frameinfo.lenx;
							value=frame[y];
							if((value >= frame[y+1])  && 	
							    (value >= frame[y-1])  &&
							    (value >= frame[y+dx]) &&
							    (value >= frame[y-dx]) &&
							    (value >= maxval)         ) {
								maxval=value;
								x_max=l;
								y_max=k;
							}
						}
					}
					/* Now we calculate the average background and
			     * spot values for this position. The box goes
			     * from i-stepy to i+stepy and j-stepx to
			     * i+stepx. The backgroiund pixels are the pixels
			     * of a one pixel frame at the border of the box, 
			     * the rest are spot pixels. For stepy=stepy=3 this
			     * gives 25 spot pixels and 24 background pixels */

					/* Borders of the box */
					bma_y = y_max+stepy;
					bmi_y = y_max-stepy;
					bma_x = x_max+stepx;
					bmi_x = x_max-stepx;
					back = 0;
					spot = 0;
					for(k=bmi_y;k<=bma_y;k++) {
						for(l=bmi_x;l<=bma_x;l++) {

							/* see above */
							y = k*frameinfo.lenx+l;
							dx = frameinfo.lenx;
							value=frame[y];

							/* if counter at border of box, the pixel
				     * is a background pixel otherwise it is a
				     * spot pixel */
							if ((k == bma_y) ||
							    (k == bmi_y) ||
							    (l == bma_x) ||
							    (l == bmi_x)    ) {
								back = back + value;
							}
							else {
								spot = spot + value;
							}
						}
					}
					/* If the average spot pixel value is larger by a 
			     * certain factor than the average background pixel value
			     * write the reflection to the output file. Check also
			     * if the maximum value is an overload value.... */

					/*
			    if (((float)spot/spot_count > noise_thresh * (float)back/back_count) &&
			        (maxval != frameinfo.overload) ) { 
			     */

					A = (double)spot;
					B = (double)back*spot_count/back_count;
					I = A - B;
					sigmaI = sqrt(A + B);
					if ( (I/sigmaI > noise_thresh) && (maxval != frameinfo.overload) ) {

						if (!near_edge(frame, frameinfo.lenx, frameinfo.lenx, x_max, y_max))
							if (!no_use(x_max, y_max))
								fprintf(out,"%d %d \n", y_max,  x_max);
					}
				}
			}
		}
	}
	fclose(out);
}

#define MINDIST	(5)

/* Test if there is a "0" pixel within +/-
 * MINDIST of the spot.
 */
near_edge(data, width, height, xpos, ypos)
unsigned short *data;
int width, height;
int xpos, ypos;
{
	int i,j,x,y;

	for(i = -MINDIST; i<= MINDIST; i++) {
		for(j = -MINDIST; j<= MINDIST; j++) {
			x = xpos + j;
			y = ypos + i;
			if ((x >= 0) && (x <= width) &&
			    (y >= 0) && (y <= height)) {
				if (data[width * y + x] == 0) {
					return 1;
				}
			}
			else {
				return 1;
			}

		}
	}
	return 0;
}

/* Special test for SN406. 11/15/97
 */
no_use(x, y)
int x, y;
{
	/*	if ((x >= 1860) && (x <= 2304) && (y >= 1196) && (y <= 1711))
		return 1; */

	return 0;
}
