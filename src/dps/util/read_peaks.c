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
 *                          File: read_peaks.c                         *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "util_def.h"
#include "read_peaks.h"


Boolean readinfo(FILE *in_file, detector_setting *det_set);

/* readinput() on a given filename this routine opens the inputfile, reads the info with readinfo()
		and closes the file. It returns the information read from the file and an
		exit value for successful reading.		*/
/*#define READ_RS_DEBUG */

Boolean readinput(char *filename, detector_setting *det_set)
{
	FILE *in_file;
	int number;
	Boolean status;

#ifdef READINPUT_DEBUG
	printf("READINPUT entered...\n");
#endif
	in_file = fopen(filename, "r");
	if (in_file == NULL) {
		printf("Can't open %s!! Exit/n", filename);
		return(FALSE);
	};
	status = readinfo(in_file, det_set);
	fclose(in_file);
#ifdef READINPUT_DEBUG
	printf("READINPUT left...\n");
#endif
	return(status);
}

/* readinfo() on a given filepointer the routine read the header 			*
 *			  information (see example ???, same format as bob) it returns	*
 *			  all the info read												*/
Boolean readinfo(FILE *in_file,  detector_setting *det_set)

{
	int i;
	char line[256];
	char command[100];
	char *rest;
	char name[100];
	char tmp_mode[20];
	float start, end;
	int line_count = 0;
	Boolean file_ok=TRUE;
	double pi = M_PI;
	Boolean command_identified;

#ifdef READ_INFO_DEBUG
	printf("READ_INFO DO_READ entered...\n");
#endif

	while (fgets(line, sizeof(line), in_file) != 0) {
		command_identified = FALSE;
		line_count++;
		if(sscanf(line, "%s", &command) <= 0) {
			command[0]='\0';
			command_identified = TRUE;
		}
		if(command[0] == '#') {
			command[0]='\0';
			command_identified = TRUE;
		}
		rest = &(line[strlen(command)]);
		upper(command);
#ifdef READ_INFO_DEBUG
		printf("line %d: %s\n", line_count, command);
#endif
		if (strcmp(command, "WAVELENGTH") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the input data file: No data after WAVELENGTH keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1: 
				det_set->lambda = start;
				break;
			default:
				printf("Error in line %d in the input data file: To many parameters after WAVELENGTH keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "RASTERSIZE") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the input data file: No data after RASTERSIZE keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1: 
				det_set->mm_per_raster = start;
				break;
			default:
				printf("Error in line %d in the input data file: To many parameters after RASTERSIZE keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "DISTANCE") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the input data file: No data after DISTANCE keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1: 
				det_set->dist_mm = start;
				break;
			default:
				printf("Error in line %d in the input data file: To many parameters after DISTANCE keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "BEAMSTOP") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start and end */
			switch(sscanf(rest, "%s %f %s %f", &command, &start,  &name, &end)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No data after BEAMSTOP keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1:
			case 2:
			case 3: 
				printf("Error in line %d in the control input file: Not enough or wrong parameters after BEAMSTOP keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 4: 
				if(strcmp(command, "R")) {
					if(strcmp(name, "S")) {
						det_set->beam_stop_r = start;
						det_set->beam_stop_s = end;
						break;
					}
					else {
						printf("Error in line %d in the control input file: Wrong parameters after BEAMSTOP  and R keyword\n",
						    line_count);
						file_ok=FALSE;
						break;
					}
				}
				if(strcmp(command, "S")) {
					if(strcmp(name, "R")) {
						det_set->beam_stop_r = end;
						det_set->beam_stop_s = start;
						break;
					}
					else {
						printf("Error in line %d in the input data file: Wrong parameters after BEAMSTOP  and S keyword\n",
						    line_count);
						file_ok=FALSE;
						break;
					}
				}
				printf("Error in line %d in the input data file: Wrong parameters after BEAMSTOP\n", line_count);
				file_ok=FALSE;
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after keyword BEAMSTOP\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "FILM_ROTATION") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the input data file: No data after FILM_ROTATION keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 1: 
				det_set->film_rot = start;
				start = det_set->film_rot/180.*pi;

				det_set->film_rot_mat[0][0]=cos(start);
				det_set->film_rot_mat[0][1]=-sin(start);
				det_set->film_rot_mat[1][0]=-(det_set->film_rot_mat[0][1]);
				det_set->film_rot_mat[1][1]=det_set->film_rot_mat[0][0];

				start=-start;

				det_set->film_rot_mat_inv[0][0]=cos(start);
				det_set->film_rot_mat_inv[0][1]=-sin(start);
				det_set->film_rot_mat_inv[1][0]=-(det_set->film_rot_mat_inv[0][1]);
				det_set->film_rot_mat_inv[1][1]=det_set->film_rot_mat_inv[0][0];
				break;
			default:
				printf("Error in line %d in the input data file: To many parameters after FILM_ROTATION keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "SSCALE") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the input data file: No data after SSCALE keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1: 
				det_set->y_scale = start;
				break;
			default:
				printf("Error in line %d in the input data file: To many parameters after SSCALE keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (command_identified == FALSE) {
			printf("Error in line %d in the control input file: Keyword %s unknown\n", line_count, command);
			file_ok=FALSE;
		}

	}
#ifdef READ_INFO_DEBUG
	printf("READ_INFO DO_READ left...\n");
#endif
	return(file_ok);
}


/* read_rs() on a given filepointer identify the type of the file and
		read all following lines with r and s*/

int read_rs(char *filename, struct rs_coord **rs)
{
	FILE *in_file;
	int i;			/* counter */
	char line[100];	/* inputline */
	float dummy;

#ifdef READ_RS_DEBUG
	printf("READ_RS entered...\n");
	printf("Try to open peaks file $s\n",filename);
#endif

	/* Try to open file */
	in_file = fopen(filename, "r");
	if (in_file == NULL) {
		printf("Can't open %s!!/n", filename);
		return(0);
	};

	/* Get the first line, identify the file and read it*/
	fgets(line, sizeof(line), in_file);
#ifdef READ_RS_DEBUG
	printf("Identification line [%s]\n",line);
#endif

	if(strncmp("DPS-PF A1.0", line, 10) == 0) {
#ifdef READ_RS_DEBUG
		printf("Identified DPS peaks file ASCII version 1.0\n");
#endif

		/* Get the number of lines existing in the file and allocate enough
	   memory to store them in the rs_coord array */
		i=0;
		while (fgets(line, sizeof(line), in_file) != 0) i++;
		*rs = realloc(*rs, i*sizeof(rs_coord));

		/* Rewind file, skip the first line and start reading */
		i=0;
		rewind(in_file);
		fgets(line, sizeof(line), in_file);

		while (fgets(line, sizeof(line), in_file) != 0) {
			sscanf(line,"%f %f", &((*rs)[i].r), &((*rs)[i].s));
#ifdef READ_RS_DEBUG
			printf("No: %d %8.3f %8.3f\n",i,(*rs)[i].r, (*rs)[i].s);
#endif
			i++;
		}
	}
	else if(strncmp("      7777", line, 10) == 0) {
#ifdef READ_RS_DEBUG
		printf("Identified DENZO peaks file\n");
#endif

		/* Get the number of lines existing in the file and allocate enough
	   memory to store them in the rs_coord array */
		i=0;
		while (fgets(line, sizeof(line), in_file) != 0) i++;
		*rs = realloc(*rs, i*sizeof(rs_coord));

		/* Rewind file, skip the first line and start reading */
		i=0;
		rewind(in_file);
		fgets(line, sizeof(line), in_file);

		while (fgets(line, sizeof(line), in_file) != 0) {
			sscanf(line,"%f %f %f", &dummy, &((*rs)[i].r), &((*rs)[i].s));
#ifdef READ_RS_DEBUG
			printf("No: %d %8.3f %8.3f\n",i,(*rs)[i].s, (*rs)[i].r);
#endif
			i++;
		}
	}
	else {
		printf("Could not identify peaks file.\n");
		return(0);
		fclose(in_file);

	}

	fclose(in_file);
#ifdef READ_RS_DEBUG
	printf("READ_RS left with %d reflections...\n",i);
#endif
	return(i);

}



/* readxyz() alternative to readpeaks() 										*
 *				on a given filename the routine returns the number of read 		*
 *				coordinates. The coordinates ar ein the the lists x, y and z.	*
 *				There is also a list with the length of the vector. 			*/

int readxyz(char *filename, struct vector *xyz)
{
	FILE *in_file;		/* Inputfile */
	char line[100];		/* Length of line to be read */
	int i;				/* Counter */
	float help1,help2;

#ifdef READXYZ_DEBUG
	FILE *help;
	help=fopen("xyz.read","w");
	printf("READXYZ entered...\n");
#endif
	in_file = fopen(filename, "r");
	if (in_file == NULL) {
		printf("Can't open %s/n", filename);
		exit(8);
	};
	i=0;
	while (i < N_XYZ) {
		if (fgets(line, sizeof(line), in_file) == 0) {
			break;
		};
		sscanf(line,"%f %f %f %f %f", &help1, &help2, &(xyz[i].x),&(xyz[i].y),&(xyz[i].z));
#ifdef READXYZ_DEBUG
		printf("x %f y %f z %f\n",xyz[i].x,xyz[i].y,xyz[i].z);
		fprintf(help,"%f %f %f\n",xyz[i].x,xyz[i].y,xyz[i].z);
#endif
		i++;
	}
#ifdef READXYZ_DEBUG
	fclose(help);
	printf("READXYZ left with %d refelctions read...\n",i-1);
#endif
	return(i-1);
}
