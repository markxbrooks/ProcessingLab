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
 *                        File: dps_index_read.c                       *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#if !defined vms && !defined(__APPLE__) && ! ( defined (__FreeBSD__) && defined (__STDC__) )
#include <malloc.h>
#endif
#include "util.h"
#include "index_def.h"
#include "dps_index_read.h"

/* read_control(): Read the Control input file and analyse it. Return the
    		  control structure and an exit value for successfull
		  reading.						*/

Boolean read_control_do_read(FILE *in_file,  index_control *control, detector_setting *det_set)

{
	int i;
	char line[256];
	char command[100];
	char *rest;
	peakfile **files;
	int peakfile_counter = 0;
	char name[100];
	char tmp_mode[20];
	float start, end;
	int line_count = 0;
	Boolean file_ok=TRUE;
	double pi = M_PI;
	Boolean command_identified;
	Boolean film_rot;
	Boolean beam_stop;

#define MAX_PEAKFILES 2000

#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL DO_READ entered...\n");
#endif

	/* Initialize values */
	control->d_max = -1.0;
	control->om_filename = NULL;
	control->resmax = -1.0;
	control->resmin = -1.0;
	control->mode=-1.0;
	det_set->lambda = -1.0;
	det_set->mm_per_raster = -1.0;
	det_set->dist_mm = -1.0;
	beam_stop=FALSE;
	film_rot = FALSE;
	det_set->y_scale = 0.0;

	files = calloc(MAX_PEAKFILES,  sizeof(peakfile *));
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
#ifdef READ_CONTROL_DEBUG
		printf("line %d: %s\n", line_count, command);
#endif
		if (strcmp(command, "CELL") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%f ", &start)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No data after CELL keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1: 
				control->d_max = start;
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after CELL keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "OMFILE") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start */
			switch(sscanf(rest, "%s ", &command)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No filename after OMFILE keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 1: 
				control->om_filename = makestring(command);
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after OMFILE keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "RESOLUTION") == 0) {
			command_identified = TRUE;
			/* Sorry for reusing start and end */
			switch(sscanf(rest, "%f %f", &start,  &end)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No data after RESOLUTION keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 1: 
				printf("Error in line %d in the control input file: Not enough or wrong parameters after RESOLUTION keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 2: 
				if (start < end) {
					control->resmax = start;
					control->resmin = end;
				}
				else {
					control->resmin = start;
					control->resmax = end;
				}
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after RESOLUTION keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
		}
		if (strcmp(command, "MODE") == 0) {
			command_identified = TRUE;
			if(sscanf(rest, "%s", &command) <= 0) command[0]='\0';
			upper(command);
			if(strcmp(command, "AUTOMATIC") == 0) {
				control->mode=INDEX_AUTOMATIC;
			}
			else control->mode=INDEX_AUTOMATIC;
		}
		if (strcmp(command, "PEAKS") == 0) {
			command_identified = TRUE;
			files[peakfile_counter] = dps_malloc("files[x]", "read_control", sizeof(peakfile));
			switch(sscanf(rest, "%s %f %f %s", &name, &start, &end)) {
			case EOF:
			case 0: 
				printf("Error in line %d in the control input file: No filename after PEAKS keyword\n", line_count);
				file_ok=FALSE;
				break;
			case 1:
			case 2:
				printf("Error in line %d in the control input file: Not enough or wrong parameters after PEAKS keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			case 3: 
				strcpy(files[peakfile_counter]->filename, name);
				if(start < end) {
					files[peakfile_counter]->phi_start=start/180.0*pi;
					files[peakfile_counter]->phi_end=end/180.0*pi;
					;
				}
				else {
					files[peakfile_counter]->phi_start=end/180.0*pi;
					files[peakfile_counter]->phi_end=start/180.0*pi;
				}
				break;
			default:
				printf("Error in line %d in the control input file: To many parameters after PEAKS keyword\n",
				    line_count);
				file_ok=FALSE;
				break;
			}
			peakfile_counter++;
			if (peakfile_counter >= MAX_PEAKFILES) {
				printf("Number of Peakfiles exceeded (> MAX_PEAKFILES)\n");
				exit(1);
			}
		}
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
						beam_stop=TRUE;
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
						beam_stop=TRUE;
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

				det_set->film_rot_mat[0][0]=cosf(start);
				det_set->film_rot_mat[0][1]=-sinf(start);
				det_set->film_rot_mat[1][0]=-(det_set->film_rot_mat[0][1]);
				det_set->film_rot_mat[1][1]=det_set->film_rot_mat[0][0];

				start=-start;

				det_set->film_rot_mat_inv[0][0]=cosf(start);
				det_set->film_rot_mat_inv[0][1]=-sinf(start);
				det_set->film_rot_mat_inv[1][0]=-(det_set->film_rot_mat_inv[0][1]);
				det_set->film_rot_mat_inv[1][1]=det_set->film_rot_mat_inv[0][0];
				film_rot = TRUE;
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
	/* Check if all necessary parameters are given */
	if (control->d_max == -1.0) {
		printf("Error in the control input file: Keyword CELL is missing\n");
		file_ok=FALSE;
	};
	if (control->om_filename == NULL) {
		printf("Error in the control input file: Keyword OMFILE is missing\n");
		file_ok=FALSE;
	};
	if ((control->resmax == -1.0) || (control->resmin == -1.0)) {
		printf("Error in the control input file: Keyword RESOLUTION is missing\n");
		file_ok=FALSE;
	};
	if (control->mode == -1.0) {
		printf("Error in the control input file: Keyword MODE is missing\n");
		file_ok=FALSE;
	};
	if (det_set->lambda == -1.0) {
		printf("Error in the control input file: Keyword WAVELENGTH is missing\n");
		file_ok=FALSE;
	};
	if (det_set->mm_per_raster == -1.0) {
		printf("Error in the control input file: Keyword RASTERSIZE is missing\n");
		file_ok=FALSE;
	};
	if (det_set->dist_mm == -1.0) {
		printf("Error in the control input file: Keyword DISTANCE is missing\n");
		file_ok=FALSE;
	};
	if (beam_stop == FALSE) {
		printf("Error in the control input file: Keyword BEAMSTOP is missing\n");
		file_ok=FALSE;
	};
	if (film_rot == FALSE) {
		printf("Error in the control input file: Keyword FILM_ROTATION is missing\n");
		file_ok=FALSE;
	};
	if (det_set->y_scale == 0.0) {
		printf("Error in the control input file: Keyword SSCALE is missing\n");
		file_ok=FALSE;
	};
	if (peakfile_counter == 0) {
		printf("Error in the control input file: Keyword PEAKS is missing\n");
		file_ok=FALSE;
	};

	control->peakfiles = calloc(peakfile_counter,  sizeof(peakfile *));
	for(i=0;i<peakfile_counter;i++) {
		control->peakfiles[i]=files[i];
	}
	control->num_peakfiles=peakfile_counter;
	free(files);
#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL DO_READ left...\n");
#endif
	return(file_ok);
}


Boolean read_control(char *filename, index_control *control, detector_setting *det_set)

{

	FILE *in_file;
	Boolean status;

#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL entered...\n");
#endif
	in_file = fopen(filename, "r");
	if (in_file == NULL) {
		printf("Can't open %s!! Exit/n", filename);
		return(FALSE);
	};
	status = read_control_do_read(in_file, control, det_set);
	fclose(in_file);
#ifdef READ_CONTROL_DEBUG
	printf("READ_CONTROL left...\n");
#endif
	return(status);

}

/* read_multiple_rs(): Reads a file with the filenames, start and end phi-values,
			and the type of the file. It then calls the appropiate 
			peak file read routine and pastes the resulting peak
			lists together. */

int read_multiple_rs(index_control *control, struct rs_coord **rs)

{
	double pi = M_PI;
	int i, j;
	int tmp_num;
	struct rs_coord *rs_temp;
	int num=0;

#ifdef READ_MULTIPLE_RS_DEBUG
	printf("READ_MULTIPLE_RS entered...\n");
#endif
	rs_temp=NULL;
	for (i=0; i < control->num_peakfiles; i++) {
		tmp_num = read_rs(control->peakfiles[i]->filename, &rs_temp);
		*rs = realloc(*rs, (num+tmp_num)*sizeof(rs_coord));
		for(j=num;j<num+tmp_num;j++) {
			(*rs)[j]=rs_temp[j-num];
			(*rs)[j].phi_start=control->peakfiles[i]->phi_start;
			(*rs)[j].phi_end=control->peakfiles[i]->phi_end;
		}
		printf("Read %d reflections from file %s\n", tmp_num,control->peakfiles[i]->filename);
		num=num+tmp_num;
	}
	free(rs_temp);
#ifdef READ_MULTIPLE_RS_DEBUG
	printf("READ_MULTIPLE_RS left with %d reflections...\n",num);
#endif
	return(num);

}
