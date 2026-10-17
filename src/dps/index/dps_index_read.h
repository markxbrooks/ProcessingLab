/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* dps_index_read.h:  Header File to dps_index_read.c  */

#ifndef DPS_INDEX_READ

#include "util.h"
#include "index_def.h"

#undef READ_CONTROL_DEBUG
#undef READ_MULTIPLE_RS_DEBUG

typedef short peaksfiletype;
#define PEAKSFILETYPE_DPS 0
typedef struct {                /* Structure for the peak file information */
                char filename[100];     /* filename */
                float phi_start;        /* Starting phi */
                float phi_end;          /* Ending phi */
                peaksfiletype type;     /* Type fo the peaks file */
        } peakfile;

typedef struct {
        Index_Mode      mode;   /* Mode of Indexing */
	char		*data_filename;  /* String for the detector data filename */
        char            *om_filename;   /* String for the OM file filename */
        float           d_max;  /* maximum expected cell edge length */
        float           resmin; /* Minimum resolution to use */
        float           resmax; /* Maximum resolution to use */
        int             num_peakfiles;  /* Number of peakfiles to be used */
        peakfile	**peakfiles;	/* The peakfiles */
} index_control;


Boolean read_control_do_read(FILE *,  index_control *, detector_setting *);
Boolean read_control(char *, index_control *, detector_setting *);
int read_multiple_rs(index_control *, struct rs_coord **rs);

#define DPS_INDEX_READ

#endif
