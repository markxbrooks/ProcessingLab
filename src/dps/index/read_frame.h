#ifndef READ_FRAME_H

#include "util_def.h"

#undef READ_FRAME_DEBUG

#define LENGTH_MAGIC 10 /* Length of the magic string */

typedef unsigned short od_type;	/* Type of optical density pixel values */
typedef int FrameType;  /* Type for the type of the fram. See floowing define too */

/* Define the different frame types (For use with FrameType) */

#define FUJI            1
#define FUJI_RAW        2    /* without applying the 10^(i/256) correction */
#define PGM_ASCII       3
#define KODAK           4
#define APS3000         5
#define RAXIS2_100      6
#define RAXIS4_100      7
#define QUANTUM1        8   /* ADSC single module */
#define QUANTUM4        9   /* ADSC 2x2 */
#define	MAR180		10
#define	MAR300		11
#define	MAR345		12


typedef struct {
	int 		lenx, leny;	/* Length of Frame in X (fast) and Y (slow) direction */
	od_type		minod;
	od_type 	maxod;		/* Maximum allowed optical density */
        FrameType 	type;		/* type of frame */	
	int		header;		/* Length of header */
	od_type		overload;
} FrameInfo;

int read_frame(char *, FrameInfo *, od_type **);

#define READ_FRAME_H
#endif
