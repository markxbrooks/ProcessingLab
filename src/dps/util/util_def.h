#ifndef UTIL_DEF_H

#ifndef _XtIntrinsic_h
typedef char Boolean;
#define TRUE 1
#define FALSE 0
#endif

#define N_XYZ 5000		/* maximum number of rec. lattic coordinates that
						   can be handled. */
#define N_RS  5000		/* maximum number of film peak coordinates that 
						   can be handled. Must be larger or equal to N_XYZ */

typedef struct {
    float lambda;	    /* Wavelength */
    float dist_mm;	    /* distance crystal detector in mm */
    float beam_stop_r;	    /* beamstop position in r direction in rasterpoints*/
    float beam_stop_s;	    /* beamstop position in s driection in rasterpoints*/
    float mm_per_raster;    /* mm per rasterstep */
    float y_scale;	    /* if pixel is larger in y than in x than y_scale > 1 */
    float film_rot;	    /* film rotation in degree, 0 means film R is parallel to
			       laboratory X */
    float film_rot_mat[2][2];	    /* Rotation matrix from film_rot */
    float film_rot_mat_inv[2][2];   /* Invers of the rotation matrix */
} detector_setting;	

typedef struct rs_coord{
	float	r;         /* R coordinate in pixels */
	float	s;         /* S coordinate in pixels */
	float 	phi_start; /* Starting phi angle of oscillation of the frame the peak was found on */
	float   phi_end;   /* Ending phi angle of oscillation of the frame the peak was found on */
} rs_coord;


#define UTIL_DEF_H

#endif
