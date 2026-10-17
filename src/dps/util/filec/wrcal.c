#include <stdio.h>
#include <string.h>
#include "filec.h"


int wrcal ( char* filename, struct readcalp *c )
{
	int istat;
	char info[2048];
	char s[80];


	/*
   * These fields MUST be in the calibration file
   */
	clrhd (info);

	puthd ("TYPE", "calibration_file", info );
	sprintf (s, "%d", c->x_size);
	puthd ("X_SIZE", s, info );
	sprintf (s, "%d", c->y_size);
	puthd ("Y_SIZE", s, info );

	sprintf (s, "%d", c->xint_start);
	puthd ("XINT_START", s, info );
	sprintf (s, "%d", c->yint_start);
	puthd ("YINT_START", s, info );
	sprintf (s, "%d", c->xint_step);
	puthd ("XINT_STEP", s, info );
	sprintf (s, "%d", c->yint_step);
	puthd ("YINT_STEP", s, info );
	sprintf (s, "%d", c->xinv_start);
	puthd ("XINV_START", s, info );
	sprintf (s, "%d", c->yinv_start);
	puthd ("YINV_START", s, info );
	sprintf (s, "%d", c->xinv_step);
	puthd ("XINV_STEP", s, info );
	sprintf (s, "%d", c->yinv_step);
	puthd ("YINV_STEP", s, info );

	sprintf (s, "%.8e", c->pscale);
	puthd ("PSCALE", s, info );

	sprintf (s, "%f", c->pixsiz);
	puthd ("PIXEL_SIZE", s, info );

	sprintf (s, "%d", c->bad_flag);
	puthd ("BAD_FLAG", s, info );

	/*
   * These fields SHOULD be in the calibrate file
   *      (MADNES required them)
   */
	sprintf (s, "%f", c->x_beam);
	puthd ("X_BEAM", s, info );
	sprintf (s, "%f", c->y_beam);
	puthd ("Y_BEAM", s, info );
	sprintf (s, "%d", c->x_size);

	sprintf (s, "%f", c->spacing);
	puthd ("SPACING", s, info );
	sprintf (s, "%f", c->x_scale);
	puthd ("X_SCALE", s, info );
	sprintf (s, "%f", c->y_scale);
	puthd ("Y_SCALE", s, info );

	/*
   * These field are not required by most routines
   */
	if (c->x_center)
	{
		sprintf (s, "%f", c->x_center);
		puthd ("X_CENTER", s, info );
	}
	if (c->y_center)
	{
		sprintf (s, "%f", c->y_center);
		puthd ("Y_CENTER", s, info );
	}
	if (c->x_pt_center)
	{
		sprintf (s, "%f", c->x_pt_center);
		puthd ("X_PT_CENTER", s, info );
	}
	if (c->y_pt_center)
	{
		sprintf (s, "%f", c->y_pt_center);
		puthd ("Y_PT_CENTER", s, info );
	}
	if (c->ratio)
	{
		sprintf (s, "%f", c->ratio);
		puthd ("RATIO", s, info );
	}
	if (c->ver_slope)
	{
		sprintf (s, "%f", c->ver_slope);
		puthd ("VER_SLOPE", s, info );
	}
	if (c->horz_slope)
	{
		sprintf (s, "%f", c->horz_slope);
		puthd ("HORZ_SLOPE", s, info );
	}
	if (c->a)
	{
		sprintf (s, "%f", c->a);
		puthd ("RADIAL_A", s, info );
	}
	if (c->a1)
	{
		sprintf (s, "%f", c->a1);
		puthd ("RADIAL_A1", s, info );
	}
	if (c->b)
	{
		sprintf (s, "%f", c->b);
		puthd ("RADIAL_B", s, info );
	}
	if (c->c)
	{
		sprintf (s, "%f", c->c);
		puthd ("RADIAL_C", s, info );
	}

	istat = wrhead ( filename, info );

	return istat;
}
