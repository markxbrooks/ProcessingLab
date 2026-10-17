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
 *                         File: ind_transform.c                       *
 *                                                                     *
 *=====================================================================*/

/*  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =*
 *                                                                     *
 * Modifications:                                                      *
 *                                                                     *
 * 31-Jul-1998 Harry Powell, MRC-LMB to recognize if it's being called *
 *             from MOSFLM; also put a #define MOSFLM_INDEX in the     *
 *             "util.h" file in this directory.                        *
 *                                                                     *
 *  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

/* rs_2_xyz() transforms 'number' rs-film coordinates to xyz rec.lattice coordinates 	*
 *			  by using the triangular mechanism. The beamstop corrections 	*
 *			  (beam_stop_r and beam_stop_s) will be applyed.Input arrays r  *
 * 			  and s, output array x, y, z. 			*/

void rs_2_xyz(int number, detector_setting det_set, struct rs_coord rs[],
struct vector xyz[],int called_from_mosflm)

{
	int i;				/* counter */

	double dist_sqr;		/* Crystal Detector Distance in mm^2 */
	double lambda_inv;		/* Inverse of wavelength */
	float r_tmp, s_tmp;		/* temporaery variables for r,s after
							   beamstop correction 	*/
	float X, Y;			/* peak position in camera coordinates in mm */
	double X_sqr, Y_sqr;		/* X,Y squared		*/
	double double_rec_length;	/* 1/(det_set.lambda*(sqrt(X^2+Y^2+dist^2))) */
	float phi;			/* Rotation angle */
	struct matrix_3x3 rotmat;	/* Rotation matrix */
	struct vector xyzo;		/* Reciprocal lattice coordinates before rotation
					   to phi == 0 degrees. */
#ifdef RS_2_XYZ_DEBUG
	FILE *help;
	help=fopen("xyz.read","w");
	printf("RS_2_XYZ entered...\n");
#endif

	dist_sqr  = (double)det_set.dist_mm*(double)det_set.dist_mm;
	lambda_inv     = 1.0/(double)det_set.lambda;

#if MOSFLM_INDEX
	printf("MOSFLM_INDEX is TRUE \n");
#else
	printf("MOSFLM_INDEX is FALSE \n");
#endif
	for(i=0;i<number;i++) {

		if (called_from_mosflm == 0) {
			/* Convert to coords relative to direct beam position */
			r_tmp = rs[i].r - det_set.beam_stop_r;
			s_tmp = rs[i].s - det_set.beam_stop_s;

			/* Convert to mm, and invert S if y_scale is negative */
			r_tmp=r_tmp*det_set.mm_per_raster;
			s_tmp=s_tmp*det_set.mm_per_raster*det_set.y_scale;
		}
		else {
			/* If called from MOSFLM, the rs[i] values are really XY[i] values!*/
			r_tmp = rs[i].r;
			s_tmp = rs[i].s;
		}

		X = det_set.film_rot_mat[0][0] * r_tmp + det_set.film_rot_mat[0][1] * s_tmp;
		Y = det_set.film_rot_mat[1][0] * r_tmp + det_set.film_rot_mat[1][1] * s_tmp;
		X_sqr = (double)X*(double)X;
		Y_sqr = (double)Y*(double)Y;
		double_rec_length = sqrt(X_sqr + Y_sqr + dist_sqr);
		double_rec_length = lambda_inv/(double_rec_length);

		/* Calculate rec. lattice coordinates after:
	   x = r/(det_set.lambda*distance_crystal-peak_on_film)
	   y = s/(det_set.lambda*distance_crystal-peak_on_film)
	   z = distance_crystal-peak_on_film/(det_set.lambda*distance_crystal-peak_on_film)-1/det_set.lambda */

		xyzo.x = (float)(double_rec_length * X);
		xyzo.y = (float)(double_rec_length * Y);
		xyzo.z = (float)(double_rec_length * (double)det_set.dist_mm - lambda_inv);

		/* Now rotate back to phi = 0 degress */
		phi=-(rs[i].phi_end+rs[i].phi_start)/2.0;
		rotmatY_from_angle(phi, &rotmat);
		mat_vec_mul_3x3(rotmat, xyzo, &(xyz[i]));
#ifdef RS_2_XYZ_DEBUG
		printf("%d r s: %f %f x y z: %f %f %f \n",i,rs[i].r,rs[i].s,xyz[i].x,xyz[i].y,xyz[i].z);
		fprintf(help,"%f %f %f \n",xyz[i].x,xyz[i].y,xyz[i].z);
		printf("phistart %f phiend %f phi %f\n",rs[i].phi_start, rs[i].phi_end,   phi);
		print_3x3(rotmat);
#endif
	};

#ifdef RS_2_XYZ_DEBUG
	fclose(help);
	printf("RS_2_XYZ finished...\n");
#endif
}

/* cut_xyz() on a given list of xyz generate a list of xyz in a given	*
 *		resolution range.		 			*/

int  cut_xyz(int number, struct vector xyz[], float resmin, float resmax,
struct vector *xyz_new)

{
	int i;		/* counter */
	float len;				/* Length of a vector */
	float max, min;

	float resmax_inv_sqr, resmin_inv_sqr;	/* (1/resmin)^2....	*/
	int number_xyz = 0;	/* Number of xyz in resolution range, index for the array */

#ifdef CUT_XYZ_DEBUG
	printf("CUT_XYZ entered with %d vectors...\n",number);
#endif
	resmax_inv_sqr = (1.0/resmax)*(1.0/resmax);
	resmin_inv_sqr = (1.0/resmin)*(1.0/resmin);

	/* Initialize variable for maximum and minimum resolution */
	min = 0.0;
	max = 100000.0;

	for(i=0;i<number;i++) {
		len = xyz[i].x * xyz[i].x + xyz[i].y * xyz[i].y + xyz[i].z * xyz[i].z;
		/* Look for reflections in selected resolution range */
		if ((len < resmax_inv_sqr) && (len > resmin_inv_sqr)) {
			cp_vectorf(xyz[i],&xyz_new[number_xyz]);
			number_xyz++;
			if(len > min) min = len;
			if(len < max) max = len;
#ifdef CUT_XYZ_DEBUG
			printf(" %f %f %f with len %f ok\n", xyz[i].x, xyz[i].y, xyz[i].z, len);
#endif
		}
		else {
#ifdef CUT_XYZ_DEBUG
			printf("rejected\n");
#endif
		}
	}
	printf("%d reciprocal lattice points between %5.1f and %5.1f Angstroem resolution\n",number_xyz,1.0/sqrt(max),1.0/
	    sqrt(min));

#ifdef CUT_XYZ_DEBUG
	printf("CUT_XYZ finished with %d vectors...\n",number_xyz);
#endif
	return(number_xyz);
}

/* xyz_2_rs_long() does the same as xyz_2_rs(), but returns also lists
                   of temporary parameters during the calculation. These 
                   might be needed in a refinement program for 
                   determinating derivatives. Unused as of version 
                   DPS 0.91						*/

void xyz_2_rs_long(struct vector *xyz, int number, detector_setting det_set,
float lambda, float osc_start,  float osc_end, struct rs_coord *rs, float *phi_r, int *sign_int,
struct vector *R)

{
	int i;			/* counter */
	float dummy = 0;

	for(i = 0; i <= number; i++) {
		xyz_2_rs_single(xyz[i], det_set, lambda, &dummy, &osc_start,  &osc_end, &rs[i], &phi_r[i], &sign_int[i],
		    &R[i]);
	}
}

/* xyz_2_rs() transforms given xyz coordinates to film coordinates 
	      using xyz_2_rs_single().					*/

void xyz_2_rs(struct vector *xyz, int number, detector_setting det_set,
float lambda, float osc_start,  float osc_end, struct rs_coord *rs)

{
	int i;			/* counter */
	float phi_r;
	int sign_int;
	struct vector R;
	float mos;

	for(i = 0; i <= number; i++) {
		xyz_2_rs_single(xyz[i], det_set, lambda, &mos, &osc_start,  &osc_end, &rs[i], &phi_r, &sign_int, &R);
	}
}

/* calc_phis(): calculates phi of intersection and  starting phi for calculation *
 *              of film coordinates from reciprocal space coordinates. Used in 
 *  		xyz_2_rs_single.					    */

int calc_phis(struct vector *xyz, float *lambda, float *rho,  float *len_sqr, 
float *phi_r1, float *phi_r2)

{
	double pi = M_PI;
	float arg;		/* Argument in acos for phi_int calculation */
	float phi_int;		/* angle origin/intersection ewald */
	float phi_0;		/* angle origion/xyz		   */

	/* First calculate the two phi's for the necessary 
		   rotation to move xyz to the two intersetion points 
		   with the ewald sphere				*/

	/* Calculate the square of the length of xyz		*/

	*rho=xyz->z*xyz->z+xyz->x*xyz->x;
	*len_sqr=xyz->y * xyz->y + *rho;
	*rho = sqrt(*rho);


	/* Calculate the angle from the origin direction to the 
		   intersection with ewald sphere. First check the
		   argument of arccosine if there is an intersection.
	  	   If not, set r and s to 0 and go to next reflection	*/
	arg = -(*lambda*(*len_sqr))/(2*(*rho));
	if ( fabs(arg) > 1 ) {
		return(2);
	}
	phi_int=acos(arg)/pi*180.;

	/* Calculate the angle from the origin direction to xyz */
	phi_0=atan2(xyz->x,xyz->z)/pi*180.;
	if (phi_0 < 0 ) phi_0 = 360 + phi_0;

	/* phi_r1 and phi_r2 are the necessary rotation to one 
		   of the intersections.				*/
	*phi_r1=phi_int - phi_0;
	*phi_r2=360.0 - phi_int - phi_0;
	return(0);

}

/* rs_from_phis(): calculate film R and S coordinates from given phi_r's
                   Used in xyz_2_rs-single.  				*/

void rs_from_phis(struct vector *xyz, float *phi_r, float *lambda, 
detector_setting *det_set, struct vector *R, struct rs_coord *rs)

{
	float X, Y;		/*  Camera coordinates of peak in mm */
	float r_tmp, s_tmp;	/* r and s before application of scalefactors and beam_stop */

	/* Now calculate the new xyz after rotation by phi_r	*/
	R->x=cos(*phi_r)*xyz->x*(*lambda)+sin(*phi_r)*xyz->z*(*lambda);
	R->y=xyz->y*(*lambda);
	R->z=-sin(*phi_r)*xyz->x*(*lambda)+cos(*phi_r)*xyz->z*(*lambda);

	/* The last step is the calculation of the film coordinates
		   using the triangular					*/
	X=R->x/(R->z+1)*det_set->dist_mm;
	Y=R->y/(R->z+1)*det_set->dist_mm;
	r_tmp = det_set->film_rot_mat_inv[0][0] * X + det_set->film_rot_mat_inv[0][1] * Y;
	s_tmp = det_set->film_rot_mat_inv[1][0] * X + det_set->film_rot_mat_inv[1][1] * Y;

	rs->r=r_tmp/(det_set->mm_per_raster)+det_set->beam_stop_r;
	rs->s=s_tmp/(det_set->mm_per_raster*det_set->y_scale)+det_set->beam_stop_s;


}

/* xyz_2_rs_single_use_mos() transforms given xyz coordinate to film coordinate by	*
 *	the algorithms described in the Ph.D. thesis of Robert 		*
 *	Bolotovski. Input is the one xyz coordinate (r/s = 0/0 for 	*
 *	not observable reflections) the detector setting and the wavelength.	*
 *	On output also the angle phi_r as well as the sign of phi_int	*
 *	is given. The later is needed for the derivatives in LS		*/

int xyz_2_rs_single_use_mos(struct vector xyz, detector_setting det_set,
float lambda, float *mos, float *osc_start,  float *osc_end,  struct rs_coord *rs, 
float *phi_r, int *sign_int,
struct vector *R)

{
	double pi = M_PI;
	float len_sqr;		/* length of a xyz vector */
	float phi_r1, phi_r2;	/* Rotation xyz/intersection ewald 1 or 2 */
	float rho;		/* factored needed in calculation, look at Thesis of Robert Bolotovsky */
	float eff_mos;
	float osc1, osc2;
	int back_value;


	if(calc_phis(&xyz, &lambda,  &rho, &len_sqr, &phi_r1, &phi_r2) == 2) {
		rs->r = 0.0;
		rs->s = 0.0;
		*sign_int=0;
		return(2);
	}
	/* Check if reflection could be observed */
	eff_mos=2*atan(tan(*mos/2)/fabs(rho/len_sqr))/pi*180;
	osc1=*osc_start - eff_mos;
	osc2=*osc_end + eff_mos;
	back_value=1;
	if (phi_r1 > osc1 && phi_r1 < osc2) {
		*phi_r = phi_r1;
		*sign_int=1;
		back_value=0;
	}
	if (phi_r2 > osc1 && phi_r2 < osc2) {
		*phi_r = phi_r2;
		*sign_int=-1;
		back_value=0;
	}
	*phi_r=*phi_r/180.*pi;

	rs_from_phis(&xyz, phi_r, &lambda, &det_set, R, rs);
	return(back_value);
}

/* xyz_2_rs_single() transforms given xyz coordinate to film coordinate by      *
 *      the algorithms described in the Ph.D. thesis of Robert          *
 *      Bolotovski. Input is the one xyz coordinate (r/s = 0/0 for      *
 *      not observable reflections) the detector setting and the wavelength.    *
 *      On output also the angle phi_r as well as the sign of phi_int   *
 *      is given. The later is needed for the derivatives in LS. The 
 *      return value is 0 if the reflection could be predicted on the   *
 *      film and 1 if not.					        */

int xyz_2_rs_single(struct vector xyz, detector_setting det_set,
float lambda, float *mos, float *osc_start,  float *osc_end,  
struct rs_coord *rs, 
float *phi_r, int *sign_int,
struct vector *R)

{
	double pi = M_PI;
	float len_sqr;		/* length of a xyz vector */
	float phi_r1, phi_r2;	/* Rotation xyz/intersection ewald 1 or 2 */
	float rho;		/* factored needed in calculation, look at Thesis of Robert Bolotovsky */
	float eff_mos;
	int back_value;
	float osc_middle, osc_range_halve;
	float rot1, rot2;

	if(calc_phis(&xyz, &lambda,  &rho, &len_sqr, &phi_r1, &phi_r2) == 2) {
		rs->r = 0.0;
		rs->s = 0.0;
		*sign_int=0;
		*mos=0.0;
		return(2);
	}
	back_value=0;
	/* Check if reflection could be observed */
	osc_middle=(*osc_start+(*osc_end))/2.0;
	osc_range_halve=fabs(*osc_end-*osc_start);
	rot1=fabs(osc_middle-phi_r1);
	rot2=fabs(osc_middle-phi_r2);
	if(rot1 < rot2) {
		*phi_r=phi_r1;
		*sign_int=1;
	}
	else {
		*phi_r=phi_r2;
		*sign_int=-1;
	}
	rot1=*phi_r-osc_range_halve;
	if (rot1 > 0.0) {
		*mos=2*tan((rot1/180.*pi)/2*fabs(rho/len_sqr))/pi*180.;
		back_value=1;
	}
	else {
		*mos=0.0;
		back_value=0;
	}
	*phi_r=*phi_r/180.*pi;

	rs_from_phis(&xyz, phi_r, &lambda, &det_set, R, rs);
	return(back_value);
}
