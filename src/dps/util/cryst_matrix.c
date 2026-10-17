/* $Id: cryst_matrix.c,v 1.5 2012/09/14 10:52:02 andrew Exp $ */
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
 *                         File: cryst_matrix.c                        *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <math.h>
#include "cryst_matrix.h"

/* om_from_vectors(): Generate orientation matrice from 3 direct cell	*
 *		vectors.						*/

void om_from_vectors(struct vector a, struct vector b, struct vector c,
struct matrix_3x3 *om)
{
	float volume;		/* Volume of the cell */
	struct vector tmp_vec;	/* temporary vector */

	/* Calculate from a,b and c, a*, b* and c* with e.g
	   a* = (b x c)/volume					*/

	volume=vec_volumef(a,b,c);
	cross_prodf(b,c,&tmp_vec);
	om->r[0].x=tmp_vec.x/volume;
	om->r[0].y=tmp_vec.y/volume;
	om->r[0].z=tmp_vec.z/volume;
	cross_prodf(c,a,&tmp_vec);
	om->r[1].x=tmp_vec.x/volume;
	om->r[1].y=tmp_vec.y/volume;
	om->r[1].z=tmp_vec.z/volume;
	cross_prodf(a,b,&tmp_vec);
	om->r[2].x=tmp_vec.x/volume;
	om->r[2].y=tmp_vec.y/volume;
	om->r[2].z=tmp_vec.z/volume;
}

/* xyz_2_hkl(): From a give list of xyz and an orientation matrix, 	*
 *		generate two lists of hkl values, on with fractional hkl*
 *		and one rounded to integer values.			*/

void xyz_2_hkl(struct matrix_3x3 om, struct vector xyz[], int number, 
struct hkl_int *hkl, struct hkl_fract *hkl_calc)

{
	int i;				/* counter */
	struct matrix_3x3 om_inv;	/* inverse of orientation matrix */

	(void)invert_3x3(om,&om_inv);
	for(i = 0; i < number; i++) {
		hkl_calc[i].h=xyz[i].x*om_inv.r[0].x+xyz[i].y*om_inv.r[1].x+xyz[i].z*om_inv.r[2].x;
		hkl_calc[i].k=xyz[i].x*om_inv.r[0].y+xyz[i].y*om_inv.r[1].y+xyz[i].z*om_inv.r[2].y;
		hkl_calc[i].l=xyz[i].x*om_inv.r[0].z+xyz[i].y*om_inv.r[1].z+xyz[i].z*om_inv.r[2].z;
		hkl[i].h=(int) rint((double)hkl_calc[i].h);
		hkl[i].k=(int) rint((double)hkl_calc[i].k);
		hkl[i].l=(int) rint((double)hkl_calc[i].l);
/*                printf("h,k,l  %d %d %d h' k' l' %f %f %f\n",hkl[i].h,hkl[i].k,hkl[i].l, *
		  hkl_calc[i].h,hkl_calc[i].k,hkl_calc[i].l);                              */
	}
}

/* hkl_2_xyz(): From a give list of hkl and an orientation matrix, 	*
 *		generate a list of xyz coordinates.			*/

void hkl_2_xyz(struct matrix_3x3 om, struct hkl_int hkl[], int number, 
struct vector *xyz)

{
	int i;		/* counter */

	for(i = 0; i <= number; i++) {
		xyz[i].x=hkl[i].h*om.r[0].x+hkl[i].k*om.r[1].x+hkl[i].l*om.r[2].x;
		xyz[i].y=hkl[i].h*om.r[0].y+hkl[i].k*om.r[1].y+hkl[i].l*om.r[2].y;
		xyz[i].z=hkl[i].h*om.r[0].z+hkl[i].k*om.r[1].z+hkl[i].l*om.r[2].z;
	}
}

/* hkl_2_xyz_single(): For a given hkl calculate the xyz using om	*/

void hkl_2_xyz_single(struct matrix_3x3 om, struct hkl_int hkl, struct vector *xyz)

{
	xyz->x=hkl.h*om.r[0].x+hkl.k*om.r[1].x+hkl.l*om.r[2].x;
	xyz->y=hkl.h*om.r[0].y+hkl.k*om.r[1].y+hkl.l*om.r[2].y;
	xyz->z=hkl.h*om.r[0].z+hkl.k*om.r[1].z+hkl.l*om.r[2].z;
}

/* dircell_from_om(): Generate direct cell constants from orientation matrix */

void dircell_from_om(struct matrix_3x3 om, struct cell_struct *cell)
{
	struct cell_struct rec_cell; /* temporary reciprocal cell */

	reccell_from_om(om, &rec_cell);
	dircell_from_reccell(rec_cell, cell);

}

/* reccell_from_om(): Generate direct cell constants from orientation matrix */

void reccell_from_om(struct matrix_3x3 om, struct cell_struct *cell)
{
	double help;

	/* First, calculate length of a*, b*, c* */
	cell->a=vec_lenf(om.r[0]);
	cell->b=vec_lenf(om.r[1]);
	cell->c=vec_lenf(om.r[2]);

	/* calculate angle by: cos gamma = (a* dot b*)/(len_a * len_b) */
	help=scalar_prodf(om.r[1],om.r[2])/(cell->b * cell->c);
	if ( help > 1 ) help = 1;
	else if ( help < -1 ) help = -1 ;
	cell->alpha = acos(help);
	help=scalar_prodf(om.r[0],om.r[2])/(cell->a * cell->c);
	if ( help > 1 ) help = 1;
	else if ( help < -1 ) help = -1 ;
	cell->beta  = acos(help);
	help=scalar_prodf(om.r[0],om.r[1])/(cell->a * cell->b);
	if ( help > 1 ) help = 1;
	else if ( help < -1 ) help = -1 ;
	cell->gamma = acos(help);


}

/* dircell_from_reccell(): Converts direct cell into reciprocal cell and */
/* vice versa. */

void dircell_from_reccell(struct cell_struct rec_cell, struct cell_struct *cell)
{
	double help;

	/* Formulas from Giacovazzo (p. 64); help = V/(abc) */

	help = sqrt(1 - pow(cos(rec_cell.alpha),2) - pow(cos(rec_cell.beta),2)
	    - pow(cos(rec_cell.gamma),2) + 2 * cos(rec_cell.alpha)
	    * cos(rec_cell.beta) * cos(rec_cell.gamma));
	cell->a = sin(rec_cell.alpha)/(rec_cell.a * help);
	cell->b = sin(rec_cell.beta)/(rec_cell.b * help);
	cell->c = sin(rec_cell.gamma)/(rec_cell.c * help);

	cell->alpha = acos((cos(rec_cell.beta) * cos(rec_cell.gamma) -
	    cos(rec_cell.alpha))/(sin(rec_cell.beta)  * sin(rec_cell.gamma)));
	cell->beta  = acos((cos(rec_cell.alpha) * cos(rec_cell.gamma) - 
	    cos(rec_cell.beta))/(sin(rec_cell.alpha) * sin(rec_cell.gamma)));
	cell->gamma = acos((cos(rec_cell.alpha) * cos(rec_cell.beta) - 
	    cos(rec_cell.gamma))/(sin(rec_cell.alpha) * sin(rec_cell.beta)));

}


/* setangles_from_reccell(): Calculate the setting angles from the given unit_cell	*
 *			and the orientation matrix. */

void setangles_from_reccell(struct cell_struct cell, struct matrix_3x3 om, 
struct setangles_struct *setangels)
{
	struct matrix_3x3 b, b_inv; /* Orthogonalization matrix and it's inverse */
	struct matrix_3x3 umatrix;  /* Combined Rotation matrices */

	/* First generate the orthogonalisation matrix from the reciprocal
	   cell parameters */

	b.r[0].x=cell.a;
	b.r[0].y=0.0;
	b.r[0].z=0.0;
	b.r[1].x=((cell.b)*cos(cell.gamma));
	b.r[1].y=((cell.b)*sin(cell.gamma));
	b.r[1].z=0.0;
	b.r[2].x=((cell.c)*cos(cell.beta));
	b.r[2].y=((cell.c)*(cos(cell.alpha)-cos(cell.beta)*cos(cell.gamma))/sin(cell.gamma));
	b.r[2].z=cell.c*pow((1.0-pow(cos(cell.alpha),2)-pow(cos(cell.beta),2)-pow(cos(cell.gamma),2) +
	    2*cos(cell.alpha)*cos(cell.beta)*cos(cell.gamma)),0.5)/sin(cell.gamma);

	/* multiplying the orientation matrix with the inverse of the 
	   orthogonalisation matrix results in the combined rotation
	   matrix */

	invert_3x3(b, &b_inv);
	mat_mul_3x3(om, b_inv, &umatrix);

	/* By analysing the umatrix we get the three rotation angles */

	setangels->y=-asin(umatrix.r[0].z);

	/* Don't know if this decission is nessecary!!! */
	if (umatrix.r[0].z > -1 ) {
		setangels->x=-atan2(umatrix.r[1].z,umatrix.r[2].z);
		setangels->z=-atan2(umatrix.r[0].y,umatrix.r[0].x);
	}
	else {
		setangels->x=-(acos(umatrix.r[1].y)+asin(umatrix.r[1].x))/2.0;
		setangels->z=-(acos(umatrix.r[1].y)-asin(umatrix.r[1].x))/2.0;
	}
}

/* setangles_from_reccell_DENZO(): THIS ROUTINE IS UNCOMPLETED AND DOES NOT WORK!!!! */

void setangles_from_reccell_DENZO(struct cell_struct cell, struct matrix_3x3 om, 
struct setangles_struct *setangels)
{
	struct matrix_3x3 b, b_inv; /* Orthogonalization matrix and it's inverse */
	struct matrix_3x3 umatrix, umatrix_orig, umatrix_den;  /* Combined Rotation matrices */
	struct matrix_3x3 om_den, denzo, help;

	float pi=acos(-1);

	printf("original om\n");
	print_3x3(om);

	/*	clear_3x3(&denzo);
	denzo.r[0].y=-1.0;
	denzo.r[1].x=1.0;
	denzo.r[2].z=1.0;

	mat_mul_3x3(denzo,om,&help);

	clear_3x3(&denzo);

	denzo.r[0].y=1.0;
	denzo.r[1].x=1.0;
	denzo.r[2].z=-1.0;



	mat_mul_3x3(denzo,help,&om_den);


	om_den.r[0].x=om.r[0].x;
	om_den.r[0].y=om.r[0].y;
	om_den.r[0].z=om.r[0].z;
	om_den.r[1].x=om.r[1].x;
	om_den.r[1].y=om.r[1].y;
	om_den.r[1].z=om.r[1].z;
	om_den.r[2].x=om.r[2].x;
	om_den.r[2].y=om.r[2].y;
	om_den.r[2].z=om.r[2].z;

	printf("denzo om\n");
	print_3x3(om_den);
*/
	reccell_from_om(om, &cell);
	B_from_reccell(cell, &b);

	invert_3x3(b, &b_inv);
	printf("b\n");
	print_3x3(b);
	printf("b_inv\n");
	print_3x3(b_inv);
	mat_mul_3x3(om, b_inv, &umatrix_orig);

	printf("umatrix orig\n");
	print_3x3(umatrix_orig);

	clear_3x3(&denzo);
	denzo.r[0].y=1.0;
	denzo.r[1].z=1.0;
	denzo.r[2].x=1.0;
	mat_mul_3x3(denzo,umatrix_orig,&umatrix_den);
	printf("denzo umatrix\n");
	print_3x3(umatrix_den);

	/*	clear_3x3(&denzo);
	denzo.r[0].z=-1.0;
	denzo.r[1].y=1.0;
	denzo.r[2].x=1.0;
	mat_mul_3x3(denzo, umatrix_den, &umatrix);
	printf("analysable umatrix\n");
*/
	umatrix = umatrix_den;
	print_3x3(umatrix);


	/* By analysing the umatrix we get the three rotation angles */

	/*	setangels->y=-asin(-umatrix.r[2].x);

	setangels->x=-atan2(-umatrix.r[2].y,umatrix.r[2].z);
	setangels->z=-atan2(-umatrix.r[1].x,umatrix.r[0].x);
*/
	setangels->y=asin(-umatrix.r[2].z);

	setangels->x=-atan2(-umatrix.r[2].y,-umatrix.r[2].x);
	setangels->z=-atan2(-umatrix.r[1].z,umatrix.r[0].z);

}

/* om_from_reccell_angles(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void om_from_reccell_angles(struct cell_struct cell, 
struct setangles_struct setangels, 
struct matrix_3x3 *om)
{
	struct matrix_3x3 b;	    /* Orthogonalization matrix and it's inverse */
	struct matrix_3x3 umatrix;  /* Combined Rotation matrices */

	/* First generate the orthogonalisation matrix from the reciprocal
	   cell parameters */

	b.r[0].x=cell.a;
	b.r[0].y=0.0;
	b.r[0].z=0.0;
	b.r[1].x=((cell.b)*cos(cell.gamma));
	b.r[1].y=((cell.b)*sin(cell.gamma));
	b.r[1].z=0.0;
	b.r[2].x=((cell.c)*cos(cell.beta));
	b.r[2].y=((cell.c)*(cos(cell.alpha)-cos(cell.beta)*cos(cell.gamma))/sin(cell.gamma));
	b.r[2].z=cell.c*pow((1.0-pow(cos(cell.alpha),2)-pow(cos(cell.beta),2)-pow(cos(cell.gamma),2) +
	    2*cos(cell.alpha)*cos(cell.beta)*cos(cell.gamma)),0.5)/sin(cell.gamma);


	/* Generate the comjbined rotation matrix from the setting angles */

	umatrix.r[0].x=cos(setangels.y)*cos(setangels.z);
	umatrix.r[0].y=cos(setangels.y)*sin(setangels.z);
	umatrix.r[0].z=-sin(setangels.y);
	umatrix.r[1].x=-cos(setangels.x)*sin(setangels.z)+sin(setangels.x)*sin(setangels.y)*cos(setangels.z);
	umatrix.r[1].y=cos(setangels.x)*cos(setangels.z)+sin(setangels.x)*sin(setangels.y)*sin(setangels.z);
	umatrix.r[1].z=sin(setangels.x)*cos(setangels.y);
	umatrix.r[2].x=sin(setangels.x)*sin(setangels.z)+cos(setangels.x)*sin(setangels.y)*cos(setangels.z);
	umatrix.r[2].y=-sin(setangels.x)*cos(setangels.z)+cos(setangels.x)*sin(setangels.y)*sin(setangels.z);
	umatrix.r[2].z=cos(setangels.x)*cos(setangels.y);

	/* Combine the matrices to the orientation matrix */
	mat_mul_3x3(umatrix, b, om);

}

/* rotmatX_from_angle(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void rotmatX_from_angle(float angle, struct matrix_3x3 *rot)

{
	clear_3x3(rot);
	rot->r[0].x=1;
	rot->r[1].y=cos(angle);
	rot->r[1].z=sin(angle);
	rot->r[2].y=-sin(angle);
	rot->r[2].z=cos(angle);
}

/* rotmatY_from_angle(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void rotmatY_from_angle(float angle, struct matrix_3x3 *rot)

{
	clear_3x3(rot);
	rot->r[1].y=1;
	rot->r[0].x=cos(angle);
	rot->r[2].x=sin(angle);
	rot->r[0].z=-sin(angle);
	rot->r[2].z=cos(angle);
}

/* rotmatZ_from_angle(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void rotmatZ_from_angle(float angle, struct matrix_3x3 *rot)

{
	clear_3x3(rot);
	rot->r[2].z=1;
	rot->r[0].x=cos(angle);
	rot->r[0].y=sin(angle);
	rot->r[1].x=-sin(angle);
	rot->r[1].y=cos(angle);
}

/* B_from_reccell(): Calculate the orthogonalization matrix from the	*
 *                   reciprocal cell and the setting angles. 		*/

void B_from_reccell(struct cell_struct cell, struct matrix_3x3 *b)

{

	b->r[0].x=cell.a;
	b->r[0].y=0.0;
	b->r[0].z=0.0;
	b->r[1].x=((cell.b)*cos(cell.gamma));
	b->r[1].y=((cell.b)*sin(cell.gamma));
	b->r[1].z=0.0;
	b->r[2].x=((cell.c)*cos(cell.beta));
	b->r[2].y=((cell.c)*(cos(cell.alpha)-cos(cell.beta)*cos(cell.gamma))/sin(cell.gamma));
	b->r[2].z=cell.c*pow((1.0-pow(cos(cell.alpha),2)-pow(cos(cell.beta),2)-pow(cos(cell.gamma),2) +
	    2*cos(cell.alpha)*cos(cell.beta)*cos(cell.gamma)),0.5)/sin(cell.gamma);

}

/* B_and_U_from_reccell_and_om(): Calculate orthogonalization matrix B 
 * 	and rotation matrix U from reccell and om			*/

void B_and_U_from_om(struct matrix_3x3 om, struct matrix_3x3 *b, struct matrix_3x3 *u)

{
	struct cell_struct reccell;
	struct matrix_3x3 b_inv;

	/* First get the reciprocal cell from the OM */
	reccell_from_om(om, &reccell);

	/* extract B from the reciprocal cell */
	B_from_reccell(reccell, b);


	/* multiplying the orientation matrix with the inverse of the 
	   orthogonalisation matrix results in the combined rotation
	   matrix */

	invert_3x3(*b, &b_inv);
	mat_mul_3x3(om, b_inv, u);

}

/* om_from_reccell_angles2(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void om_from_reccell_angles2(struct cell_struct cell, 
struct setangles_struct setangels, 
struct matrix_3x3 *om)
{
	struct matrix_3x3 b;	    /* Orthogonalization matrix and it's inverse */
	struct matrix_3x3 umatrix,  rotX, rotY, rotZ;  /* Combined Rotation matrices */

	/* First generate the orthogonalisation matrix from the reciprocal
	   cell parameters */

	B_from_reccell(cell, &b);

	rotmatX_from_angle(setangels.x, &rotX);
	rotmatY_from_angle(setangels.y, &rotY);
	rotmatZ_from_angle(setangels.z, &rotZ);

	mat_mul_3x3(rotZ, rotY, om);
	mat_mul_3x3(*om, rotX, &umatrix);
	mat_mul_3x3(umatrix, b, om);

}

/* om_from_B_and_angels(): Calculate the orientation matrix from the	*
 *                           reciprocal cell and the setting angles. 	*/

void om_from_B_and_angels(struct matrix_3x3 b, 
struct setangles_struct setangels, 
struct matrix_3x3 *om)

{
	struct matrix_3x3 umatrix,  rotX, rotY, rotZ;  /* Combined Rotation matrices */
	rotmatX_from_angle(setangels.x, &rotX);
	rotmatY_from_angle(setangels.y, &rotY);
	rotmatZ_from_angle(setangels.z, &rotZ);

	mat_mul_3x3(rotZ, rotY, om);
	mat_mul_3x3(*om, rotX, &umatrix);
	mat_mul_3x3(umatrix, b, om);


}

/* reduce_one_vector(): Routine for reducing one of the three lattice vectors; 
		       this routine is called by reduce_cell; input are the two 
		       reference vectors and the vector to reduce; output is the
		       changed reduced vector and an integer if the vector was 
		       changed						*/

int reduce_one_vector(struct vector a, struct vector b, struct vector *c)

{
        static const float EPS = 10e-05;
 	struct vector c_per,  c_par; /* perpendicular and parallel part of
						  work vector (c*)		*/

	struct vector n_vec;    /* normal to non-work vector plane (per to a* b* */
	struct vector I, m;	    /* I is perpendicular to the plane a*,n;
			       m st perpendicular to the plane b*, n;
				       length is 1				    */

	float p, q;		    /* parts of c in a* and b* direction */
	float p_new, q_new;	    /* New values of p and q */
	float tmp_float;

	/* Calculating n_vec */

	/* Getting direction */
	cross_prodf(a, b, &n_vec);
	/* Getting length */
	tmp_float=vec_lenf(n_vec);
	/* Normalizing */
	n_vec.x=n_vec.x/tmp_float;
	n_vec.y=n_vec.y/tmp_float;
	n_vec.z=n_vec.z/tmp_float;

	/* Getting c_per */
	/* Getting length */
	tmp_float=scalar_prodf(*c, n_vec);
	/* Calculation c_per */
	c_per.x=n_vec.x*tmp_float;
	c_per.y=n_vec.y*tmp_float;
	c_per.z=n_vec.z*tmp_float;

	/* Calculating c_par from c_per */
	c_par.x=c->x-c_per.x;
	c_par.y=c->y-c_per.y;
	c_par.z=c->z-c_per.z;

	/*Calculating p and q */
	/* Calculating I, perpendicular to the n-b plane */

	cross_prodf(n_vec, b, &I);
	/*Normalizing*/
	tmp_float=vec_lenf(I);
	I.x=I.x/tmp_float;
	I.y=I.y/tmp_float;
	I.z=I.z/tmp_float;

	/* Calculating m, perpendicular to the n-a plane */

	cross_prodf(n_vec, a, &m);
	/*Normalizing*/
	tmp_float=vec_lenf(m);
	m.x=m.x/tmp_float;
	m.y=m.y/tmp_float;
	m.z=m.z/tmp_float;

	/*Calculating p and q */
	p=scalar_prodf(c_par, I)/scalar_prodf(a, I);
	q=scalar_prodf(c_par, m)/scalar_prodf(b, m); /* In contrast to Kim's paper !! scalar_prodf(b,I) */

	p_new=fmod(p, 1.0);
	if (p_new >  0.5) p_new = p_new - 1.0;
	if (p_new < -0.5) p_new = p_new + 1.0;

	q_new=fmod(q, 1.0);
	if(q_new >  0.5) q_new = q_new - 1.0;
	if(q_new < -0.5) q_new = q_new + 1.0;

/*	if (p != p_new || q != q_new) { */
       if (fabs(p-p_new) > EPS || fabs(q-q_new) > EPS) {
		c_par.x=p_new*a.x + q_new*b.x;
		c_par.y=p_new*a.y + q_new*b.y;
		c_par.z=p_new*a.z + q_new*b.z;

		c->x=c_par.x+c_per.x;
		c->y=c_par.y+c_per.y;
		c->z=c_par.z+c_per.z;

		return(1);
	}
	else {
		return(0);
	}
}


/* reduce_cell(): Routine for cell reduction. input om. output new om.
		  programmed after the publication of Sansoo Kim, J. Appl
		  Cryst. 1989						*/

void reduce_cell(struct matrix_3x3 *om)
{
  int changed = 1;
      /* following needed to break out of an occasional infinite loop with 
           gcc/g77 v 2.96 */
      int changed_1 = 0;
      int changed_2 = 0;
      int changed_3 = 0;
      int same_1 = 0;
      int same_2 = 0;
      int same_3 = 0;

	while (changed != 0) {
        same_1 = changed_1;
        same_2 = changed_2;
        same_3 = changed_3;
	changed=reduce_one_vector(om->r[0], om->r[1], &(om->r[2]));
	changed_1 = changed;
	changed=changed+reduce_one_vector(om->r[2], om->r[0], &(om->r[1]));
	changed_2 = changed;
	changed=changed+reduce_one_vector(om->r[1], om->r[2], &(om->r[0]));
	changed_3 = changed;
	if (same_1 == changed_1 && same_2 == changed_2 && same_3 == changed_3)changed = 0;
  }
  sort_cell_vectors(om);
}


/* sort_cell_vectors(): Sort the vectors of a cell after their length */
void sort_cell_vectors(struct matrix_3x3 *om)

{
  	int i, j;
	struct vector help;
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			if(vec_lenf(om->r[j]) < vec_lenf(om->r[j+1])) {
				help=om->r[j+1];
				om->r[j+1]=om->r[j];
				om->r[j]=help;
			}
		}
	}
}


/* cell_from_cell(): Transform direct cell into reciprocal cell and vice versa. */

void cell_from_cell(struct cell_struct i, struct cell_struct *o)

{
	float vol;	/* Volume of input cell */
	float abc;	/* product i.a*i.b*i.c */
	double tmp;

	abc=i.a*i.b*i.c;
	vol=abc*sqrt(1-cos(i.alpha)*cos(i.alpha)
	    -cos(i.beta)*cos(i.beta)
	    -cos(i.gamma)*cos(i.gamma)
	    +2*cos(i.alpha)*cos(i.beta)*cos(i.gamma));

	o->a=i.b*i.c*sin(i.alpha)/vol;
	o->b=i.a*i.c*sin(i.beta)/vol;
	o->c=i.a*i.b*sin(i.gamma)/vol;

	tmp=vol/(abc*sin(i.beta)*sin(i.gamma));
	if (tmp >=1.0) tmp=1.0;
	o->alpha=asin(tmp);
	tmp=vol/(abc*sin(i.alpha)*sin(i.gamma));
	if (tmp >=1.0) tmp=1.0;
	o->beta =asin(tmp);
	tmp=vol/(abc*sin(i.alpha)*sin(i.beta));
	if (tmp >=1.0) tmp=1.0;
	o->gamma=asin(tmp);
}

/* print_cell_from_om(): Prints cell parameters generated from OM */

void print_cell_from_om(struct matrix_3x3 om)

{
	double pi = M_PI; 
        struct cell_struct cell, cell_rec; /* Cell and reciprocal cell */
	struct setangles_struct setangels; /* Setting angles */
	reccell_from_om(om,&cell_rec);
	dircell_from_reccell(cell_rec, &cell);
	print_cell(cell);
	setangles_from_reccell(cell_rec, om, &setangels);
	/*	printf("crystal rotx %f roty %f rotz %f\n",setangels.x*180/pi, setangels.y*180/pi, setangels.z*180/pi); */


}

/* print_cell(): Prints cell */

void print_cell(struct cell_struct cell)

{
	double pi = M_PI;

	/*	printf("unit cell %7.2f  %7.2f  %7.2f %6.2f %6.2f %6.2f\n",cell.a, cell.b, cell.c, cell.alpha*180/pi, cell.beta*180/ 
		pi, cell.gamma*180/pi); */


}
/* resort_cell_vectors(): Sort the vectors of a reciprocal cell according */
/* to the lengths of corresponding direct cell edges */
void resort_cell_vectors(struct matrix_3x3 *om, struct cell_struct *cell)

{
	int i, j;
	struct vector help;
	float a[3], rtemp;

	a[0] = cell->a;
	a[1] = cell->b;
	a[2] = cell->c;
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			if(a[j] > a[j+1]) {
				help=om->r[j+1];
				om->r[j+1]=om->r[j];
				om->r[j]=help;
				rtemp = a[j];
				a[j] = a[j+1];
				a[j+1] = rtemp;
			}
		}
	}
}
