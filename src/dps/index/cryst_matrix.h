/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* cryst_matrix.h:  Headerfile for crystal_matrix  */

#include "vector_math.h"


#ifndef CRYSTL_MATRIX_H


struct cell_struct {
	double a,b,c;                 /* cell axis length in angstroem */
	double alpha,beta,gamma;    /* cell angles in radians        */
};

struct setangles_struct {
	double x,y,z;    /* Setting angles of rotations around x, y and z axex */
};

struct hkl_int {
	int h;
	int k;
	int l;
};

struct hkl_fract {
	float h;
	float k;
	float l;
};

void om_from_vectors(struct vector, struct vector, struct vector,
		     struct matrix_3x3 *);
void dircell_from_om(struct matrix_3x3, struct cell_struct *);
void reccell_from_om(struct matrix_3x3, struct cell_struct *);
void dircell_from_reccell(struct cell_struct, struct cell_struct *);
void setangles_from_reccell(struct cell_struct, struct matrix_3x3, 
		struct setangles_struct *);
void setangles_from_reccell_DENZO(struct cell_struct, struct matrix_3x3,
                struct setangles_struct *);
void hkl_2_xyz_single(struct matrix_3x3, struct hkl_int, struct vector *);
void xyz_2_hkl(struct matrix_3x3, struct vector *, int, struct hkl_int *,
	       struct hkl_fract *);	
void hkl_2_xyz(struct matrix_3x3, struct hkl_int *, int, struct vector *);	
void om_from_reccell_angles(struct cell_struct, 
			    struct setangles_struct, 
			    struct matrix_3x3 *);
void rotmatX_from_angle(float angle, struct matrix_3x3 *rot);
void rotmatY_from_angle(float angle, struct matrix_3x3 *rot);
void rotmatZ_from_angle(float angle, struct matrix_3x3 *rot);
void B_from_reccell(struct cell_struct cell, struct matrix_3x3 *b);
void B_and_U_from_om(struct matrix_3x3, struct matrix_3x3 *, struct matrix_3x3 *);
 
void om_from_reccell_angles2(struct cell_struct,
                            struct setangles_struct,
                            struct matrix_3x3 *);

void om_from_B_and_angels(struct matrix_3x3 b, 
			  struct setangles_struct setangels, 
			  struct matrix_3x3 *om);
void cell_from_cell(struct cell_struct , struct cell_struct *);
void reduce_cell(struct matrix_3x3 *);
void sort_cell_vectors(struct matrix_3x3 *);
void print_cell_from_om(struct matrix_3x3);
void print_cell(struct cell_struct);

#define CRYSTL_MATRIX_H
#endif
