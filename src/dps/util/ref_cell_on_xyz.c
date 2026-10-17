/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 *======================================================================*/

/* ind_transform.c:  Transform peak coordinates  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

struct matrix_3x3 detinvB[CELL_N_PAR] = {
	{{{1, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {1, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}},
	{{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {0, 1, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}}
};



/* ref_cell_on_xyz(): Refines the given unit cell in respect to the reciprocal
		      space coordinates. Used in index where we do not want to 
		      have the detector geometry information
*/

void ref_cell_on_xyz(int number_xyz, struct vector *xyz, struct matrix_3x3 *om)

{
	int i, j, k;
	float rhs[CELL_N_PAR];
	double a[CELL_N_PAR][CELL_N_PAR];
	struct matrix_3x3 B, U, invU, invB, detinvOM, detinvOM_i, detinvOM_j;
	struct vector deth, deth_i, deth_j;
	struct hkl_fract *hkl;
	struct hkl_int *hkl_calc;
	struct vector *w;
	int num_w;
	int cycle;

	hkl=dps_malloc("REF_CELL_ON_XYZ", "HKL", number_xyz*sizeof(struct hkl_fract));
	hkl_calc=dps_malloc("REF_CELL_ON_XYZ", "HKL_CALC", number_xyz*sizeof(struct hkl_int));
	w=dps_malloc("REF_CELL_ON_XYZ", "W", number_xyz*sizeof(struct vector));

	cycle=0;
	xyz_2_hkl(*om, xyz, number_xyz, hkl_calc, hkl);
	num_w=0;
	for(i=0;i<number_xyz;i++) {
		if ((fabs(hkl[i].h-hkl_calc[i].h) > CELL_REF_MAX_DIFF_HKL) ||
		    (fabs(hkl[i].k-hkl_calc[i].k) > CELL_REF_MAX_DIFF_HKL) ||
		    (fabs(hkl[i].l-hkl_calc[i].l) > CELL_REF_MAX_DIFF_HKL) ) {
			w[i].x=0.0;
			w[i].y=0.0;
			w[i].y=0.0;
		}
		else {
			w[i].x=1.0;
			w[i].y=1.0;
			w[i].z=1.0;
			++num_w;
		}
	}
	while( cycle < 10 ) {
		xyz_2_hkl(*om, xyz, number_xyz, hkl_calc, hkl);
		B_and_U_from_om(*om, &B, &U);
		invert_3x3(B, &invB);
		invert_3x3(U, &invU);
		printf("before:\n");
		printf("om\n");
		print_3x3(*om);
		print_cell_from_om(*om);
		printf("B:\n");
		print_3x3(B);
		printf("U:\n");
		print_3x3(U);
		printf("invB:\n");
		print_3x3(invB);
		printf("invU:\n");
		print_3x3(invU);
		printf("Number of Reflections: %d\n", num_w);
		for(i=0;i<CELL_N_PAR;i++) {
			rhs[i]=0.0;
			mat_mul_3x3(detinvB[i], invU, &detinvOM);
			for(k=0;k<number_xyz;k++) {
				mat_vec_mul_3x3(detinvOM, xyz[k], &deth);
				rhs[i]=rhs[i]+w[k].x*deth.x*(hkl[k].h-hkl_calc[k].h)
				    +w[k].y*deth.y*(hkl[k].k-hkl_calc[k].k)
				    +w[k].z*deth.z*(hkl[k].l-hkl_calc[k].l);
			}
			for(j=0;j<CELL_N_PAR;j++) {
				a[i][j]=0.0;
				mat_mul_3x3(detinvB[i], invU, &detinvOM_i);
				mat_mul_3x3(detinvB[j], invU, &detinvOM_j);
				for(k=0;k<number_xyz;k++) {
					mat_vec_mul_3x3(detinvOM_i, xyz[k], &deth_i);
					mat_vec_mul_3x3(detinvOM_j, xyz[k], &deth_j);

					a[i][j]=a[i][j]+w[k].x*deth_j.x*deth_i.x
					    +w[k].y*deth_j.y*deth_i.y
					    +w[k].z*deth_j.z*deth_i.z;
				}
				printf("% .5e ", a[i][j]);
			}
			printf("\n");
		}
		gaussj(&a[0][0],CELL_N_PAR,&rhs[0],1);
		printf("shifts after cycle %d:\n",cycle);
		for(i=0;i<CELL_N_PAR;i++)
			printf("%f ",rhs[i]);
		invB.r[0].x=invB.r[0].x-rhs[0];
		invB.r[1].x=invB.r[1].x-rhs[1];
		invB.r[2].x=invB.r[2].x-rhs[2];
		invB.r[1].y=invB.r[1].y-rhs[3];
		invB.r[2].y=invB.r[2].y-rhs[4];
		invB.r[2].z=invB.r[2].z-rhs[5];
		invert_3x3(invB, &B);
		printf("after:\n");
		printf("B:\n");
		print_3x3(B);
		printf("invB:\n");
		print_3x3(invB);
		mat_mul_3x3(U, B, om);
		printf("om\n");
		print_3x3(*om);
		print_cell_from_om(*om);
		++cycle;
	}
}

struct matrix_3x3 detinvOM[OM_N_PAR] = {
	{{{1, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {1, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}},
	{{{0, 1, 0}, {0, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {0, 1, 0}}},
	{{{0, 0, 1}, {0, 0, 0}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 1}, {0, 0, 0}}},
	{{{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}}
};


void ref_om_on_xyz(int number_xyz, struct vector *xyz, struct matrix_3x3 *om)

{
	int i, j, k;
	float rhs[OM_N_PAR];
	double a[OM_N_PAR][OM_N_PAR];
	struct matrix_3x3 invOM;
	struct vector deth, deth_i, deth_j;
	struct hkl_fract *hkl;
	struct hkl_int *hkl_calc;
	struct vector *w;
	int num_w;
	int cycle;

	hkl=dps_malloc("REF_OM_ON_XYZ", "HKL", number_xyz*sizeof(struct hkl_fract));
	hkl_calc=dps_malloc("REF_OML_ON_XYZ", "HKL_CALC", number_xyz*sizeof(struct hkl_int));
	w=dps_malloc("REF_OM_ON_XYZ", "W", number_xyz*sizeof(struct vector));

	cycle=0;
	while( cycle < 10 ) {
		xyz_2_hkl(*om, xyz, number_xyz, hkl_calc, hkl);
		invert_3x3(*om, &invOM);
		printf("before:\n");
		printf("om\n");
		print_3x3(*om);
		printf("om_inv\n");
		print_3x3(invOM);
		num_w=0;
		for(i=0;i<number_xyz;i++) {
			if ((fabs(hkl[i].h-hkl_calc[i].h) > CELL_REF_MAX_DIFF_HKL) ||
			    (fabs(hkl[i].k-hkl_calc[i].k) > CELL_REF_MAX_DIFF_HKL) ||
			    (fabs(hkl[i].l-hkl_calc[i].l) > CELL_REF_MAX_DIFF_HKL) ) {
				w[i].x=0.0;
				w[i].y=0.0;
				w[i].y=0.0;
			}
			else {
				w[i].x=1.0;
				w[i].y=1.0;
				w[i].z=1.0;
				++num_w;
			}
		}
		printf("Number of Reflections: %d\n", num_w);
		for(i=0;i<OM_N_PAR;i++) {
			rhs[i]=0.0;
			for(k=0;k<number_xyz;k++) {
				mat_vec_mul_3x3(detinvOM[i], xyz[k], &deth);
				rhs[i]=rhs[i]+w[k].x*deth.x*(hkl[k].h-hkl_calc[k].h)
				    +w[k].y*deth.y*(hkl[k].k-hkl_calc[k].k)
				    +w[k].z*deth.z*(hkl[k].l-hkl_calc[k].l);
			}
			for(j=0;j<OM_N_PAR;j++) {
				a[i][j]=0.0;
				for(k=0;k<number_xyz;k++) {
					mat_vec_mul_3x3(detinvOM[i], xyz[k], &deth_i);
					mat_vec_mul_3x3(detinvOM[j], xyz[k], &deth_j);

					a[i][j]=a[i][j]+w[k].x*deth_j.x*deth_i.x
					    +w[k].y*deth_j.y*deth_i.y
					    +w[k].z*deth_j.z*deth_i.z;
				}
				printf("% .5e ", a[i][j]);
			}
			printf("\n");
		}
		gaussj(&a[0][0],OM_N_PAR,&rhs[0],1);
		printf("shifts after cycle %d:\n",cycle);
		for(i=0;i<OM_N_PAR;i++)
			printf("%f ",rhs[i]);
		invOM.r[0].x=invOM.r[0].x-rhs[0];
		invOM.r[1].x=invOM.r[1].x-rhs[1];
		invOM.r[2].x=invOM.r[2].x-rhs[2];
		invOM.r[0].y=invOM.r[0].y-rhs[3];
		invOM.r[1].y=invOM.r[1].y-rhs[4];
		invOM.r[2].y=invOM.r[2].y-rhs[5];
		invOM.r[0].z=invOM.r[0].z-rhs[6];
		invOM.r[1].z=invOM.r[1].z-rhs[7];
		invOM.r[2].z=invOM.r[2].z-rhs[8];
		invert_3x3(invOM, om);
		printf("after:\n");
		printf("om\n");
		print_3x3(*om);
		print_cell_from_om(*om);
		printf("om_inv\n");
		print_3x3(invOM);
		++cycle;
	}
	free(hkl);
	free(hkl_calc);
	free(w);
}
