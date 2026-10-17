/*=======================================================================
  $Id: ind_choose.c,v 1.20 2014/09/01 11:50:51 harry Exp $
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
 *                          File: ind_choose.c                         *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "index.h"
#include "vector_math.h"
#include "ind_choose.h"
#include "ind_gen_direct.h"
#include "util.h"
#include "celldb.h"

// a and b must be normalized
void find_vector_by_angle(struct vector a, struct vector b, float angle_a, float angle_b, struct vector *out) {
  struct vector x100 = {1, 0, 0}, y010 = {0, 1, 0};
  struct vector x;
  cross_prodf(x100, a, &x);
  if (vec_lenf(x) < 0.01) {
    cross_prodf(y010, a, &x);
  }
  float xlen = vec_lenf(x);
  x.x /= xlen; x.y /= xlen; x.z /= xlen;
  struct vector y;
  cross_prodf(a, x, &y);
  float ylen = vec_lenf(y);
  y.x /= ylen; y.y /= ylen; y.z /= ylen;
  
  float r, best_angle_diff = 999;
  for (r = 0; r < 2 * M_PI; r += 2 * M_PI / 360) {
    struct vector c;
    c.x = a.x / tan(angle_a) + cos(r) * x.x + sin(r) * y.x;
    c.y = a.y / tan(angle_a) + cos(r) * x.y + sin(r) * y.y;
    c.z = a.z / tan(angle_a) + cos(r) * x.z + sin(r) * y.z;
    float clen = vec_lenf(c);
    c.x /= clen; c.y /= clen; c.z /= clen;
    
    float angle_diff = fabs(acos(scalar_prodf(c, b)) - angle_b);
    //    printf("%f %f, ", r, angle_diff);
    if (angle_diff < best_angle_diff) {
      best_angle_diff = angle_diff;
      *out = c;
    }
  }	  
}  

/* choose_non_lin(): Choose the first three non linear depended 	*
 *		direct lattice vectors.					*/
int choose_non_lin(	int dumpch,float sinthmin, float volmin, float hkldev, float hklmin,
                        int number_xyz, struct vector xyz[],
			int t_vectors, struct tvector t[], int max_om, float tarcell[],
			struct matrix_3x3 *om )
{
  int i,j,k,imax,jmax,kmax;		/* counter */
  float triplevec;      /* Triple vector product of unit vectors */
  float	volume;		/* volume of the cell */
  float	tmp;		/* temporary variable */
  int om_count = 0;   /* Found orientation matrices */
  int a,b,c;		/* Number of vectors for generating om after
			   sorting */
  float vec_scale[N_TVECTOR],sin_vector[4],sin_vector_max,dtor,rtod;
  struct vector newvec[N_TVECTOR];
  int tmp_indexed,max_indexed;
  float triplevecmin;
  struct hkl_int	*hkl;		/* Integer hkl, same index as xyz */
  struct hkl_fract *hkl_calc;	/* Fractional hkl, same index as xyz */
  struct matrix_3x3 junk;
  FILE *chooseout;
  hkl = dps_malloc("choose_best_om", "hkl", sizeof(struct hkl_int)*number_xyz);
  hkl_calc = dps_malloc("choose_best_om", "hkl_calc", sizeof(struct hkl_fract)*number_xyz);

  /* hrp 04-03-2000 start This finds the solution which has 
     the largest product of the sines of the three cell angles, 
     to allow a cutoff to around half this value to reject 
     really bad solutions. Also the maximum number of reflections
     indexed by any solution is found, and the really bad ones 
     rejected at this stage.
  */

  max_indexed = 0;
  sin_vector_max = 0.0;
  triplevecmin = 0.2;
  imax = jmax = kmax = 0;
  if (dumpch != 0) 
      printf("\n In choose_non_lin, sinthmin %f volmin %f hkldev %f hklmin %f triplevecmin %f t_vectors %d\n",
	     sinthmin,  volmin,  hkldev,  hklmin, triplevecmin, t_vectors);

  // newvec: normalized. vec_scale: rec. space length
  for(i=0;i<=t_vectors;i++){
    vec_scale[i] = 1.0/vec_lenf(t[i].vec);
    newvec[i].x = vec_scale[i]*t[i].vec.x;
    newvec[i].y = vec_scale[i]*t[i].vec.y;
    newvec[i].z = vec_scale[i]*t[i].vec.z;
  }
  for(i=2;i<=t_vectors;i++) {
    //    printf("i=%d", i);
    for(j=0;j<=i-2;j++) {
      //      printf("j=%d", j);
      for(k=j+1;k<=i-1;k++) {
	if (dumpch != 0) printf("i=%d j=%d k=%d\n", i, j, k);
	// i - a, j - b, k - c
	// then 0 -- gamma, 1 -- alpha, 2 -- beta
	sin_vector[0] = sin(acos((newvec[i].x*newvec[j].x)+
				 (newvec[i].y*newvec[j].y)+
				 (newvec[i].z*newvec[j].z)));

	sin_vector[1] = sin(acos((newvec[i].x*newvec[k].x)+
				 (newvec[i].y*newvec[k].y)+
				 (newvec[i].z*newvec[k].z)));

	sin_vector[2] = sin(acos((newvec[k].x*newvec[j].x)+
				 (newvec[k].y*newvec[j].y)+
				 (newvec[k].z*newvec[j].z)));

	sin_vector[3] = fabs(sin_vector[0]*sin_vector[1]*
			     sin_vector[2]);
        /* calculate triple vector product */
	triplevec = vec_volumef(newvec[i],newvec[j],newvec[k]);
        if (triplevec>triplevecmin) {
	  if (sin_vector[3]>sin_vector_max){
	    sin_vector_max = sin_vector[3];
	  }

	  om_from_vectors(t[i].vec,t[j].vec,t[k].vec,&om[0]);
	  
	  /* Generate fractional and integer (rounded) HKL */
	  xyz_2_hkl(om[0], xyz, number_xyz, hkl, hkl_calc);
	  /* Calculate the number of reflections indexable with this OM. */
	  /*	tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, MAX_DIFF_HKL); */
	  tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
	  if (dumpch != 0) {
	    printf(" number spots indexed %d/%d for vectors %d %d %d triplesin %f triplevec %f \n", 
		   tmp_indexed,number_xyz,i,j,k,sin_vector[3],triplevec);
	  }
	  if (dumpch != 0)
	    print_3x3(om[0]);
	  if(max_indexed<tmp_indexed){
	    max_indexed=tmp_indexed;
	    imax=i;
	    jmax=j;
	    kmax=k;
	  }
	}
      }
    }
  }
  /* Debug output */
  if (dumpch != 0) 
	  printf("max value of triple sine product %f  max number spots indexed %d for vectors %d %d %d\n", 
                sin_vector_max,max_indexed,imax,jmax,kmax);

  max_indexed = max_indexed*hklmin;
  if (max_indexed == 0) max_indexed = number_xyz * hklmin * 0.9; // FIXME: ad hoc criteria

  /* hrp 04-03-2000 end insert */
  if (dumpch != 0) 
	  printf("\n max_indexed now set to %d for tests\n",max_indexed);
  if (dumpch != 0) 
      printf("\n now testing triple sin and volume\n");

  struct vector tmp_a, tmp_b, tmp_c;
  struct matrix_3x3 orig_cell, new_cell;
  float cell_tolerance = 10; // angstrom
  float angle_tolerance = 5; // degree
  int has_prior = (tarcell[0] != 0);

#define N_VECCOMB 50
  short vec_comb[N_TVECTOR][N_VECCOMB] = {};

  if (has_prior) {
    abc_to_cell(tarcell, &orig_cell);

    if (dumpch!= 0) printf("\nPreliminary Analysis of Basis vectors.\n\n");
    for (i = 0; i <= t_vectors; i++) {
      if (dumpch!= 0) printf("Vector #%d len = %6.2f:\n", i, 1 / vec_scale[i]);
      
      int a, b, c, cnt = 0;
      struct vector tmp;
      for (a = 0; a <= 2; a++) { // fix sign
	for (b = -2; b <= 2; b++) { 
	  for (c = -2; c <= 2; c++) {
	    if (a == 0 && b < 0) continue;
	    if (a == 0 && b == 0 && c < 0) continue; // avoid redundancy
	    if (fabs(a) + fabs(b) + fabs(c) > 4) continue;
	    tmp.x = a * orig_cell.r[0].x + b * orig_cell.r[1].x + c * orig_cell.r[2].x;
	    tmp.y = a * orig_cell.r[0].y + b * orig_cell.r[1].y + c * orig_cell.r[2].y;
	    tmp.z = a * orig_cell.r[0].z + b * orig_cell.r[1].z + c * orig_cell.r[2].z;
	    float tmp_len = vec_lenf(tmp);
	    if (fabs(1 / vec_scale[i] - tmp_len) < cell_tolerance / 3 && cnt < N_VECCOMB) {
	      vec_comb[i][cnt++] = 100 * (a + 2) + 10 * (b + 2) + c + 2;
	      if (dumpch!= 0) printf("   close to (%+1d, %+1d, %+1d) len:%6.2f\n", a, b, c, tmp_len);
	    }
	  }
	}
      }
    }
    
  }

  int ntried = 0, nfit = 0;
  int cnt;

  int i_op, j_op;
  
  // TAKANORI modified loop condition
  if (has_prior) {
    printf(" Fitting candidate vectors into supplied cell:\n|");
    for(i=1;i<=t_vectors;i++) {
      printf("-");
    }
    printf("|\n ");
  }
  for(i=1;i<=t_vectors;i++) {
    if (has_prior) {
      printf(".");
      fflush(stdout);
    }
    for(j=0;j<=i-1;j++) {
      if (!has_prior) { // make loop equivalent to the original when no prior is available
	if (i == 1) continue;
	if (j == i - 1) continue;
      } else {
	for (i_op = 0; vec_comb[i][i_op] != 0; i_op++) {
	  for (j_op = 0; vec_comb[j][j_op] != 0; j_op++) {
	    if (dumpch != 0) printf("i=%d j=%d i_op=%d vec_comb=%d j_op=%d vec_comb=%d\n", i, j, i_op, vec_comb[i][i_op], j_op, vec_comb[j][j_op]);
	    int reindex[9];
	    reindex[0] = vec_comb[i][i_op] / 100 - 2;
	    reindex[1] = vec_comb[i][i_op] % 100 / 10 - 2;
	    reindex[2] = vec_comb[i][i_op] % 10 - 2;
	    reindex[3] = vec_comb[j][j_op] / 100 - 2;
	    reindex[4] = vec_comb[j][j_op] % 100 / 10 - 2;
	    reindex[5] = vec_comb[j][j_op] % 10 - 2;
	    cnt = ops_index[reindex[0] + 2][reindex[1] + 2][reindex[2] + 2][reindex[3] + 2][reindex[4] + 2][reindex[5] + 2];

	    for (;cnt != 0 && ops[cnt] != -1; cnt++) {
	      reindex[6] = ops[cnt] / 100 - 2;
	      reindex[7] = ops[cnt] % 100 / 10 - 2;
	      reindex[8] = ops[cnt] % 10 - 2;
		
	      if (dumpch != 0) printf("expand cell by: %d %d %d, %d %d %d, %d %d %d cnt=%d\n", 
                reindex[0], reindex[1], reindex[2], reindex[3], reindex[4], reindex[5], reindex[6], reindex[7], reindex[8], cnt);
	      new_cell.r[0].x = reindex[0] * orig_cell.r[0].x + reindex[1] * orig_cell.r[1].x + reindex[2] * orig_cell.r[2].x;
	      new_cell.r[0].y = reindex[0] * orig_cell.r[0].y + reindex[1] * orig_cell.r[1].y + reindex[2] * orig_cell.r[2].y;
	      new_cell.r[0].z = reindex[0] * orig_cell.r[0].z + reindex[1] * orig_cell.r[1].z + reindex[2] * orig_cell.r[2].z;
	      new_cell.r[1].x = reindex[3] * orig_cell.r[0].x + reindex[4] * orig_cell.r[1].x + reindex[5] * orig_cell.r[2].x;
	      new_cell.r[1].y = reindex[3] * orig_cell.r[0].y + reindex[4] * orig_cell.r[1].y + reindex[5] * orig_cell.r[2].y;
	      new_cell.r[1].z = reindex[3] * orig_cell.r[0].z + reindex[4] * orig_cell.r[1].z + reindex[5] * orig_cell.r[2].z;
	      new_cell.r[2].x = reindex[6] * orig_cell.r[0].x + reindex[7] * orig_cell.r[1].x + reindex[8] * orig_cell.r[2].x;
	      new_cell.r[2].y = reindex[6] * orig_cell.r[0].y + reindex[7] * orig_cell.r[1].y + reindex[8] * orig_cell.r[2].y;
	      new_cell.r[2].z = reindex[6] * orig_cell.r[0].z + reindex[7] * orig_cell.r[1].z + reindex[8] * orig_cell.r[2].z;

	      cell_to_abc(&new_cell, tarcell);
	      if (dumpch != 0) printf("got cell: %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", tarcell[0], tarcell[1], tarcell[2], tarcell[3], tarcell[4], tarcell[5]);		
	      
	      // TAKANORI: try to generate 3rd vector from tarcell
	      // TODO: make tolerance modifiable
		
	      // alpha(3) - between b(1) and c(2)
	      // beta(4)  - between a(0) and c(2)
	      // gamma(5) - between a(0) and b(1)
	      ntried++;
	    
	      float calc_angle = vec_angle(newvec[i], newvec[j]) / M_PI * 180.0;
	      if (fabs(calc_angle - tarcell[5]) > angle_tolerance) {
		//		printf(" angle tolerance violated (%f should be %f). rejecting %d %d\n", calc_angle, tarcell[5], i, j);
  		continue;
	      }
	      nfit++;

	      struct vector z;
	      find_vector_by_angle(newvec[i], newvec[j], tarcell[3] / 180 * M_PI, tarcell[4] / 180 * M_PI, &z);

	      z.x *= tarcell[2]; z.y *= tarcell[2]; z.z *= tarcell[2];
	      if (dumpch != 0) printf("a = vec #%d, b = vec #%d, c from prior. alpha (b-c): %3.0f beta(a-c): %3.0f, ",i, j, vec_angle(z, newvec[i]) / M_PI * 180, vec_angle(z, newvec[j]) / M_PI * 180);

	      struct vector basis[3];
	      basis[0] = t[i].vec; basis[1] = t[j].vec; basis[2] = z;
	      if (t[i].len > t[j].len) {
		if (t[i].len > tarcell[2]) {
		  a=0;
		  b = t[j].len >= tarcell[2] ? 1 : 2;
		  c = t[j].len > tarcell[2] ? 2 : 1;
		}
		else {
		  a=2;
		  b=0;
		  c=1;
		}
	      }
	      else {
		if (t[j].len > tarcell[2]) {
		  a=1;
		  b = t[i].len >= tarcell[2] ? 0 : 2;
		  c = t[i].len > tarcell[2] ? 2 : 0;
		}
		else {
		  a=2;
		  b=1;
		  c=0;
		}
	      }
		
	      om_from_vectors(basis[2],basis[1],basis[0],&om[om_count]);
	      xyz_2_hkl(om[om_count], xyz, number_xyz, hkl, hkl_calc);
	      tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
	      if (dumpch != 0) printf("%d/%d indexed. ", tmp_indexed, number_xyz);
		
	      if(tmp_indexed >= max_indexed){
		++om_count;
		if (dumpch != 0) printf("Keeping this");
		if (om_count > (max_om-1)) return(om_count-1);
	      }
	      if (dumpch != 0) printf("\n");
	    }
	  }
	}
	continue;
      }
	
      // When no prior cell is available. 
      for(k=j+1;k<=i-1;k++) {
	  
	sin_vector[0] = sin(acos((newvec[i].x*newvec[j].x)+
				 (newvec[i].y*newvec[j].y)+
				 (newvec[i].z*newvec[j].z)));
	  
	sin_vector[1] = sin(acos((newvec[i].x*newvec[k].x)+
				 (newvec[i].y*newvec[k].y)+
				 (newvec[i].z*newvec[k].z)));
	  
	sin_vector[2] = sin(acos((newvec[k].x*newvec[j].x)+
				 (newvec[k].y*newvec[j].y)+
				 (newvec[k].z*newvec[j].z)));
	  
	sin_vector[3] = fabs(sin_vector[0]*sin_vector[1]*
			     sin_vector[2]);
	if (dumpch != 0) printf("i,j,k %d %d %d triplesin %f \n",i,j,k,sin_vector[3]);
	if(sin_vector[3]>sinthmin*sin_vector_max){ 
	  volume = vec_volumef(t[i].vec,t[j].vec,t[k].vec);
	  tmp=t[i].len * t[j].len * t[k].len;
          if (dumpch != 0) printf("  volume %f tmp %f \n",volume,tmp);
	  /*	  if ( (volume > tmp*LIN_DEP_FACTOR) ) { */
	  
	  // make sure edge length is i >= j >= k
	  if ( (volume > tmp*volmin) ) {
	    if (t[i].len > t[j].len) {
	      if (t[i].len > t[k].len) {
		a=i;
		b = t[j].len >= t[k].len ? j : k;
		c = t[j].len > t[k].len ? k : j;
	      }
	      else {
		a=k;
		b=i;
		c=j;
	      }
	    }
	    else {
	      if (t[j].len > t[k].len) {
		a=j;
		b = t[i].len >= t[k].len ? i : k;
		c = t[i].len > t[k].len ? k : i;
	      }
	      else {
		a=k;
		b=j;
		c=i;
	      }
	    }
	    
	    om_from_vectors(t[c].vec,t[b].vec,t[a].vec,&om[om_count]);
	    /* Generate fractional and integer (rounded) HKL */
	    xyz_2_hkl(om[om_count], xyz, number_xyz, hkl, hkl_calc);
	    /* Calculate the number of reflections indexable with this OM. */
	    /*	    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, MAX_DIFF_HKL); */
	    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
	    
	    if (dumpch != 0) 
	      printf("number %d chosen %d %d %d vol %f tmp %f rel %f indexed %d max %d\n", 
		     om_count,c,b,a,volume,tmp,tmp/volume,tmp_indexed,max_indexed);
	    
	    if(tmp_indexed >= max_indexed){
	      ++om_count;
	      if (dumpch != 0) printf("stored this matrix\n");
	      if (om_count > (max_om-1)) {
		printf("\n");
		return(om_count-1);
	      }
	    }
	  }
	  /*	  else {
		  if (dumpch != 0) 
		  printf("vectors i,j,k %d %d %d rejected, on volume: actual %f minimum %f\n",
		  i,j,k,volume,tmp*volmin);
		  }*/
	  
	}
	/*else {
          if (dumpch != 0) 
	  printf("vectors i,j,k %d %d %d rejected, triplesine only %f \n", 
	  i,j,k,sin_vector[3]);
	  }*/
      }
    }
  }
  

  if (has_prior) printf("\n%d/%d combinations fitted\n", nfit, ntried);
  if (dumpch != 0) printf("CHOOSE_NON_LINEAR left with %d matrices ...\n",om_count);
  return(om_count);
}

/* choose_best_om(): Given a list of possible orientation matrices, this routine returns the	*
 *		     number of the 'best' OM. This is decided on the number of reflections	*
 *		     that can be indexed within a certain margin around inter (MAX_DIFF_HKL)	*
 *		     and on the size of the unit cell. In principle the smallest unit cell	*
 *		     with the highest number of reflections is choosen.				*
 * Prerequirements: The list of orientation matrices should start with a OM generated by the	*
 *		    three shortest vectors (See routine index() ).				*/

int choose_best_om(int dumpch, float hkldev, int num_om, struct matrix_3x3 om[], int number_xyz, struct vector xyz[], int tarcellflg)

{
  int i,dpserrflg;			/* Counter */
  int best_indexed;	/* Number of indexed reflections for the currently chosen OM */
  int tmp_indexed;	/* Temporarely variable, same usage as best_indexed. */
  float best_volume;	/* Volume of the currently choosen OM */
  float tmp_volume;	/* Temporary variable, same usage as best_volume. */
  struct hkl_int	*hkl;		/* Integer hkl, same index as xyz */
  struct hkl_fract *hkl_calc;	/* Fractional hkl, same index as xyz */
  int selected;		/* Finally selected best OM. */
  struct matrix_3x3 junk;

  if (dumpch != 0) 
    printf("\n Choose_best_om entered, number of matrices: %d\n",num_om);

  /* Allocate memory for the hkl arrays. */
  hkl = dps_malloc("choose_best_om", "hkl", sizeof(struct hkl_int)*number_xyz);
  hkl_calc = dps_malloc("choose_best_om", "hkl_calc", sizeof(struct hkl_fract)*number_xyz);

  /* Set starting values, so that we definitely reach the statements of the 'if' statement */
  best_indexed = 0;
  best_volume = 1.0/vec_volumef(om[0].r[0], om[0].r[1], om[0].r[2]);

  /* Go through all OM's and find the 'best' */

  for(i=0;i<num_om;i++) {
    /* Generate fractional and integer (rounded) HKL */
    xyz_2_hkl(om[i], xyz, number_xyz, hkl, hkl_calc);
    /* Calculate the number of reflections indexable with this OM. */
    /*    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, MAX_DIFF_HKL); */
    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
    /* Calculate the volume of this OM. */
    tmp_volume = 1.0/vec_volumef(om[i].r[0], om[i].r[1], om[i].r[2]);
    /* Decide whether this OM is 'better' than the best found until now */
    /* zeroth: has it actually indexed any reflections? */
    if (tmp_indexed<1) {
      dpserrflg = 1;
      /* er - this was return(1) which failed if there was no solution better than 0! 
	 Harry 02-10-2002 */
      return(0);
    }
    /* First: Is the number of indexable reflections higher? */
    if (dumpch != 0) {
    printf("OM %d: indexed: %d volume %f\n",i,tmp_indexed,tmp_volume);
    print_3x3(om[i]);
    junk = om[i];
    print_cell_from_om(junk);
    xyz_2_hkl(junk, xyz, number_xyz, hkl, hkl_calc);
    /*    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, MAX_DIFF_HKL); */
    tmp_indexed = num_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
    bin_ok_hkl(hkl_calc, number_xyz, hkl, hkldev);
    printf("Number of indexed reflections %d\n\n",tmp_indexed);
    }

    if (best_indexed < tmp_indexed) {
      /* If yes, is the volume not too much higher. This should only avoid	    *
       * finding a higher order cell ( face/body centered etc.) which would have  *
       * twice the volume							    */
      /*hrp 01-03-2000*/

      if (tarcellflg == 1 || best_volume*VOLUME_FACTOR > tmp_volume) { // TAKANORI: disabled this check

	/* Ok, OM is better... */
	best_indexed = tmp_indexed;
	selected = i;
	best_volume=tmp_volume;

      }
    }
  }
  if (dumpch != 0) printf("CHOOSE_BEST_OM left with %d selected...\n",selected);

  return(selected);
}
