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
 *                        File: ind_sort_merge.c                       *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <malloc.h> */
#include "util.h"
#include "index.h"
#include "ind_sort_merge.h"

/* t_vector_exchange(): Exchange two t_vectors				*/

void t_vector_exchange(struct tvector *a, struct tvector *b)

{
	struct tvector tmp;

	/* a to tmp */
	cp_vectorf(a->vec,&(tmp.vec));
	tmp.fit=a->fit;
	tmp.len=a->len;
	tmp.theta1=a->theta1;
	tmp.theta2=a->theta2;
	/* b to a */
	cp_vectorf(b->vec,&(a->vec));
	a->fit=b->fit;
	a->len=b->len;
	a->theta1=b->theta1;
	a->theta2=b->theta2;
	/* tmp to b */
	cp_vectorf(tmp.vec,&(b->vec));
	b->fit=tmp.fit;
	b->len=tmp.len;
	b->theta1=tmp.theta1;
	b->theta2=tmp.theta2;
}


/* sort_tvectors_on_l(): Sort a given list of tvectors on theirlength	*/

void sort_tvectors_on_l(int t_vectors, struct tvector t[])

{
	int i,j;		/* counters */

	/* Sort by exchange						*/
	for(i=0;i<= t_vectors;i++) {
		for(j=0; j<= t_vectors-1; j++) {
			if (( t[j].len/t[j].fit > t[j+1].len/t[j+1].fit ) ) {
				t_vector_exchange(&t[j],&t[j+1]);
			}
		}
	}


}

/* sort_tvestors_on_n(): Sort a given list of tvectors on the number	*
 *			 of fitting reflections				*/

void sort_tvectors_on_n(int t_vectors, struct tvector t[])

{
	int i,j;		/* counters */

	/* Sort by exchange						*/
	for(i=0;i<= t_vectors;i++) {
		for(j=0; j<= t_vectors-1; j++) {
			if (( t[j].fit < t[j+1].fit ) ) {
				t_vector_exchange(&t[j],&t[j+1]);
			}
		}
	}


}


/* sort_tvestors_on_n_and_len(): Sort a given list of tvectors on the 	*
 *			number of fitting reflections and the on the 	*
 *			length of the vector.				*/

void sort_tvectors_on_n_and_l(int t_vectors, struct tvector t[])

{
	int i,j;		/* counters */

	/* Sort by exchange						*/
	for(i=0;i<= t_vectors;i++) {
		for(j=0; j<= t_vectors-1; j++) {
			if (( t[j].fit < t[j+1].fit ) ||
			    ( ( t[j].fit == t[j+1].fit ) &&
			    ( t[j].len > t[j+1].len  ) ) ) {
				t_vector_exchange(&t[j],&t[j+1]);
			}
		}
	}


}
/* Harry's new version!! 26.10.2006
 * merge_tvectors_same_n(): Merge a list of t_vectors. Merge	*
 *		vectors with same length and same or opposite direction *
 *		that have approx the same number of fitting reflections.	 *
 *		Input a list with t_vectors, output			*
 *		a list of merged t_vectors and the number		*
 *		of merged t_vectors					*
 *              Modified 23 Apr 07 to select the vector with the largest peak *
 *              rather than averaging the vectors AGWL */

int merge_tvectors_same_n(int t_vectors, struct tvector t[])

{
	int i,j;		/* counters */
	struct tvector opposite;
	int res_vectors;	/* Number of merged vectors and index of
				   res and res_count			*/
	int emptyslots;
	float match_pos,match_opp;
	float tolerance = 0.01; /* <1% difference is sufficient to say the vectors are the same */
        double twopi = 2.0*M_PI;

	emptyslots = 0;
	for(i = 0; i < t_vectors ; i++) {
	  if (t[i].fit > 0.1) {
	    for(j = i+1 ; j <= t_vectors; j++) {
	    /* check opposite direction, too */
	      opposite.vec.x = -t[j].vec.x;
	      opposite.vec.y = -t[j].vec.y;
	      opposite.vec.z = -t[j].vec.z;
	      opposite.len = t[j].len;
	      opposite.fit = t[j].fit;
	      opposite.theta1 = t[j].theta1;
	      opposite.theta2 = t[j].theta2;

	      /* I think that this is necessary because we're using floats, not doubles */

              if (scalar_prodf(t[i].vec,t[j].vec) == 0){
               match_pos = 0.0;
               match_opp = 0.0;
              } else {
	      match_pos = fabs((scalar_prodf(t[i].vec,t[j].vec) - (t[i].len*t[j].len))/scalar_prodf(t[i].vec,t[j].vec));
	      match_opp = fabs((scalar_prodf(t[i].vec,opposite.vec) - (t[i].len * opposite.len))/scalar_prodf(t[i].vec,opposite.vec));
              }
	      if ((match_pos <= tolerance || match_opp <= tolerance )&& t[j].len > 0.1) {
	      		  /* Try just choosing the larger peak, which will be the first one */
		  /* printf(" match between %d and %d\n",i,j); */
		emptyslots++;
		t[j].vec.x = 0.0;
		t[j].vec.y = 0.0;
		t[j].vec.z = 0.0;
		t[j].len = 0.0;
		t[j].fit = 0.0;
		t[j].theta1 = 0.0;
		t[j].theta2 = 0.0;
	      }
	      /* next line end of j loop */
	    }
	  }
	}
	res_vectors = t_vectors - emptyslots;
	/* now shuffle up the array */
	j = 0;
	for(i = 0 ; i <= res_vectors ; i++) {
	  if(t[i].fit == 0.0) {
	    while ( t[j].fit == 0 ) {
	      j++;
	    }
	    /* printf("setting vector %d to vector %d\n",i,j);   */
	    t[i].vec.x = t[j].vec.x;
	    t[i].vec.y = t[j].vec.y;
	    t[i].vec.z = t[j].vec.z;
	    t[i].len = t[j].len;
	    t[i].fit = t[j].fit;
	    t[i].theta1 = t[j].theta1;
	    t[i].theta2 = t[j].theta2;
	    t[j].vec.x = 0.0;
	    t[j].vec.y = 0.0;
	    t[j].vec.z = 0.0;
	    t[j].len = 0.0;
	    t[j].fit = 0.0;
	    t[j].theta1 = 0.0;
	    t[j].theta2 = 0.0;
	  }
	  j++;
	}
	for (i = res_vectors + 1; i <= t_vectors; i++){
	  t[i].vec.x = 0.0;
	  t[i].vec.y = 0.0;
	  t[i].vec.z = 0.0;
	  t[i].len = 0.0;
	  t[i].fit = 0.0;
	  t[i].theta1 = 0.0;
	  t[i].theta2 = 0.0;
	}

	/* return the new t_vector count */

	return(res_vectors);
}


