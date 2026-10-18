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
 *                              dps_bravais                            *
 *                                                                     *
 *                        Written by Ingo Steller                      *
 *                                                                     *
 *                         File: bravais_goof.c                        *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "bravais.h"
#include "bravais_goof.h"

/* calc_tensor(): Calculate tensor from OM after International Tables A,
 *		  9.3.2,  p. 741 (1). 
 */

void calc_tensor(struct matrix_3x3 *om, Tensor *tensor)
{
	tensor->A=scalar_prodf(om->r[0], om->r[0]);
	tensor->B=scalar_prodf(om->r[1], om->r[1]);
	tensor->C=scalar_prodf(om->r[2], om->r[2]);
	tensor->D=scalar_prodf(om->r[1],om->r[2]);
	tensor->E=scalar_prodf(om->r[0],om->r[2]);
	tensor->F=scalar_prodf(om->r[0],om->r[1]);
}

/* goof_general(): Calculate a goof value for the general conditions for 
 *		   Type-I and Type-II cell,  9.3.2,  pp. 743,  (2a), (2b) for
 *		   TypeI; (4a), (4b), (4c) for Type-II.
 */

float goof_general(int type, Tensor *tensor)
{
	float goof;
	float help;

	/* The comments at the end of the lines are equivalent to the number of
       the formula in International tables Volume A,  9.3.2,  p. 743
    */
	switch(type) {
	case 1: 
		goof=0.0;
		help=tensor->A - tensor->B; /* (2a): a.a <= b.b */
		if(help > 0 ) {
			goof=goof+help;
		}
		help=tensor->B - tensor->C; /* (2a): b.b <= c.c */
		if(help > 0 ) {
			goof=goof+help;
		}
		help=2*fabs(tensor->D)-tensor->B; /* (2a): 2|b.c| <= b.b */
		if(help > 0) {
			goof=goof+help;
		};
		help=2*fabs(tensor->E)-tensor->A; /* (2a): 2|a.c| <= a.a */
		if(help > 0.0)  {
			goof=goof+help;
		};
		help=2*fabs(tensor->F)-tensor->A; /* (2a): 2|a.b| <= a.a */
		if(help > 0.0)  {
			goof=goof+help;
		};
		if(tensor->D < 0.0)  {  /* (2b): b.c >0 */
			goof=goof-tensor->D;
		}
		if(tensor->E < 0.0)  {  /* (2b): a.c >0 */
			goof=goof-tensor->E;
		}
		if(tensor->F < 0.0)  {  /* (2b): a.b >0 */
			goof=goof-tensor->F;
		}
		break;
	case 2: 
		goof=0.0;
		help=tensor->A - tensor->B; /* (2a): a.a <= b.b */
		if(help > 0 ) {
			goof=goof+help;
		}
		help=tensor->B - tensor->C; /* (2a): b.b <= c.c */
		if(help > 0 ) {
			goof=goof+help;
		}
		help=2*fabs(tensor->D)-tensor->B; /* (4a): 2|b.c| <= b.b */
		if(help > 0) {
			goof=goof+help;
		};
		help=2*fabs(tensor->E)-tensor->A; /* (4a): 2|a.c| <= a.a */
		if(help > 0.0)  {
			goof=goof+help;
		};
		help=2*fabs(tensor->F)-tensor->A; /* (4a): 2|a.b| <= a.a */
		if(help > 0.0)  {
			goof=goof+help;
		};
		help = 2*(fabs(tensor->D)+fabs(tensor->E)+fabs(tensor->F))
		    - (tensor->A + tensor->B); /* (4b): 2*(|b.c|+|a.c|+|a.b|) <= (|a.a+|b.b|) */
		if(help > 0) {
			goof=goof+help;
		}
		if(tensor->D > 0.0)  {  /* (4c): b.c < 0 */
			goof=goof+tensor->D;
		}
		if(tensor->E > 0.0)  {  /* (4c): a.c < 0 */
			goof=goof+tensor->E;
		}
		if(tensor->F > 0.0)  {  /* (4c): a.b < 0 */
			goof=goof+tensor->F;
		}
		break;
	default: 
		printf("Internal Error: Unknown type %d in goof_general!!\n", type);
		exit(1);
	}
	return(goof);
}

/* goof_group(): Calculate a goof value for the special conditions for the
 *		 groups in Table 9.3.1,  pp. 746
 */

float goof_group(int group, Tensor *tensor)
{
	float goof;
	float help;

	switch(group) {
	case 1: /* A=B=C */
		goof=fabs(tensor->A-tensor->B); /* A=B */
		goof=goof+fabs(tensor->B-tensor->C); /* B=C */
		break;
	case 2: /* A=B */
		goof=fabs(tensor->A-tensor->B); /* A=B */
		break;
	case 3: /* B=C */
		goof=fabs(tensor->B-tensor->C); /* B=C */
		break;
	case 4: 
		goof=0.0; /* No condition */
		break;
	default: 
		printf("Internal Error: Unknown group %d in goof_general!!\n", group);
		exit(1);
	}
	return(goof);
}

/* get_goof: Calculate a goof value for the goodness of fit of a given OM
 *	     to a lattice type. Lattice type is one off the 44 lattice types
 *	     in International Tables A,  9.3,  pp. 741 (see p 746 for a table
             of these transformations; note that #17, mC is incorrect in IT - 
             the entry should read 1,-1,0/-1,-1,0/-1,0,-1).
 */

float get_goof(int lattice, struct matrix_3x3 *om)

{
	float goof;
	Tensor tensor;
	float help;

	/* Generate Tensor from OM */
	calc_tensor(om, &tensor);
#ifdef BRAVAIS_GOOF_DEBUG
	printf("tensor: A %f B %f C %f D %f E %f F %f\n", tensor.A, tensor.B, tensor.C, tensor.D, tensor.E, tensor.F);
#endif
	switch(lattice) {
	case 1:	
	  /* cF */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(2.0*tensor.D-tensor.A); /* 2D = A */
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 2:	
	  /* hR */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		goof=goof+fabs(tensor.D-tensor.F); /* D = F */
		break;
	case 3:		
	  /* cP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 4:		
	  /* hR */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		goof=goof+fabs(tensor.D-tensor.F); /* D = F */
		break;
	case 5:		
	  /* cI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(-3.0*tensor.D-tensor.A); /* -3D = A */
		goof=goof+fabs(-3.0*tensor.E-tensor.A); /* -3E = A */
		goof=goof+fabs(-3.0*tensor.F-tensor.A); /* -3F = A */
		break;
	case 6:		
	  /* tI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		/* nonono! goof=goof+fabs(tensor.E-tensor.F);  E = F */
		break;
	case 7:		
	  /* tI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		goof=goof+fabs(tensor.E-tensor.F); /* E = F */
		break;
	case 8:		
	  /* oI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(1, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		break;
	case 9:		
	  /* hR */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(2.0*tensor.D-tensor.A); /* 2D = A */
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 10:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		break;
	case 11:	
	  /* tP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 12:	
	  /* hP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(-2.0*tensor.F-tensor.A); /* -2F = A */
		break;
	case 13:	
	  /* oC */
		goof=goof_general(2, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+goof_group(2, &tensor);
		break;
	case 14:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		/* no! goof=goof+fabs(tensor.E-tensor.F);  E = F */
		break;
	case 15:	
	  /* tI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(-2.0*tensor.D-tensor.A); /* -2D = A */
		goof=goof+fabs(-2.0*tensor.E-tensor.A); /* -2E = A */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 16:	
	  /* oF */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		goof=goof+fabs(tensor.D-tensor.E); /* D = E */
		/* no! goof=goof+fabs(tensor.E-tensor.F);  E = F */
		break;
	case 17:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(2, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		break;
	case 18:	
	  /* tI */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(4.0*tensor.D-tensor.A); /* 4D = A */
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 19:	
	  /* oI */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 20:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(tensor.E-tensor.F); /* E = F */
		break;
	case 21:	
	  /* tP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 22:	
	  /* hP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(-2.0*tensor.D-tensor.B); /* -2D = B */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 23:	
	  /* oC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 24:	
	  /* hR */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		goof=goof+fabs(-3.0*tensor.E-tensor.A); /* -3E = A */
		goof=goof+fabs(-3.0*tensor.F-tensor.A); /* -3F = A */
		break;
	case 25:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(3, &tensor);
		goof=goof+fabs(tensor.E-tensor.F); /* E = F */
		break;
	case 26:	
	  /* oF */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(4.0*tensor.D-tensor.A); /* 4D = A */
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 27:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 28:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(2.0*tensor.E-tensor.A); /* 2E = A */
		goof=goof+fabs(tensor.F-2*tensor.D); /* F = 2D */
		break;
	case 29:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.E-2*tensor.D); /* E = 2D */
		goof=goof+fabs(2.0*tensor.F-tensor.A); /* 2F = A */
		break;
	case 30:	
	  /* mC */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(2.0*tensor.D-tensor.B); /* 2D = B */
		goof=goof+fabs(2.0*tensor.E-tensor.F); /* 2E = F */
		break;
	case 31:	
	  /* aP */
		goof=goof_general(1, &tensor);
		goof=goof+goof_group(4, &tensor);
		break;
	case 32:	
	  /* oP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 33:	
	  /* mP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 34:	
	  /* mP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		break;
	case 35:	
	  /* mP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 36:	
	  /* oC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(-2.0*tensor.E-tensor.A); /* -2E = A */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 37:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(-2.0*tensor.E-tensor.A); /* -2E = A */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 38:	
	  /* oC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.D);  /* D=0 */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(-2.0*tensor.F-tensor.A); /* -2F = A */
		break;
	case 39:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(-2.0*tensor.F-tensor.A); /* -2F = A */
		break;
	case 40:	
	  /* oC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(-2.0*tensor.D-tensor.B); /* -2D = B */
		goof=goof+fabs(tensor.E);  /* E=0 */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 41:	
	  /* mC */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(-2.0*tensor.D-tensor.B); /* -2D = B */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 42:	
	  /* oI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(-2.0*tensor.D-tensor.B); /* -2D = B */
		goof=goof+fabs(-2.0*tensor.E-tensor.A); /* -2E = A */
		goof=goof+fabs(tensor.F);  /* F=0 */
		break;
	case 43:	
	  /* mI */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		goof=goof+fabs(2*fabs(tensor.D+tensor.E+tensor.F)
		    -(tensor.A+tensor.B));	/* 2|D+E+F| = A+B */
		goof=goof+fabs(fabs(2*tensor.D+tensor.F)
		    -tensor.B);	/* |2D+F| = B */
		break;
	case 44:	
	  /* aP */
		goof=goof_general(2, &tensor);
		goof=goof+goof_group(4, &tensor);
		break;
	default: 
		printf("Internal Error: Unknown lattice type %d in get_goof!!\n", lattice);
		exit(1);
	}
	return(goof);
}


