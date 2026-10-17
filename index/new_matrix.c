#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/util.h"

/* new_matrix_(): This routine allows FORTRAN program to use DPS routines to
                  to take the regularized cell and calculate a new UB matrix
                  using the old setting angles. I'm not sure if we need this, 
		  but MOSFLM seems to like things that are very regular.
*/

/* this should either ve a void function - ie returning nothing
   or have a return statement in - I expect that the former is 
   appropriate GW 11.11.02 */

void new_matrix_(float *a, float fcell[6])  

{
  struct matrix_3x3 om;
  struct cell_struct cell;
  struct cell_struct reccell;
  struct setangles_struct setangels;

 float rtod = 180.0/3.141592654;
    om.r[0].x = a[0];
    om.r[0].y = a[1];
    om.r[0].z = a[2];
    om.r[1].x = a[3];
    om.r[1].y = a[4];
    om.r[1].z = a[5];
    om.r[2].x = a[6];
    om.r[2].y = a[7];
    om.r[2].z = a[8];

    cell.a = fcell[0];
    cell.b = fcell[1];
    cell.c = fcell[2];
    cell.alpha = fcell[3]/rtod;
    cell.beta  = fcell[4]/rtod;
    cell.gamma = fcell[5]/rtod;

    reccell_from_om(om, &reccell);
    setangles_from_reccell(reccell, om, &setangels);
    dircell_from_reccell(cell, &reccell);
    om_from_reccell_angles(reccell, setangels, &om);
    dircell_from_om(om,&cell);
    a[0] = om.r[0].x;
    a[1] = om.r[0].y;
    a[2] = om.r[0].z;
    a[3] = om.r[1].x;
    a[4] = om.r[1].y;
    a[5] = om.r[1].z;
    a[6] = om.r[2].x;
    a[7] = om.r[2].y;
    a[8] = om.r[2].z;

    fcell[0] = cell.a;
    fcell[1] = cell.b;
    fcell[2] = cell.c;
    fcell[3] = (cell.alpha)*rtod;
    fcell[4] = (cell.beta)*rtod;
    fcell[5] = (cell.gamma)*rtod;

}


