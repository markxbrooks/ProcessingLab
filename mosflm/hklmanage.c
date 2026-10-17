/* hklmanage.c
 * originally by G.Winter
 * 30th September 2002
 * 
 * A group of three routines which are designed to allow storage of the
 * HKL values of the reflections between suibroutine calls, without all of 
 * the messing around with the image space - they are copied into a freshly
 * malloc'd array, which should greatly simplify matters. There are three 
 * routines:
 * 
 * call storehkl(ihklstr, ntot, status)
 * call retirevehkl(ihklstr, ntot, status)
 * call freehkl(status)
 * 
 * integer status, ntot
 * integer*2 ihklstr(mcolstr, nrefstr)
 * 
 * 
 * 
 * $Id: hklmanage.c,v 1.5 2004/06/02 14:49:33 harry Exp $
 */

#include <stdio.h>
#include <stdlib.h>

static short int *hklstore = NULL;
static int hklcount = 0;

/* storehkl
 * Store the hkl values from the array (input) into the newly allocated
 * array. 
 * 
 * 
 */

int storehkl_(short int *hklarray,
	      int *_ntot,
	      int *istat)
{
  int i, size, ntot;
  
  ntot = *_ntot;

  /* fprintf(stderr, "Storing %d reflections\n", *_ntot); */

  if(hklstore != NULL) {
    free(hklstore);
    hklstore = NULL;
  }

  /* allocate storage */

  hklstore = (short int *) malloc (sizeof(short int) * ntot * 6);

  /*  fprintf(stderr, "Mallocing some space %d bytes\n", */

  if(hklstore == NULL) {
    *istat = -1;
    return -1;
  }

  size = 6 * ntot;

  for(i = 0; i < size; i++) {
    hklstore[i] = hklarray[i];
  }

  hklcount = ntot;

  return 0;
}

/* retrievehkl
 * 
 * Retrieve the HKL values from the array (output) - the number
 * of reflections was stored.
 * 
 * 
 * 
 */ 

int retrievehkl_(short int *hklarray,
		 int *ntot,
		 int *istat)
{
  int i, size;

  if(hklstore == NULL) {
    *istat = -1;
    return -1;
  }

  if(hklcount == 0) {
    *istat = 0;
    return 0;
  }

  size = 6 * hklcount;

  for(i = 0; i < size; i++) {
    hklarray[i] = hklstore[i];
  }

  *ntot = hklcount;

  /*  fprintf(stderr, "getting %d reflections\n", hklcount); */

  return 0;
}

/* freehkl
 * 
 * Free up the allocated space which was used for the storage
 * of the reflections - this will make them impossible to retrieve.
 * 
 * 
 */ 

int freehkl_(int *istat)
{
  if(hklstore == NULL) {
    *istat = 0;
    return 0;
  }

  free(hklstore);
  hklstore = NULL;
  hklcount = 0;

  return 0;
}
      
    
    
  
