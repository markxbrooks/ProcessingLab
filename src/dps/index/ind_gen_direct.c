/* $Id: ind_gen_direct.c,v 1.32 2014/07/18 14:25:52 harry Exp $ */
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
 *                                                                     *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#if !defined vms && !defined(__APPLE__) && ! ( defined (__FreeBSD__) && defined (__STDC__) )
#include <malloc.h>
#endif 
#include "util.h"
#include "ind_gen_direct.h"
/* char *strcat( char *, const char *); */


/*------------------------------*
 * Global memory for FFT        *
 *------------------------------*/
static	int	last_fft_size	= 	0;
static	int	last_fft_dimen	= 	0;
static float 	*fft_help= NULL; /* Help array for storing arithm. tables etc for the FFT.
				    Also static, because this array needs only to be changed 
				    if max_int changes. */
static float 	*data  	= NULL;	 /*  Array for holding the binned data.
				     Also static because the size only 
				     changes if max_int changes. */

void  DeallocateFFT(void ) 
{
  if ( fft_help != NULL ){
    free(fft_help);
    fft_help = NULL;
  }
  if ( data != NULL ){
    free(data);
    data = NULL;
  }
  last_fft_dimen = 0;
  last_fft_size = 0;
	
}
int AllocateFFT( int n ) 
{

  /* 
   * Quick return if possible 
   */

  if ( last_fft_size == n ) {
    return 0;
  }
  /*
   * Only allocate memory if we are going to a bigger FFT
   */
  if ( n > last_fft_dimen ){
    size_t	size;

    last_fft_dimen = n;

    /* 
     * Free the old arrays 
     */
    if ( fft_help != NULL ){
      free(fft_help);
      fft_help = NULL;
    }
    if ( data != NULL ){
      free(data);
      data = NULL;
    }

    size = (2*last_fft_dimen + 15)*sizeof(float);
    fft_help = dps_malloc("get_fitting", "fft_help", size );

    size = last_fft_dimen*sizeof(float);
    data = dps_malloc("get_fitting", "data", size );

  }
  /* 
   * Call the initialization for the FFT routine 
   */
  if ( n != last_fft_size ){
    last_fft_size = n;
    rffti_(&last_fft_size, &fft_help[0]);
  }

  return last_fft_size;
}


/* get_projections_unsorted() Generates an array with the projections	*
 *		of all reflections in the xyz array onto the vector	*
 *		n_vector. Store the smallest projection in the proj[0]	*
 *		and the largest in proj[number_xyz-1].			*/

void get_projections_unsorted(struct vector xyz[], int number_xyz, 
			      struct vector n_vector, float *proj)
{
  int i;		/* Counter */
  int min_i, max_i;	/* Indices of the array elements with the maximum
			   and minimum values	*/
  float max, min;	/* Minimum and maximum values of the projections. */


  /* For initializing the variables for finding the minimum
     and maximum we have to handle the first coordinate outside
     the loop. */

  /* The projection is the scalar product of the two vectors	*/
  proj[0]=scalar_prodf(xyz[0],n_vector);

  /* Initialize variables for maximum, minimum and the corresponding
     number of the element in the proj array.	*/
  max=proj[0];
  min=proj[0];
  min_i=0;
  max_i=0;

  /* That's the main loop, looping through the reflections, generating
     the projections and getting the new maxima or minima if necessary    */
  for( i = 1; i < number_xyz; i++) {
    /* generate the projection and store it in the array */
    proj[i] = scalar_prodf(xyz[i],n_vector);
    /* Is the value bigger than the old max? */
    if(proj[i]>max) {
      max=proj[i];
      max_i=i;
      continue;
    }
    /* IS the value smaller than the old min? */
    if(proj[i]<min) {
      min=proj[i];
      min_i=i;
    }
  }
  /* Now we have to exchange the array element with the minimum
     value against the array element in proj[number_xyz-1] and
     the maximum value in proj[0]. */

  /* First save the values of proj[number_xyz-1] and proj[0] in the
     array elements of the maximum and minimum value. */
  proj[min_i]=proj[number_xyz-1];
  proj[max_i]=proj[0];
  /* Now store the minimum and maximum value in proj[number_xyz-1]
     and proj[0]. */
  proj[number_xyz-1]=min;
  proj[0]=max;
}


/* get_fitting(): This routine takes a list of projections, bins them and   *
 *		  does a FFT on the bins. After calculating the absolute    *
 *		  value of the complex result, the highest value between    *
 *		  the first minimum after 0 and a value corresponding to    *
 *		  maximum expected cell edge length is searched. The heights*
 *		  of this peak and the corresponding d_spacing value is	    *
 *		  returned.						    */

void get_fitting(float grid, float stretch, struct vector xyz[], int number_xyz, float d_max, 
		 float d_min, struct vector *n_vector, float proj[], float *d_experiment, 
		 float *n_fit, int peakmode, int npkfftmin, int nmultfft, int debug)

{
  /* on entry, 
     "debug" normally = 0. 
     For writing projections and FFTs to a file, = 2. 
     For normal debug statements, = 1. 
     For both, = 3 */

  int i,j;		    /* Counter */
  float extent;	    /* The range covered by the projections (max-min) */
  float shift;	    /* The shift necessary to bring the min value to zero */
  int max;		    /* Number of bins needed */
  float step;		    /* Stepsize from bin to bin */
  float fit;		    /* Heights of the highest peak in the FFT and temporary
			       variable during calculation. */
  float d_exp;	    /* Corresponding d-spacing for fit. */
  double mean_fft, rms_fft;
  double pi = M_PI;
  double	dtmp;
  unsigned int	fft_mask;
  int		fft_size;	/* FFT size -- Assumed to be a power of 2 */
  int np;

  int 	itemp, itest, overtones,width;
  int first_peak;
  float 	ftmp,overtonesfit,addwidth;
  float	df;
  int 	i0, i1, i2, i3;
  float 	w0, w1, w2, w3;
  int  harrydum;
  static int harryname = 0;
  char extension[6] = ".out";
  char charryname[13]  = "fft.";
  char rharryname[16]  = "revfft.";
  char pharryname[14] = "proj.";
  char pandrewname[18] = "unbinned.";
  char cincrement[5];   /*filename for putting the FFTs in...*/
  FILE     *fftout;
  FILE   *projout;
  FILE *revfftout;     	    /* Files for debugging output */
  FILE *projunbinned;
  
  /* set debug level to 2 if non zero */
  if (debug > 0) debug = 2;

  if(debug >= 2 && harryname <= N_TVECTOR) {
    harryname++;
    harrydum = harryname;
    cincrement[4] = '\0';
    if(harrydum >= 1000) {
      cincrement[0] = (harrydum/1000)+48;
      harrydum = harrydum - (1000*(harrydum/1000));
    }
    else {
      cincrement[0] = '0';
    }
    if(harrydum >= 100) {
      cincrement[1] = (harrydum/100)+48;
      harrydum = harrydum - (100*(harrydum/100));
    }
    else {
      cincrement[1] = '0';
    }
    if(harrydum >= 10) {
      cincrement[2] = (harrydum/10)+48;
      harrydum = harrydum - (10*(harrydum/10));
    }
    else {
      cincrement[2] = '0';
    }
    if(harrydum >= 0) {
      cincrement[3] = (harrydum)+48;
    }
    else {
      cincrement[3] = '0';
    }
    strcat(strcat(charryname,cincrement),extension);
    strcat(strcat(pharryname,cincrement),extension);
    strcat(strcat(rharryname,cincrement),extension);
    strcat(strcat(pandrewname,cincrement),extension);

    fftout = fopen(charryname,"w");

   //   outputting the projections is often of little use 
    projout = fopen(pharryname,"w");
    projunbinned = fopen(pandrewname,"w");

#ifdef WRITE_REVERSE
    revfftout = fopen(rharryname,"w");
#endif
   // outputting the projections is often of little use 

    if(harryname > N_TVECTOR)harryname=0;
  }

  
  /* First we have to bin the data. Therefore we have to know the maximum
     and minimum values of the projections. Also the projections list
     can be unsorted,  the first and last value have to correspond to the
     highest and lowest values.
  */
  extent=proj[0]-proj[number_xyz-1];

  /* We have to start at zero and therefore move the origin. This makes
     no difference for the FFT. 
  */
  shift=-proj[number_xyz-1];

  /* max is the number of bins we need. This is calculated in a way to
     assure that the introduced error by binning is not too large. The 
     value of grid gives the minimum number of bins per maximum expected
     cell edge length. The value of stretch can be used to get a higher 
     resolution for the returned d-spacing. With stretch != 1,  only
     part of the array will be used for the binned data. The rest is 
     padded with zeros. In the FFT this will result in a higher resolution
     in the returned values. Empirically we found,  that values of 5 for
     grid in 1 for stretch are sufficient. 
  */

  max=extent*d_max*grid*stretch;
  if(max>32768)max=32768;

  /* Most of the FFT routines work best with arrays of the length of the 
     power of two. Because we are binning anyway, we can modify stretch to
     match this requirement. 
  */

  fft_size = 1;
  while(fft_size<max){
    fft_size = fft_size*2;
  }

  /* ... modifying stretch to fit the array. */

  stretch=stretch*(fft_size/max);
  
  if(debug == 2) {
    printf("max %d fft_size %d \n", max, fft_size);
  }
  
  
  if(debug == 2) {
    printf("fft_size %d stretch %f\n", fft_size, stretch);
  }
  
  /* Calculating the stepsize between two consecutive bins. This value is
     needed to calculate the bin number of a projection. 
  */
  step = extent/((float)(fft_size-1)/stretch);
  
  if(debug == 2) {
    printf("max %d step %f \n", max, step);
  }
  
  /* In the first call or when the length of the bin array changes, we
     have to reallocate memory for the array and initialize the FFT routine. */

  if ( AllocateFFT( fft_size ) < 0 ){
    printf("Failed to set up FFT of size %d \n", fft_size );
    exit(1);
  }

  
  /* Zero the data array */
  for(i=0;i<fft_size;i++) {
    data[i]=0.0f;
  }
  
  /* Go through the list of reflections and calculate the bin number; then
     increment the number in the bin by one. 
  */
  /*
   * Modified GMcM to use cubic spline function weights:
   *     f(x) =   (1/6) ( 4 - x^2 + 3*|x|^3 )       |x| <= 1
   *              (1/6) ( 2 - |x| )^3               1 < |x| < 2
   *              0                                 |x| > 2
   * rather than delta function. Also we might as well drop
   * the (1/6).
   *
   * Array is made periodic by using fft_mask but this implementation
   * will only work for a fft_size that is a power of 2
   */
  /* original code was - reinstated 27.05.2010 to try to separate bins better
   for finding multiple lattices */
  /*   for(i=0;i<number_xyz;i++) {
     ++data[(int)((proj[i]+shift)/step)];
     }
  */

  fft_mask = fft_size - 1;
  for(i=0;i<number_xyz;i++) {
    ftmp = (proj[i]+shift)/step;
    i1 = (int)ftmp;
    df = ftmp - i1;
    i0 = ((i1 - 1) + fft_size) & fft_mask;
    i2 = (i1 + 1) & fft_mask;
    i3 = (i1 + 2) & fft_mask;
    w3 = df*df*df;
    w0 = 1.0 - 3.0*df + 3.0*df*df - w3;
    w1 = 4.0 - df*df*( 6.0 - 3.0*df );
    w2 = 1.0 + df*( 3.0 + df*( 3.0 - 3.0*df ));
    data[i0] += w0/6.0;
    data[i1] += w1/6.0;
    data[i2] += w2/6.0;
    data[i3] += w3/6.0;
  //harry trying to output the unbinned projections
  if(debug >= 2) {
    fprintf(projunbinned, "%d %f\n",i, ftmp);
    }
  }
  
  /*harry trying to output the projection */
    if(debug >= 2 && harryname <= N_TVECTOR) {
   for(i=1;i<(fft_size/2)-1;i++) {
     fprintf(projout, "%d %f\n",i, data[i]);
    }
   }
  
  /* Apply the FFT */
  rfftf_(&fft_size, &data[0], &fft_help[0]);
  /* Calculate the absolute value of the complex numbers. The complex numbers
     are stored as:
     c_real[1], c_imag[1], c_real[2], c_imag[2].....
     
     The absolute values are finally stored the following way:
     abs[1], abs[2], ....
     
     !! We do not take the square root,  because we only want to compare the
     size of the values. This saves a lot of time!! */
  /* Harry changed next line from for(i=1;i<fft_size;i=i+2) */



  data[0] = data[0]*data[0];			/* Zero frequency part is real */
  for(i=1; i < fft_size/2; i++ ) {
    data[i] = data[2*i-1]*data[2*i-1] + data[2*i]*data[2*i];
  }
  data[fft_size/2] = data[fft_size-1]*data[fft_size-1];	/* Nyquist frequency part is real */

  /* write FFT to file for inspection */
  if(debug >= 2 && harryname <= N_TVECTOR) {
      /*    printf("fft_size = %d\n",fft_size); */
    for(i=1;i<(fft_size/2)-1;i++) {
      fprintf(fftout, "%f %f\n",i/(fft_size*step),sqrt(data[i]));
    }
    /* this would print out a running average over 3 points, do not do this now
    fprintf(fftout, "\n"); 
    for(i=2;i<(fft_size/2)-2;i++) {
      fprintf(fftout, "%f %f\n",i/(fft_size*step),sqrt((data[i]+data[i-1]+data[i+1])/3.0f));
      } */
  }


  /* We need now to find the highest maximum corresponding to cell edges 
     between 0 and the maximum expected cell edge length. We also have to 
     avoid sitting on the slope of the maximum at zero. Therefore,  we
     start at 0, move along until we find a local minimum and look 
     in the area between the minimum and the maximum expected cell edge
     length for the highest peak.
  */


  /* First we need the value at 0 */


  
  /* The following loop finds the minimum, assign the value at the minimum to
     fit and sets i to  start with the next value */
  /*	i=d_max*(fft_size-1)*step*0.1; */
  /* first find the Mean value of the FFT */

  mean_fft = 0.0;
  for(i=0; i <= fft_size/2; i++){
    mean_fft = mean_fft + data[i];
  }
  mean_fft = mean_fft/((double)(fft_size/2+1));
  /* then find the standard deviation */
  rms_fft = 0.0;
  for(i=0; i <= fft_size/2; i++){
    // comment this line for RMS     dtmp = (data[i] - mean_fft);
    dtmp = data[i];
    rms_fft = rms_fft + dtmp*dtmp;
  }
  rms_fft = sqrtf(rms_fft/((double)(fft_size/2)));

  if(debug >= 2 && harryname <= N_TVECTOR) {
    // need to square root here to scale to the values plotted to the o/p files
    fprintf(fftout, "\n"); 
    fprintf(fftout, "%f %f\n",2.0f*i/(fft_size*step),sqrtf(mean_fft)); 
    fprintf(fftout, "0.0 %f\n",sqrtf(mean_fft));
    fprintf(fftout, "\n"); 
    fprintf(fftout, "%f %f\n",2.0f*i/(fft_size*step),sqrtf(rms_fft)); 
    fprintf(fftout, "0.0 %f\n",sqrtf(rms_fft)); 
  }
  /* these are the values for the squared FFT */

  if (peakmode==0) {
  /*
   * Looking for the first minimum
   */

  fit=data[0];
  i = 1;
  while( 2*i < fft_size-3 && data[i] < fit) {
    fit=data[i];
    i++;
  }
  if ( 2*i >= fft_size ){
    printf("Stepped over end of data array looking for first minimum %d\n", i);
    fflush(stdout);
    /*	exit(1); */
  }
	
  if(debug == 1 || debug == 3) {
    printf("minimum %d with %f\n", i-1, data[i-1]);
  }
  i++;
  itest = i;

  /* This loop, running from the minimum+1 to the bin corresponding to the
     maximum cell edge length,  finds the maximum in this area. It assigns
     the maximum value to fit and stores the corresponding bin number
     in d_exp. 
  */

  d_exp=i; /* If the bin after the minimum is already the global maximum. */
  fit=data[i];
  
  /*   
   * Added to avoid big origin peaks 
   */
  /*	while( 2*i < fft_size-3 && data[i] > mean_fft*mean_fft/(rms_fft) ){ */
  while( 2*i < fft_size-3 && (data[i]+data[i+1])/2 > mean_fft ){
    i++;
  }
  if ( 2*i >= fft_size && itest != i){
    printf("Stepped over end of data array while avoiding origin peak %d\n", i);
    fflush(stdout);
    /*		exit(1); */
  }
  fit=data[i];
  /* now find first big peak */
  first_peak = 0;
  while( i < d_max*(fft_size-1)*step) {
    // rms_fft = RMS value here
    //    if (data[i] > rms_fft){
      if(data[i] > fit) {
	fit=data[i];
	d_exp=i;
	first_peak = 1;
	if(debug == 1 || debug ==3) {
	  printf("%d ", i);
	}
      }
      //    }
    // rms_fft = RMS value here
      //    if (data[i] <= rms_fft && first_peak == 1){
      //      break;
      //    }
    i++;
  }
  if(debug == 1 || debug == 3) {
    printf("d_exp before %f \n", d_exp);
  }

#ifdef WRITE_REVERSE
  /* hrp 18.10.2006 using the maximum from this fft, reverse transform - start */
  if(debug >= 2 && harryname <= N_TVECTOR) {
    printf("shift is %f\n",shift);
    for( i=0; i < REVERSE_SIZE ;i++) { 
      // still needs the phase change angle...
      fprintf(revfftout, "%d %f\n",i,cos(((float)(i-50)*2.0*pi*d_exp/fft_size)+shift));
      proj[i] = cos(((float)(i-50)*2.0*pi*d_exp/fft_size)+shift);
    }
  }
  /* hrp 18.10.2006 using the maximum from this fft, reverse transform - end */
#endif
  /* 
   * Calculating the d_spacing from the bin number 
   */
  d_exp =((fft_size-1)*step)/d_exp;
  if(1.0/d_exp < d_min){
    fit = 0;
  }
  if(debug == 1 || debug == 3) {
    printf("max_search %f fit %f d_exp %f = 1/%f\n", 
            fabsf(d_max*(fft_size-1)*step), fit, d_exp,1/d_exp);
  }
  /* Preparing the return values */
  *d_experiment = d_exp;
  *n_fit=fit;

  }
  else {
  /* Call alternative FFT peakfinding here */
      np = fft_size/2;
      findlen_(&data[0],&np,&npkfftmin,&nmultfft,&debug,&fit,&d_exp);
       if (d_exp != 0.0) d_exp =((fft_size-1)*step)/d_exp;
      *d_experiment = d_exp;
      *n_fit=fit;
  }
  /* close the output files */

  if(debug >= 2) {
#ifdef WRITE_REVERSE
    fclose(revfftout);
#endif
    fclose(fftout);
    fclose(projout);
    fclose(projunbinned);
  }
}

/* ref_n_vector() Refines a direct space vector recursively and returns *
 *		  the refined vector.				        */

void ref_n_vector(int dumpref, int dumpfft[], struct vector xyz[],  int number_xyz, float d_max, 
                  float d_min, float proj[], float theta1_start, float theta2_start, 
                  float theta1_step, float theta2_step, int maxcyc, float convrg,
		  struct vector *n_vector,  float *n_fit, float *d_exp, float *theta1_ref, 
                  float *theta2_ref)

{
  float max_fit, work_fit;
  int i, j;
  int i_max, j_max;
  int peakmode, npkfft, nmfft;
  float work_exp, max_exp;
  float r;
  float theta1, theta2;
  float theta1_step_2, theta2_step_2;
  float theta1_max, theta2_max;
  struct vector work_vector, max_vector;
  float theta2_step_r;
  float pi=acosf(-1.0);
  int idebug;
  max_fit=0;
  theta1_step=theta1_step/2;
  theta2_step=theta2_step/2;

  if (dumpref != 0) 
    printf("step1 %f step2 %f theta1_start %f theta2_start %f cycle %d\n", 
            theta1_step, theta2_step, theta1_start,  theta2_start, maxcyc);

  for(i=-2;i<=2;i++) {
    theta1=theta1_start-i*theta1_step;
    work_vector.z=cosf(theta1);
    if(theta1<0) theta1=-theta1;
    r=sinf(theta1);

  if (dumpref != 0) printf("radius %f \n", r);

    if(r<0.000001) {
      theta2_step_r=0;
    }
    else {
      if(2.0*pi/5.0>theta2_step) {
	theta2_step_r=theta2_step/r;
      }
      else {
	theta2_step_r=(2.0*pi)/5.0;
      }
    }
    for(j=-2;j<=2;j++) {
      theta2=theta2_start+j*theta2_step_r;
      work_vector.x=sinf(theta2)*r;
      work_vector.y=cosf(theta2)*r;
      get_projections_unsorted(xyz, number_xyz, work_vector, proj);

      /* In this stage, use peakmode=0 to find first maximum in FFT */
      peakmode=0;
      npkfft=0;
      nmfft=0;
      idebug=dumpfft[1];
      get_fitting(5.0f, 1.0f,  xyz, number_xyz, d_max, d_min, &work_vector, 
         proj, &work_exp, &work_fit, peakmode, npkfft, nmfft, idebug); 

      /*      get_fitting(5.0f, 1.0f,  xyz, number_xyz, d_max, d_min, &work_vector, 
                          proj, &work_exp, &work_fit, dumpfft[0]); */


      if (dumpref != 0)
        printf("work_vector: %f %f %f fit: %f d_exp %f thetas: %f %f \n", 
                work_vector.x, work_vector.y, work_vector.z,work_fit, work_exp, theta1, theta2);

      if(work_fit>max_fit){
	cp_vectorf(work_vector, &max_vector);
	max_fit=work_fit;
	max_exp=work_exp;
	theta1_max=theta1;
	theta2_max=theta2;
	i_max=i;
	j_max=j;
      }

    }
  }

      if (dumpref != 0)
         printf("new_vector: %f %f %f fit: %f d_exp %f thetas: %f %f count: %d %d\n", 
         max_vector.x, max_vector.y, max_vector.z,
	 max_fit, max_exp, theta1_max, theta2_max, i_max, j_max);

      maxcyc-- ;

  if(((i_max<=-2 || i_max>=2) || (j_max<=-2 || j_max>=2)) && ( theta1_max >= 0.000001)) {
    theta1_step=theta1_step*1.5;
    theta2_step=theta2_step*1.5;
  }
  /* hrp 13.11.2006 0.1 degrees ~= 0.002 radians, sufficient for convergence */  
  if(((theta1_step>convrg) || (theta2_step>convrg)) && (maxcyc>0)){
    ref_n_vector(dumpref, dumpfft, xyz,  number_xyz, d_max, d_min, proj, theta1_max, theta2_max, 
    theta1_step, theta2_step, maxcyc, convrg, n_vector, n_fit, d_exp, theta1_ref, theta2_ref);
  }
  else {
    cp_vectorf(max_vector, n_vector);
    *d_exp=max_exp;
    *n_fit=max_fit;
    *theta1_ref=theta1_max;
    *theta2_ref=theta2_max;
    return;
  }
}
/* compare(): compares two floating point numbers to determine the bigger;  *
 * needed for qsort used in gen_direct_lattice_vectors()                    *
 *                                                                          */
int compare(const void *a, const void *b )
{
  return *(float *)b - * (float *)a;
}
/* gen_direct_lattice_vectors(): Generates direct lattice vectors of unity  *
 *	length,  calculates the projections of the reciprocal lattice points*
 *	onto this direction and uses FFT to find a periodicity in the	    *
 *	projections. Then find the most propable directions,  refine	    *
 *	the directions and return a list of direct space vectors with the   *
 *	heights of the highest peak in the FFT. The return value is the	    *
 *	number of found direct space vectors.				    */

int  gen_direct_lattice_vectors(int max_tvectors, int numvectref, 
                                int maxcyc, int inicyc, float convrg,
                                int numormat, float stepscal,
                                int numrlp, int numsptdps, int npkfftmin, 
                                int nmultfft, int idumpth, float theta1d, float theta2d,
                                int dumpfft[], int dumpvec, int dumpref,
                                int dumpgen, float d_max, float d_min,
				struct vector xyz[], int *number_xyz,
				struct tvector t_vector[])

{

  int i, j, iii;				/* Counters	*/
  double pi = M_PI;

  float *proj;	/* Array for the projection on one direct space direction */
  float beamprojdata[REVERSE_SIZE];
  float theta1_area;			/* Area to be covered by theta1 */

  float theta1, theta2;		/* Actually sphere coordinates */
  float theta1end,theta2end;               /* Final theta1,theta2 value */
  float theta1_step, theta2_step;	/* Stepsize for the sperical coordinates */
  float r;				/* Radius for the calculation of theta2 */
  float theta2_step_r;		/* Modified theta2_step value for the circle
				   with radius r. */

  int i_no, j_no;			/* Number of values in the array for storing
					   the heights of the highest peaks in the
					   FFT of a direct lattice directions. */
  float *circle = NULL;			/* Pointer to this array */
  unsigned short *intcircle = NULL;			/* Pointer to this array */
  float mean;				/* Mean value of the values of the array */
  float max;				/* Maximum value of the values of the array */
  int points;				/* Filled points in the array */

  struct vector work_vector;		/* direct space vector direction */
  float d_exp;			/* d-spacing found for the projections
				   on work vector in the FFT. */
  float n_fit;			/* Heights of the highest peak in the
				   FFT (at d_exp). */

  float low;				/* Actual threshold for identifying maxima
					   in the array */
  int mean_factor;			/* Factor for modifying threshold */
  float old_low;			/* Low value of the round before. */
  float small,last_small,big,last_big;
  float minfit = 999999.9;
  int iminfit = 9999;
  int no_more_dots = 0;
  int k, l;				/* Counters */
  int yes;				/* Boolean if there is a higher maximum in the
					   neighborhood... */
  int size;
  int t_vectors;			/* Number of direct space vectors, index for the array */
  int *local_maxima = NULL;             /* store t_vectors X [i,j] pairs */
  float *local_peak = NULL;             /* store t_vectors maxima */
  int icut,ipoint,isign,icheck;
  float theta1_ref;
  float theta2_ref;
  int peakmode;                         /* Determines algorithm used to find first peak in FFT */
  int nfinref;

  float vec_scale[N_TVECTOR],sin_vector[4],sin_vector_max;
  struct vector newvec[N_TVECTOR];



  if (dumpgen != 0) printf(" GEN_DIRECT_LATTICE_VECTORS entered...\n");

#ifdef WRITE_2D_ARRAY
  FILE *out; 
  printf(" GEN_DIRECT_LATTICE_VECTORS entered...\n");
  out=fopen("circle_001.img","w");
#endif

  /* Do some output */
  printf(" Generating direct lattice vectors (This may take some time):\n");

  /* Allocate some memory for the projections list... */
  proj = dps_malloc("gen_direct_lattice_vectors", "proj", sizeof(float)*(*number_xyz+1));

  /* The direct space vectors are generated by applying sphere	    *
   * coordinates. The two angles are theta1 and theta2. For a more	    *
   * detailed description please have a look in the publication	    */

  /* Area to cover for theta1 pi/2 == semi-sphere */
  theta1_area=pi/2;

  /* Step size for theta1. 0.03 was found empirically to work. */
  /* hrp 24.10.2006 calculate this dynamically, but never less than 0.03 */
  theta1_step = (28.55f/d_max > 0.03f) ? 28.55f/d_max : 0.03f;

  /* Modify step size by stepscal, defaults to 1.0 */
  theta1_step = stepscal*theta1_step;

  /* Choose the same stepsize for theta2 */
  theta2_step=theta1_step;

  /* The heights of the highest peak in the FFT for one direction is	    *
   * stored in an array with the indices i and j. i_max and j_max are the  *
   * maximum values of i and j to cover theta1_area with the stepsize	    */
  i_no=(int)(theta1_area/theta1_step+2);
  j_no=(int)(2*pi/theta2_step+2);

  if (dumpgen != 0) {
     printf("Theta 1 and 2 step size %f\n",theta1_step);
     printf("max_i max_j %d %d\n", i_no, j_no);
     printf("number_xyz: %d d_max: %f\n",*number_xyz,d_max);
  }


  /* Allocate memory for the array */
  /*	 free(circle); harry removes this as it's implicated in crashes after indexing failure with a too-small cell specified*/

  if ( i_no*j_no > 0 ){
    circle = dps_malloc("gen_direct_lattice_vectors", "circle", i_no*j_no*sizeof(float));
    intcircle = dps_malloc("gen_direct_lattice_vectors", "intcircle", i_no*j_no*sizeof(unsigned short));
    local_maxima = dps_malloc("gen_direct_lattice_vectors", "local_maxima",2*max_tvectors*sizeof(int));
    local_peak = dps_malloc("gen_direct_lattice_vectors", "local_peak",max_tvectors*sizeof(float));
  } else {
    printf("ERROR: i_no %d j_no %d \n", i_no, j_no );
    exit(1);
  }

  /* Harry initializes "circle" here as it causes a crash on a second call*/
  for(i=0;i<=i_no*j_no-1;i++)
    {
      circle[i] = 0.0f;
      intcircle[i] = 0;
    }
  /* Start values for some statistics. These will be used in finding the  *
   * maxima in the array afterwards.					    */
  mean=0.0;
  max=0.0;
  points=0;

  /* Start value for theta1. theta 1 will run from 0 to theta1_area */
  theta1=-theta1_step;

  /* Some more output */
  printf ("|");
  for(i=0; i<i_no; i++) printf("-");
  printf("|\n ");

  theta1end = theta1_area - theta1_step;
  /* start of theta1 loop */
  for(i=0;theta1<=theta1end;i++) {
    theta1=theta1+theta1_step;

    printf(".");
    fflush(stdout);

    /*    if (dumpgen != 0) printf("theta1 %f i %d\n",theta1, i); */


    /* The z-coordinate of the direct lattice vector is only dependent
       on theta1 */
    work_vector.z=cosf(theta1);

    /* Calculating the radius of the circle, theta2 will be generated
       in. */
    r=sinf(theta1);
    /* If theta1 == 0, we have only one point for theta2. */
    if(r<=0.00001) {
      /* These values will result in exact one step in theta2 */
      theta2_step_r=2*pi;
      theta2=-3*pi;
    }
    else {
      /* calculating the stepsize for theta2. */
      theta2_step_r= theta2_step/r;
      theta2=-pi-theta2_step_r;
    }
    /* start of theta2 loop */
    /*    for(j=0;theta2<pi;j++){*/
    theta2end = pi - theta2_step_r;
    for(j=0;theta2<theta2end;j++){
      theta2=theta2+theta2_step_r;

      /*     if (dumpgen !=0) printf("theta2 %f j %d\n", theta2, j); */


      /* If requested, calculate fft at specific theta1, theta2 values */
      if (j==0 && i==0 && idumpth!=0) {
        work_vector.z=cosf(theta1d);
        r=sinf(theta1d);
        work_vector.x=sinf(theta2d)*r;
        work_vector.y=cosf(theta2d)*r;
        get_projections_unsorted(xyz, *number_xyz, work_vector, proj);
        peakmode=0;
        get_fitting(5.0f, 1.0f, xyz, *number_xyz, d_max, d_min, &work_vector, proj, &d_exp, &n_fit, 
                  peakmode, npkfftmin, nmultfft, dumpfft[0]);
        printf("special theta1 %f theta2 %f fit: %f d_exp %f work_vector: %f %f %f \n\n", 
                 theta1d, theta2d, n_fit, d_exp, work_vector.x, work_vector.y, work_vector.z);
        work_vector.z=cosf(theta1);
        r=sinf(theta1);
      }

      /* Calculating the two other components of the direct lattice
	 vector */
      work_vector.x=sinf(theta2)*r;
      work_vector.y=cosf(theta2)*r;

      /*if (dumpgen !=0) printf("work_vector.x %f work_vector.y %f\n",work_vector.x,work_vector.y); */


      /* Now we can calculate the projections of the reciprocal
	 lattice points on the direct lattice vector. */
      get_projections_unsorted(xyz, *number_xyz, work_vector, proj);
      /* Finding a periodicity in the projections using FFT */
      /* In this stage, use peakmode=0 to find the first maximum in FFT */
      peakmode=0;
      get_fitting(5.0f, 1.0f, xyz, *number_xyz, d_max, d_min, &work_vector, proj, &d_exp, &n_fit, 
                  peakmode, npkfftmin, nmultfft, dumpfft[0]);

      if (dumpgen != 0) {
	  printf("theta1 %f theta2 %f fit: %f d_exp %f work_vector: %f %f %f \n", 
                 theta1, theta2, n_fit, d_exp, work_vector.x, work_vector.y, work_vector.z);

        fflush(stdout);
	/* printf("xyz: %f %f %f number_xyz: %d d_max: %f\n",xyz->x,xyz->y,xyz->z,*number_xyz,d_max); */
	/*        fflush(stdout); */
      }      



      /* This should test WRITE_2D_ARRAY rather than dumpgen, so comment */
      if (dumpgen != 0) {
      /*	fprintf(out, "%f %f %f %d %d\n", sinf(theta2)*theta1, cosf(theta2)*theta1, sqrt(n_fit),i,j);*/
      }


      /* Store n_fit, the heights of the highest peaks in the FFT in
	 the array. */
      circle[i*j_no+j]=n_fit;
      intcircle[i*j_no+j]=(unsigned short)sqrt(n_fit);

      /* For the statistics of mean and maximum value */
      mean += n_fit;
      if(n_fit>max) max=n_fit;
      ++points;
    }
    /* end of theta2 loop */

    /* Padding of unused fields (Some platforms (OSF)) will have problems otherwise) */
    for(j=j+1;j<j_no;j++) circle[i*j_no+j]=0.0;
  }
  /* end of theta1 loop */
#ifdef WRITE_2D_ARRAY
  fwrite((const void *)intcircle,sizeof(intcircle[0]),i_no*j_no,out);
  printf("i_no = %d, j_no = %d\n",i_no,j_no);
  fclose(out);
#endif
  printf(".\n");

  /* Calculating the mean value */
  mean=mean/points; 

  /* At this point we have an array (circle) with the heights of the highest
     peak of the FFT for the calculated direct space vectors. Now we have
     to find the largest of them and refine their directions. To minimise
     the number of direct space directions to be refined, we are looking
     for local maxima. This is done by looking at the array points around
     a point that is above a certain threshold. If one of the surrounding 
     points is higher,  we do not refine the direction of the actual
     direction. Also the array is distorted and the above described 
     algorithm will not work for small theta1's, it seems to be working.   */
  /* Initializing some variables */
  t_vectors=0;
  mean_factor=1; 
  low=max;       

  if (dumpgen != 0) printf("mean %f max %f\n", mean, max);

  /* We decrease our threshold for finding maxima until we have at least
     max_tvectors vectors refined. */
  printf(" Refining  %d direct space vectors:\n", max_tvectors);

  printf("|");
  for(i=0;i<max_tvectors;i++) printf("-");
  printf("|\n ");
  /* need to initialize old_low _before_ the next loop in case we repeat this step!! */
  old_low = -0.1;
  for (k = 0; k < REVERSE_SIZE; k++) {
    beamprojdata[k] = 0.0;
  }
  for(i=0;i<max_tvectors;i++) {
    local_maxima[2*i] = 0;
    local_maxima[2*i+1] = 0;
    local_peak[i] = 0.0;
  }  
  if (dumpgen != 0){
    printf("List of %d raw vectors (before refinement)\n",max_tvectors);
  }
  while(t_vectors < max_tvectors) {
    /* We need the old low value (in the first round = max) to not refine
       the same directions in every round. */
    if (low > old_low && old_low > 0.0) break;
    old_low=low; 
    /* If this is the first cycle, need to increas old_low or it will reject the highest peak */
    if (mean_factor==1) old_low=old_low+0.1;

    /* The threshold low will be varied from max to mean until we find
       enough maxima.  Sometimes "max" below is reduced to a value less than "mean" which 
       results in infinite looping - a kludge is to break out from here if this happens. 
       This has had 100% success rate in the one case I've tried it on... (HRP, 25.09.2002) */
    if(mean>max)break; 
    low=((float)mean_factor*mean+max)/(2.0+(float)mean_factor); 

    if (dumpgen != 0) printf("NEW MEAN_FACTOR!!! %d low %f \n", mean_factor, low);

    minfit = 999999.9;
    iminfit = 9999;
    /* looping through the array - theta1 stepping */
    for(i=0;i<i_no;i++) {

    /* looping through the array - theta2 stepping */
      for(j=0;j<j_no;j++) {

	/* Is the value between our threshold (low) and the threshold 
	   of the last round (old_low). If yes, the check the 
	   neighborhood... */
	if(circle[i*j_no+j]>low && circle[i*j_no+j]<old_low) { 
	  /* Sorry for reusing max. Here we store the value of the
	     central point in it. */
	  max=circle[i*j_no+j];

	  if (dumpgen != 0) printf("peak found, value %f i %d j %d \n",max,i,j);
						  
	  /* If yes stays 0, we do not have a higher value in the
	     neighboorhood. */
	  yes=0;

	  /* Looking one step in the array in each direction */
	  for(k=i-1;k<=i+1;k++) {
	    for(l=j-1;l<=j+1;l++) {
	      if(k >= 0 && k<i_no && l >= 0 && l < j_no) {
					      
		/* Did we find a higher value? */
		if(circle[k*j_no+l]>max) {
		  /* code folded from here */
		    if (dumpgen != 0) printf("bigger value at i %d j %d value %f\n",k,l,circle[k*j_no+l]);
		  yes=1;
		  /* unfolding */
		}
	      }
	    }
	  }

	  /* we have found a local maximum - store it and move onto the next one */
	  if (! yes ) {
	    if (t_vectors < max_tvectors) {
	      local_maxima[(2*t_vectors)] = i;
	      local_maxima[(2*t_vectors)+1] = j;
	      local_peak[t_vectors] = circle[i*j_no+j];
	      /* is this the smallest local_peak so far? */
	      if(circle[i*j_no+j] < minfit) {
		minfit = circle[i*j_no+j];
		iminfit = t_vectors;
	      }
	      if (dumpgen != 0){
                  theta1=i*theta1_step;
		  r=sinf(i*theta1_step);
                  if(r<0.000001) {
                    theta2=-pi;
                  }
                  else {
                    /* calculating the stepsize for theta2. */
                    theta2_step_r= theta2_step/r;
                    theta2=-pi+j*theta2_step_r;
                  }
		printf("New vector %3d theta1, theta2 = %6.3f %6.3f fit = %6.0f\n",t_vectors,
		       theta1, theta2, circle[i*j_no+j]);
	      }
	      t_vectors++;
	    } else {
	    /* hrp 26.10.2006 - if we already have max_tvectors but have not finished the 
	       search, we find the smallest one currently stored and replace its contents
	       if this one is bigger */
	      	      t_vectors = max_tvectors;
	      /*---------------------------------*/
	      if (circle[i*j_no+j] > minfit){
	      if (dumpgen != 0){
                  theta1=i*theta1_step;
		  r=sinf(theta1);
                  if(r<0.000001) {
                    theta2=-pi;
                  }
                  else {
                    /* calculating the stepsize for theta2. */
                    theta2_step_r= theta2_step/r;
                    theta2=-pi+j*theta2_step_r;
                  }
		printf("New vector %3d theta1, theta2 = %6.3f %6.3f fit = %6.0f  replacing  original value %6.0f\n",iminfit,
		       theta1, theta2, circle[i*j_no+j],local_peak[iminfit]);
	      }
	      	local_maxima[2*iminfit] = i;
	      	local_maxima[2*iminfit+1] = j;
	      	local_peak[iminfit] = circle[i*j_no+j];
		/* find the next minimum in the local_maxima array */
	      	minfit = 999999.9;
	      	iminfit = 9999;
	      	for (k = 0; k < max_tvectors; k++) {
	      	  if(local_peak[k] < minfit) {
	      	    minfit = local_peak[k];
	      	    iminfit = k;
	      	  }
	      	}
	      }
	      /*---------------------------------*/
	    }
	  }
	}
	/* end looping through the array - theta2 stepping */
      }
      /* end looping through the array - theta1 stepping */
    }
    /* We now looped through the array. If we have to do a next round we
       need the mean_factor increased for a lower threshold. */
    mean_factor++; 
  }
  if (dumpgen != 0){
    printf("Edited list of %d raw vectors (before refinement)\n",max_tvectors);
    for (k=0; k < t_vectors; k++) {
      i = local_maxima[2*k];
      j = local_maxima[2*k+1];
      theta1=i*theta1_step;
      r=sinf(i*theta1_step);
      if(r==0) {
        theta2=-pi;
      }
      else {
        /* calculating the stepsize for theta2. */
        theta2_step_r= theta2_step/r;
        theta2=-pi+j*theta2_step_r;
      }

      printf("    Vector %3d theta1, theta2 = %6.3f %6.3f fit = %6.0f\n",k,
	     theta1, theta2,circle[i*j_no+j]);
    }
  }
  /* we should now have the top max_tvector peaks based on FFT peak values
     at the grid points of the search. Now step through them and 
     partially refine each one (only inicyc cycles). 
     Then we sort the list on FFT peak height and 
     do further refinement on only the top numvectref vectors  */

  for (k=0; k < t_vectors; k++) {
    i = local_maxima[2*k];
    j = local_maxima[2*k+1];
    /* theta1 can be calculated from i and theta1_step */
    theta1=i*theta1_step;
    
    /* Theta2 will be calculated from j, theta2_step and theta1 */
    if(theta1==0) {
      theta2=0;
    }
    else {
      /* printf("theta1 = %f and sinf(theta1) = %f\n",theta1,
	 sinf(theta1));*/
      theta2=j*theta2_step/sinf(theta1)-pi;
    }

    if (dumpgen != 0) printf("\nVector %d initial theta1 %6.3f theta2 %6.3f value %6.0f\n",
                              k, theta1, theta2,circle[i*j_no+j] );
    
    
    /* Refine the direction recursively... */
    ref_n_vector(dumpref, dumpfft, xyz,  *number_xyz, d_max, d_min, proj, theta1, theta2, theta1_step,
		 theta2_step, inicyc, convrg, &work_vector, &n_fit, &d_exp, &theta1_ref, &theta2_ref);
    
    if (dumpgen != 0) printf("Refined vector %d : theta1, theta2: %6.3f %6.3f vector: %8.5f %8.5f %8.5f fit: %7.1f d_exp %8.6f theta1ref %6.3f theta2ref %6.3f\n", k,theta1,theta2,work_vector.x,
	   work_vector.y, work_vector.z, n_fit, d_exp, theta1_ref, theta2_ref);

    
    /* We have to recalculate the d-spacing on the final
       vector (d_exp), because it got lost
       in the refinement recursion. For this purpose we
       have to recalculate the directions and do the FFT. */
    get_projections_unsorted(xyz, *number_xyz, work_vector, proj);
    peakmode=0;
    get_fitting(5.0f, 1.0f, xyz, *number_xyz, d_max,d_min, &work_vector, proj, &d_exp,
		&n_fit, peakmode, npkfftmin, nmultfft, dumpfft[1]);

#ifdef WRITE_REVERSE
      // hrp 21.10.2006 
      for (l = 0; l < REVERSE_SIZE; l++){
	beamprojdata[l] += proj[l];
      }
#endif
    /* Now we can store the vector, scaled by the inverse 
       d-spacing found in the FFT. */
    
    t_vector[k].vec.x=work_vector.x/d_exp;
    t_vector[k].vec.y=work_vector.y/d_exp;
    t_vector[k].vec.z=work_vector.z/d_exp;
    
    /* Store also the heights of the highest peak in the 
       FFT. We do the square-root now, which bring n_fit
       onto a scale that makes it comparable to the
       number of reflections we put in. */
    t_vector[k].fit=sqrtf(n_fit);
    
    /* Calculate also the length of the vector and store it */
    t_vector[k].len=vec_lenf(t_vector[k].vec);
    t_vector[k].theta1=theta1_ref;
    t_vector[k].theta2=theta2_ref;
    
    if (dumpgen != 0) printf("final vector %f %f %f fit: %f len: %f d_exp: %f\n", t_vector[k].vec.x,
	   t_vector[k].vec.y, t_vector[k].vec.z, t_vector[k].fit,
	   t_vector[k].len, d_exp);


    printf(".");
    fflush(stdout);

    /* finished looping over vectors for initial refinement */    
  }
  

  /* Now sort vectors */


  /* Debug output of vector list */
  if (dumpvec != 0) {
    printf("t_vectors after initial refinement\n");
    print_tvector_list(t_vectors-1,&t_vector[0]);
  }

  /* Sort vectors in order of decreasing FFT peak height */

  sort_tvectors_on_n_and_l(t_vectors-1, &t_vector[0]);

  /* Debug output of sorted vector list */
  if (dumpvec != 0) {
    printf("t_vectors after first sort and initial refinement \n");
    print_tvector_list(t_vectors-1,&t_vector[0]);
  }

  /* Merge t_vectors with same length, same or opposite direction and
 similar number of fitting reflections. For more detail look at the
 comments in ind_sort_merge.c					*/

  t_vectors = merge_tvectors_same_n(t_vectors-1, &t_vector[0]);


  /* Debug output of merged vector list */
  if (dumpvec != 0) {
    printf("t_vectors after merge and initial refinement\n");
    print_tvector_list(t_vectors-1,&t_vector[0]);
  }


  if (t_vectors<numvectref) numvectref=t_vectors;

  /* now test for orthogonality of the top numvectref images. If there is 
     insufficient orthogonality, use all vectors in final refinement */

  sin_vector_max = 0.0;
  for(i=0;i<numvectref;i++){
    vec_scale[i] = 1.0/vec_lenf(t_vector[i].vec);
    newvec[i].x = vec_scale[i]*t_vector[i].vec.x;
    newvec[i].y = vec_scale[i]*t_vector[i].vec.y;
    newvec[i].z = vec_scale[i]*t_vector[i].vec.z;
  }
  for(i=2;i<numvectref;i++) {
    for(j=0;j<=i-2;j++) {
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
	if (sin_vector[3]>sin_vector_max){
	  sin_vector_max = sin_vector[3];
	}
      }
    }
  }

  if (dumpvec != 0) {
    printf("\n sin_vector_max for %d (numvectref) vectors: %f\n",numvectref,sin_vector_max);
  }
  if (sin_vector_max < 0.8) {
      numvectref=t_vectors;
  }

  if (dumpvec != 0) {
    printf("\nnumvectref now %d \n",numvectref);
  }

  /* Now do additional refinement of the top numvectref vectors */


  /* Calculate appropriate step size (depends on how many initial cycles performed */

  for (k=0; k <inicyc; k++) theta1_step = 0.5*theta1_step;
  theta2_step=theta1_step;


  for (k=0; k < numvectref; k++) {

      theta1=t_vector[k].theta1;
      theta2=t_vector[k].theta2;
      n_fit=t_vector[k].fit;
    if (dumpgen != 0) printf("\nVector %d initial theta1 %6.3f theta2 %6.3f value %6.0f\n",k, theta1, theta2,n_fit );
    
    /* Refine the direction recursively... */
    ref_n_vector(dumpref, dumpfft, xyz,  *number_xyz, d_max, d_min, proj, theta1, theta2, theta1_step,
		 theta2_step, maxcyc-inicyc, convrg, &work_vector, &n_fit, &d_exp, &theta1_ref, &theta2_ref);
    
    if (dumpgen != 0) printf("Refined vector %d : theta1, theta2: %6.3f %6.3f vector: %8.5f %8.5f %8.5f fit: %7.1f d_exp %8.6f theta1ref %6.3f theta2ref %6.3f\n", k,theta1,theta2,work_vector.x,
	   work_vector.y, work_vector.z, n_fit, d_exp, theta1_ref, theta2_ref);

    
    /* We have to recalculate the d-spacing on the final
       vector (d_exp), because it got lost
       in the refinement recursion. For this purpose we
       have to recalculate the directions and do the FFT. */
    get_projections_unsorted(xyz, *number_xyz, work_vector, proj);
    /* MOSFLM uses peakmode=1 (findlen_) here; that path is rejecting all
       vectors (d_exp==0) when called from processing-lab AutoIndexer on the
       insulin oracle .spt even though peakmode=0 already found ~78Å cells.
       Keep peakmode=0 for the final FFT peak pick until findlen parity is fixed. */
    peakmode=0;
    get_fitting(5.0f, 1.0f, xyz, *number_xyz, d_max,d_min, &work_vector, proj, &d_exp,
		&n_fit, peakmode, npkfftmin, nmultfft, dumpfft[2]);

    /* Now we can store the vector, scaled by the inverse 
       d-spacing found in the FFT. */
    /* Reject vectors for which no satisfactory peak found */
    if (d_exp==0) {
      t_vector[k].fit=0;
    }
    else {
      t_vector[k].vec.x=work_vector.x/d_exp;
      t_vector[k].vec.y=work_vector.y/d_exp;
      t_vector[k].vec.z=work_vector.z/d_exp;
    
      /* Store also the heights of the highest peak in the 
        FFT. We do the square-root now, which bring n_fit
        onto a scale that makes it comparable to the
         number of reflections we put in, unless using robust peakmode in which case
         the sqrt has already been done. */
      if (peakmode==1) {
        t_vector[k].fit=n_fit;
      }
      else t_vector[k].fit=sqrtf(n_fit);
    
      /* Calculate also the length of the vector and store it */
      t_vector[k].len=vec_lenf(t_vector[k].vec);
      t_vector[k].theta1=theta1_ref;
      t_vector[k].theta2=theta2_ref;
    }
    if (dumpgen != 0) printf("final vector %f %f %f fit: %f len: %f d_exp: %f\n", t_vector[k].vec.x,
	   t_vector[k].vec.y, t_vector[k].vec.z, t_vector[k].fit,
	   t_vector[k].len, d_exp);

    /* finished looping over vectors for final refinement */    
  }

    /* Now sort vector list again to eliminate those for which no peak was found */

  sort_tvectors_on_n_and_l(numvectref-1, &t_vector[0]);

  /* Debug output of sorted vector list */
  if (dumpvec != 0) {
    printf("t_vectors after finding true repeats \n");
    print_tvector_list(numvectref-1,&t_vector[0]);
  }
  nfinref=0;
  for (k=0; k < numvectref; k++) {
      if (t_vector[k].fit > 0.0) nfinref++ ;
  }				      
  numvectref=nfinref;


  printf("\n");
	
  /* Free the memory from the array */
  if ( circle != NULL ){
    free(circle);
    circle = NULL;
  }



#ifdef WRITE_REVERSE
  FILE *beamprojfile;
  beamprojfile = fopen("beamproj.file","w");
  for (k = 0; k < REVERSE_SIZE; k++){
    fprintf( beamprojfile,"%d %f\n",k,beamprojdata[k]);
  }
  fclose (beamprojfile);
#endif

  if (dumpgen != 0) printf(" GEN_DIRECT_LATTICE_VECTORS left with %d t_vectors...\n", numvectref);

  return(numvectref-1);
  
}
