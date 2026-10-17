#include <math.h>

#include "fxt.h"
#include "fxtaux.h"


void 
fht_spectrum(double *f, int ldn)
//
//   power_spectrum made with fht
//
//   thanks to Erwin Giling for pointing out 
//   that my code was wrong ! 
{
	fht(f,ldn);

	const long n=(1<<ldn);
	const long k=n/2;
	for (long i=1,j=n-1; i<k; i++,j--)
	{
		double a,b;

		a = f[i]+f[j];  
		// real part 
		    b = f[i]-f[j];  
		// imag part

		    f[i] = sqrt(a*a+b*b)*0.5;
		//	f[j] = atan2(b,a);  
		// if you want phases 
	}

	f[0]=sqrt(f[0]*f[0]+f[n/2]*f[n/2]);  
	// *0.5 ???
}
/* ============= end FHT_SPECTRUM ============ */


void 
fft_spectrum(double *f, int ldn)
//
//   power_spectrum made with fft
//
{
	split_radix_real_complex_fft(f,ldn);
	//    fht_real_complex_fft(f,ldn);

	const long n=(1<<ldn);
	for (long i=1,j=n-1; i<n/2; i++,j--)
	{
		double a,b;

		a = f[i];
		b = f[j];

		f[i] = sqrt(a*a+b*b);
		//	f[j] = atan2(b,a);  
		// if you want phases
	}

	f[0]=sqrt(f[0]*f[0]+f[n/2]*f[n/2]);
}
/* ============= end FFT_SPECTRUM ============ */

#if defined TEST

#include <stdio.h>
#define ldn 4
#define n   (1<<ldn)
int
main()
{
	long i;
	double f[n];

	for(i=0; i<n; ++i) f[i]=i;
	fft_spectrum(f,ldn);
	for(i=0; i<n/2; ++i)
		printf(" f[%ld]=%g \n",i,f[i]);
	printf("\n\n");


	for(i=0; i<n; ++i) f[i]=i;
	fht_spectrum(f,ldn);
	for(i=0; i<n/2; ++i)
		printf(" f[%ld]=%g \n",i,f[i]);
	printf("\n\n");

	return 0;
}

#endif /* defined TEST */
