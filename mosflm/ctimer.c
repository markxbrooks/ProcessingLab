/* $Id: ctimer.c,v 1.4 2012/01/24 09:34:32 ccb Exp $ microsecond timer routine */
#if (defined _MSC_VER) || (defined __WIN32)
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#ifndef NULL
#define NULL (void *)0
#endif

void ctimer_(double *time_a) {
  
#if (defined _MSC_VER) || (defined __WIN32)
  struct timeb b_time;
  ftime(&b_time);
  *time_a = b_time.time + b_time.millitm/1000.0;
#else
  struct timeval a_time ;
  gettimeofday(&a_time,NULL) ;
  *time_a = a_time.tv_sec + a_time.tv_usec/1000000.0;
#endif
}

/* ****************************************************** */
