/* $Id: datetime.c,v 1.2 2012/04/18 09:24:15 ccb Exp $ */
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "../config.h"
#ifdef __WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif     
#define SIZE 256
     
int datetime_(char *dateandtime, char *username, char *compiler, 
	      char *cversion, char *mosflm_exec, char *builder,
	      char *build_date)
     {
       time_t curtime;
       struct tm *loctime;
#ifdef __WIN32
       unsigned long userbuffer = SIZE;
       char locname[SIZE+1];
#else
       char *locname[SIZE];
#endif
       memset(dateandtime,'\0',SIZE);
       memset(username,'\0',SIZE);
       memset(compiler,'\0',SIZE);
       memset(cversion,'\0',SIZE);
       memset(mosflm_exec,'\0',SIZE);
       memset(builder,'\0',SIZE);
       memset(build_date,'\0',SIZE);
     
#ifdef  __WIN32
       if (GetUserName(locname,&userbuffer) != 0 ) {
	 strcpy(username,locname);
#else
       if ((*locname = getlogin()) != NULL ) {
	 strcpy(username,*locname);
#endif
       } else {
	 strcpy(username,"non chieda il mio nome!");
       }
       /* Get the current time. */
       curtime = time (NULL);
     
       /* Convert it to local time representation. */
       loctime = localtime (&curtime);
       /* Print it out in a nice format. */
       strftime (dateandtime, SIZE, "%A, %B %d %Y at %H:%M", loctime);
       /* Do not change the next three lines! They are modified by the 
	  "build" script to contain the expected values */
/* Fortran compiler; */
 strncpy(compiler,FORTRAN_COMPILER,SIZE-1);
/* Compiler version; */
 strncpy(cversion,FORTRAN_VERSION,SIZE-1);
/* Executable type; */
 strncpy(mosflm_exec,FMOSFLM_EXEC,SIZE-1);
/* User name of builder */
 strncpy(builder,BUILD_USER,SIZE-1);
/* build date */
 strncpy(build_date,BUILD_DATE,SIZE-1);
      return 0;
     }
