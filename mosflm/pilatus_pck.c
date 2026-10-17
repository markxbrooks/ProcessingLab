/****************************************************************************\
**    
$Id: pilatus_pck.c,v 1.1 2012/01/24 09:34:32 ccb Exp $                                                                      **
\****************************************************************************/
// Forced reading of a cbf image not using CBF_lib
// Only for CBF_BYTE_OFFSET or CBF_NONE compression
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int pilatus_pck_(char *image, char *odfile, int *nrec, int *iylen)
{

  FILE *cbf_filepointer;
  //struct stat statbuf;
  char *buf, *p;
  int size, *ierr, diff, nerr, binsize, cval , countx, county,
    twiceimgwidth;
  short *dp;
  double atime,btime;
	unsigned long fileLen;
union {	char *cp;
		short *sp;
		int *ip;
		} t;

// read entire image into buf and perform sanity checks.
 ctimer_(&atime);
 size = *nrec * *iylen;
 twiceimgwidth = *iylen * 2;
 memset (image,'\0',size);
  cbf_filepointer = fopen (odfile, "rb");
  p = 0;
  if (!cbf_filepointer)
  {
    fprintf (stderr, " Couldn't open the CBF file %s\n", odfile);
    *ierr = 1;
    return 1;
  } else {  

    fseek(cbf_filepointer, 0, SEEK_END);
    fileLen=ftell(cbf_filepointer);
    fseek(cbf_filepointer, 0, SEEK_SET);

    buf = (char *) malloc(fileLen+1);
    fread(buf, fileLen, 1, cbf_filepointer);
    
  }

// then:

// is it uncompressed? if so, just copy the image data
if ( (p=strstr(buf, "CBF_NONE")) )
	{
	if ( !(p=strstr(buf, "X-Binary-Size:")) )
		{
		printf("*** Bad header - cannot find X-Binary-Size\n");
		free(buf);
		return -1;
		}
	p += strlen("X-Binary-Size:");
	sscanf(p, "%d", &binsize);
	if ( !(p=strstr(p, "\f\032\004\325")) )
		{
		printf("*** Bad header - cannot find data start\n");
		free(buf);
		return -1;
		}
	p += strlen("\f\032\004\325");
	memcpy(image, p, binsize);
	free(buf);
	return 0;
	}

// find the start of the data
if ( !(p=strstr(buf, "CBF_BYTE_OFFSET")) )
	{
	printf("*** Bad header - cannot find CBF_BYTE_OFFSET\n");
	free(buf);
	return -1;
	}
if ( !(p=strstr(buf, "X-Binary-Size:")) )
	{
	printf("*** Bad header - cannot find X-Binary-Size\n");
	free(buf);
	return -1;
	}
p += strlen("X-Binary-Size:");
sscanf(p, "%d", &binsize);
if ( !(p=strstr(p, "\f\032\004\325")) )
	{
	printf("*** Bad header - cannot find data start\n");
	free(buf);
	return -1;
	}
p += strlen("\f\032\004\325");
 countx = 0;
 dp = (short *)image + size - *iylen;
// remember Mosflm uses 2 bytes per pixel
// size = 2*size; not used in current code
t.cp = p;
cval = 0;
 countx = 0;
 /*
 ** note that the somewhat odd mapping is to maintain consistency 
 ** with the original CBFlib code that is really very odd indeed
 */

 while (t.cp-p < binsize)		// correction here
   {
     /*
     ** If the file is corrupted with a lot of zeroes, 
     ** then (abs(*t.cp)<=127) is true too many times,
     ** and 'dp' tries to go beyond the buffer end.
     */
     if ((char *)dp < (image))
       break;
     if (abs(*t.cp)<=127) {
       diff = (int)*t.cp++;
     }
     else	                	// *t.cp == 0x80
       {
	 t.cp++;
	 if (*t.sp==0)			// flag to reset cval
	   {
	     cval = 0;
	     continue;
	   }
	 if (abs(*t.sp)<=32767){
	   
	   diff = (int)*t.sp++;
	 }
	 else				// *t.sp == 0x8000
	   {
	     t.sp++;
	     diff = *t.ip++;
	   }
       }
     cval = cval + diff;
     /* 
	super-packing into 16 bits from Andrew 
     */
     if(cval <= 32767) {
       *dp++ = cval;
     } else if (cval < 131072) {
       *dp++ = -((cval+8)/16) + 2045;
     } else if (cval < 3538000) {
       *dp++ = -((cval+64)/128) - 5124;
     } else {
       *dp++ = -32768;
     }
     countx++;
     if(countx >= *iylen){
       countx = 0;
       dp -= twiceimgwidth;
     }
     
   }

free(buf);
ctimer_(&btime);
printf("reading image took %f\n",btime-atime);
return 0;
}
