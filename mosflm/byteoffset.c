/* $Id: byteoffset.c,v 1.2 2018/09/19 13:50:20 andrew Exp $
   Routine to read byte offset style 4-byte image data from Pilatus
   specifically for Rigaku distributed detectors with d*Trek headers

   This should only read the image data, and should skip over header

   information
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


int byteoffset_(int *ierr, int *headersize, int *mode, int *nrec, int *iylen, int *invert, short image[64000000], char *imagefilename) {

  FILE *filein;

  clock_t a,b;


  int id, index, iindex,imagebytes;

  int i,ii,jjj,colrow, first,second,dirsta[2],dirend[2],dirinc[2],tmp_pixel;

  int maxpix, minpix, minx, miny, maxx, maxy ;
  int *localimage;
  char *header;
  char localfilename[480];
  /* the following declarations are for transferring data between 
     this C routine and the FORTRAN calling routine */
  minpix = 999999;
  maxpix = 0;
  /* start time for file read */
  a = clock ();
                      
    /* Read the file */
  i = 0;
  while ( imagefilename[i] != ' '){
    localfilename[i] = (char)imagefilename[i];
    i++;
  }
  localfilename[i] = '\0';
  filein = fopen (localfilename, "rb");
  if (!filein)
  {
    fprintf (stderr, " byteoffset couldn't open the image file >>>%s<<<\n", localfilename);
    *ierr = 1;
    return (*ierr);
  }
  
  header = malloc(sizeof(char) * *headersize);
  fread(header,sizeof(char),*headersize,filein);
  
  imagebytes = *nrec * *iylen;
  localimage = malloc(sizeof (int) * imagebytes);
  int* plocalimage = localimage;
  fread(localimage,sizeof(int),imagebytes,filein);

  //  free(localimage);
  //  return 0;
       /* invert the image in X */
  for ( i = 0; i < imagebytes; i++ ) {
    if (*invert==0) {
	index = i;
    }
    else {
	index = imagebytes - (2*(i/ *iylen)+1) * *iylen + i;
    }
    /*    if ( i < 20 ) {
	printf("mapping pixel %d in original to pixel %d \n",i,index);
    }
    if ( i > imagebytes - 20 ) {
	printf("mapping pixel %d in original to pixel %d \n",i,index);
	} */
    if (*mode < 10) {
      if (*localimage <= 32767) {
        image[index] = *localimage;
      } 
      else if (*localimage < 131072) {
        image[index] = -((*localimage+8)/16) + 2045;
      }
      else if (*localimage < 3538000) {
        image[index] = -((*localimage+64)/128) - 5124;
      }
      else {
        image[index] = -32766;
      }
    }
    else if (*mode==10) {
      if (*localimage == -1) {
	/* Preserve flagged pixels, -1 and -2 */
        image[index] = *localimage;
      }
      else if (*localimage > -1) {
	    image[index] = image[index] + *localimage;
      } 
    }
    else if (*mode==11) {
      if (*localimage == -1) {
	/* Preserve flagged pixels, -1 and -2 */
        image[index] = *localimage;
      }
      else if (*localimage > -1) {
	    image[index] = image[index] - *localimage;
      } 
    }
    localimage++;
  }
  free(header);
  free(plocalimage);
  fclose(filein);
  return 0;

}
