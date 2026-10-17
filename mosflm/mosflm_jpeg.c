/* 
 *  $Id: mosflm_jpeg.c,v 1.74 2015/05/27 10:37:41 harry Exp $
 *

Specification for the image requests/responses
15th April 2002
G. Winter
 
1. All coordinates will be given in terms of the display, in which the origin
   is in the top lefthand corner of the image, and `x' is perpendicular to the
   lower edge of the image:

   +----->x
   |
   |
   |
   V
   y

2. Creation of the jpeg images should handle the possibility that the limits 
   may be a little back-to-front. For instance, the possibility that the x
   coordinate of the tlc might be greater than tat of the brc.

3. The section of the image should be specified in terms of the tlc and brc,
   and nothing else. This will make for simple tracing of the variables. Thus,
   in every subroutine/function, the input must be of the form of two arrays
   each of two integers. These will be in the pixel coordinates, rotated to
   match the screen orientation.

4. When referring to the x and y ordinates in an ordered pair, the order will
   be (x, y).

5. The jpeg writing routine should be:

   subroutine write_jpeg(destination,
			 nrec, 
			 iylen, 
			 image, 
			 factor, 
			 tlc, 
			 brc,
			 newline) 

   Where destination is the file descriptor of the socket, nrec and iylen are
   the dimensions of the image in memory, image is the pointer to the image
   array, factor is the desired `zoom' factor, tlc and brc are as discussed
   and newline is an integer to say whether to append a newline to the image
   or no (relevent for embedding the image in an XML document).

   In this subroutine, the relative orientation of the tlc and brc are 
   important, since this should still be in the display coordinates.	       
 */
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if (defined _MSC_VER) 
/* no need to do this for __MINGW_ (Mingw gcc) this is just about unique!*/
#define S_IRWXU S_IREAD|S_IWRITE|S_IEXEC
#else
#include <unistd.h>
#endif
#if (defined _MSC_VER) || (defined __MINGW_H) || (defined __WIN32)
#include <winsock.h>
#endif
#include "jconfig.h"
#include "jpeglib.h"
#include "jerror.h"
#include "mosflm_jpeg.h"

typedef struct triple{
  unsigned char r;
  unsigned char g;
  unsigned char b;
} triple;

/* TGGB added typdef greylevel to be used instead of triple for greyscale image */

typedef unsigned char greylevel ;

static int jpeg_quality=90;
static char imagename[80] = "";
static unsigned char base64convert[64];

/* structures used internally */

typedef struct mosflm_jpeg_mgr {
  struct jpeg_destination_mgr pub;
  unsigned char cmpbuffer[54000]; /* 54000 bytes at a time will be written here */
  unsigned char encbuffer[72000]; /* then base64 encoded to 72000 bytes to write */
  int socket;
} mosflm_jpeg_mgr;

/* from TkPhotoImageBlock definition from Tk libs */

typedef struct {
  unsigned char *pixelPtr;
  int width;
  int height;
  int pitch;
  int pixelSize;
  int offset[4];
}  TkPhotoImageBlock;

/* utility functions, for base64 encoding */
/* this should be called before any encoding is done */

int base64init() {
  int i;

  for(i = 0; i < 26; i++) {
    base64convert[i] = 65 + i;
    base64convert[i + 26] = 97 + i;
  }

  for(i = 0; i < 10; i++) {
    base64convert[i + 52] = 48 + i;
  }

  base64convert[62] = 43;
  base64convert[63] = 47;

  return 0;
}

int base64get(
	      int i
	      )
{
  return base64convert[i];
}

/*
** Translation Table as described in RFC1113
*/
static const char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

int base64enc (unsigned char in[3], unsigned char out[4], int len )
{
  out[0] = cb64[ in[0] >> 2 ];
  out[1] = cb64[ ((in[0] & 0x03) << 4) | ((in[1] & 0xf0) >> 4) ];
  out[2] = (unsigned char) (len > 1 ? cb64[ ((in[1] & 0x0f) << 2) | ((in[2] & 0xc0) >> 6) ] : '=');
  out[3] = (unsigned char) (len > 2 ? cb64[ in[2] & 0x3f ] : '=');
  return 0 ;
}

int base64encline(
		  unsigned char *in,
		  unsigned char *out,
		  int nbytes
		  )
{
  int i, n, remainder;

  n = nbytes / 3;
  remainder = nbytes % 3;

  for(i = 0; i < n; i++) {
    base64enc(&(in[3 * i]), &(out[4 * i]), 3);
  }

  if(remainder != 0) {
    base64enc(&(in[3 * n]), &(out[4 * n]), remainder);
  }

  return 0;
}

/* this routine is called to initiate the jpeg writing process */

void mosflm_init_destination(
			     j_compress_ptr cinfo
			     )
{
  struct mosflm_jpeg_mgr *dest = (struct mosflm_jpeg_mgr *) cinfo->dest;
  /* the storage here is statically allocated */
  dest->pub.next_output_byte = dest->cmpbuffer;
  memset(dest->cmpbuffer, 0, 54000);
  dest->pub.free_in_buffer = 54000;
}

/* this routine is called each time the buffer gets full */

boolean mosflm_empty_output_buffer(
				   j_compress_ptr cinfo
				   )
{
  struct mosflm_jpeg_mgr *dest = (struct mosflm_jpeg_mgr *) cinfo->dest;
  struct stat socket_or_file;
  int file_or_socket = 0;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(dest->socket,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }

  base64encline(dest->cmpbuffer, dest->encbuffer, 54000);
#if (defined __MINGW_H) || (defined __WIN32)
  if (file_or_socket == 0) {
    write(dest->socket, dest->encbuffer, 72000);
  } else {
    send(dest->socket, dest->encbuffer, 72000,0);
  }
#else
  write(dest->socket, dest->encbuffer, 72000);
#endif
  memset(dest->cmpbuffer, 0, 54000);
  memset(dest->encbuffer, 0, 72000);
  dest->pub.next_output_byte = dest->cmpbuffer;
  dest->pub.free_in_buffer = 54000;
  return TRUE;
}

boolean mosflm_empty_output_bin(
				j_compress_ptr cinfo
				)
{
  struct mosflm_jpeg_mgr *dest = (struct mosflm_jpeg_mgr *) cinfo->dest;
  struct stat socket_or_file;
  int file_or_socket = 0;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(dest->socket,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }
#if (defined __MINGW_H) || (defined __WIN32)
  if (file_or_socket == 0) {
    write(dest->socket, dest->cmpbuffer, 54000);
  } else {
    send(dest->socket, dest->encbuffer, 54000,0);
  }
#else
  write(dest->socket, dest->cmpbuffer, 54000);
#endif
  dest->pub.next_output_byte = dest->cmpbuffer;
  dest->pub.free_in_buffer = 54000;
  return TRUE;
}

/* and this routine gets called when the writing process is finished */

void mosflm_term_destination(
			     j_compress_ptr cinfo
			     )
{
  int i, n;
  struct mosflm_jpeg_mgr *dest = (struct mosflm_jpeg_mgr *) cinfo->dest;
  struct stat socket_or_file;
  int file_or_socket = 0;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(dest->socket,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }
  i = 54000 - dest->pub.free_in_buffer;
  n = 4 * (i / 3);
  if((4 * i) > (3 * n)) {
    n += 4;
  }

  base64encline(dest->cmpbuffer, dest->encbuffer, i);
#if (defined __MINGW_H) || (defined __WIN32)
  if (file_or_socket == 0) {
    write(dest->socket, dest->encbuffer, n);
  } else { 
    send(dest->socket, dest->encbuffer, n,0);
  }
#else
  write(dest->socket, dest->encbuffer, n);
#endif
}

void mosflm_term_destination_bin(
				 j_compress_ptr cinfo
				 )
{
  int i, n;
  struct mosflm_jpeg_mgr *dest = (struct mosflm_jpeg_mgr *) cinfo->dest;
  struct stat socket_or_file;
  int file_or_socket = 0;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(dest->socket,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }
  i = 54000 - dest->pub.free_in_buffer;
#if (defined __MINGW_H) || (defined __WIN32)
  if (file_or_socket == 0) {
    write(dest->socket, dest->cmpbuffer, i);
  } else {
    send(dest->socket, dest->cmpbuffer, i,0);
  }
#else
  write(dest->socket, dest->cmpbuffer, i);
#endif
}


void mosflm_write_ppm(
		      int sock,
		      int nx,
		      int ny,
		      triple *newimage,
		      int newline,
		      int colours,
		      int binary
		      )
{
  greylevel *image_byte, *temp_image_byte, *encbuffer;
  //  triple *image_byte3, *temp_image_byte3, *encbuffer3;
  char header_string[22];
  char twentyeightbytes[28];
  double mostime_a,mostime_b;
  /* 
     for PPM (colour), image_code = 'P3'(ascii)  or 'P6' (binary)
     for PGM (grey),   image_code = 'P2'   "     or 'P5'     "
     ascii not used!
     this code now allows coloured images
  */
  int nxy = nx * ny;
  int nxy43;
  int i,j,stringlength;
  struct stat socket_or_file;
  int file_or_socket = 0;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(sock,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }
  base64init();
  nxy43 = 4 * (nxy / 3);
  if((4 * nxy) > (3 * nxy43)) {
    nxy43 += 3;
  }
  if (colours == 1){
    /*
      encode the red channel for grey level
    */
    j = 5;  
    image_byte = (greylevel *) malloc (sizeof(greylevel) * nxy + 2);
    encbuffer = (greylevel *) malloc (sizeof(greylevel) * nxy43);
    temp_image_byte = image_byte;
    for (i = 0; i < nxy; i++){
      *temp_image_byte++ = newimage++->r;
    }
  } else {
    /* colour image */
    j = 6;
    nxy *=3;
    nxy43 *=3;
    image_byte = (greylevel *) malloc (sizeof(greylevel) * nxy + 2);
    encbuffer = (greylevel *) malloc (sizeof(greylevel) * nxy43);
    
  }
  stringlength = 22;
  i = snprintf(header_string,stringlength,"P%d %6d %6d 255 \n",j,nx,ny);
  if (binary == 0){
    base64encline((greylevel *)header_string,(greylevel *)twentyeightbytes,21);
#if (defined __MINGW_H) || (defined __WIN32)
    if (file_or_socket == 0) {
      write( sock, twentyeightbytes, 28*sizeof(char));
    } else {
      send(sock, twentyeightbytes, 28*sizeof(char), 0);
    }
#else
    write( sock, twentyeightbytes, 28*sizeof(char));
#endif
    base64encline(image_byte, encbuffer, nxy);
    i = 0;
    if (colours == 1) {
      while ( i < nxy43 ){
	j = nxy43 - i > 72000 ? 72000 : nxy43 - i;
#if (defined __MINGW_H) || (defined __WIN32)
	if (file_or_socket == 0) {
	  write(sock,&encbuffer[i],j*sizeof( greylevel));
	} else {
	  send(sock,&encbuffer[i],j*sizeof( greylevel),0);
	}
#else
	write(sock,&encbuffer[i],j*sizeof( greylevel));
#endif
	i += 72000;
      }
    } else {
      while ( i < nxy43 ){
	j = nxy43 - i > 72000 ? 72000 : nxy43 - i;
#if (defined __MINGW_H) || (defined __WIN32)
	if (file_or_socket == 0) {
	  write(sock,&newimage[i],j*sizeof( greylevel));
	} else {
	  send(sock,image_byte,j*sizeof( greylevel),0);
	}
#else
	write(sock,&newimage[i],j*sizeof( greylevel));
#endif
	i += 72000;
      }
    }
  } else {
    if (colours == 1) {
#if (defined __MINGW_H) || (defined __WIN32)
      if (file_or_socket == 0) {
	write( sock, header_string, 21*sizeof( char));
	write( sock, image_byte, nxy*sizeof( greylevel));
      } else {
	send( sock, header_string, 21*sizeof( char),0);
	send( sock, image_byte, nxy*sizeof( greylevel),0);
      }
#else
      write( sock, header_string, 21*sizeof( char));
      write( sock, image_byte, nxy*sizeof( greylevel));
#endif
    } else {
#if (defined __MINGW_H) || (defined __WIN32)
      if (file_or_socket == 0) {
	write( sock, header_string, 21*sizeof( char));
	write( sock, newimage, nxy*sizeof( greylevel));
      } else {
	send( sock, header_string, 21*sizeof( char),0);
	send( sock, image_byte, nxy*sizeof( greylevel),0);
      }
#else
      write( sock, header_string, 21*sizeof( char));
      write( sock, newimage, nxy*sizeof( greylevel));
#endif
    }
  }

  if(newline != 0) {
#if (defined __MINGW_H) || (defined __WIN32)
    if (file_or_socket == 0) {
      write( sock, "\n", 1);
    } else {
      send( sock, "\n", 1,0);
    }
#else
    write( sock, "\n", 1);
#endif
  }
  free(image_byte);
  free(encbuffer);
}

int mosflm_write_jpeg(
		      int sock,
		      int nx,
		      int ny,
		      triple *image,
		      int newline,
		      int colours,
		      int binary
		      )
{
  struct mosflm_jpeg_mgr dest;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  struct stat socket_or_file;
  int file_or_socket = 0;
  JSAMPROW row_pointer[1];
  int i,nxy;
  greylevel *image1,*tmpimage1;
  base64init();
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  cinfo.dest = (struct jpeg_destination_mgr *) &dest;
  /* trying to see if sock is a socket or file */
  memset(&socket_or_file,0,sizeof(socket_or_file));
  fstat(sock,&socket_or_file);
  switch (socket_or_file.st_mode & S_IFMT) {
  case S_IFREG: 
    /* it's a file */
    file_or_socket = 0;
    break;
  default: 
    /* it must be a socket since it isn't a file */
    file_or_socket = 1;
    break;
  }

  dest.socket = sock;

  dest.pub.init_destination = mosflm_init_destination;
  if(binary == 0) {
    dest.pub.empty_output_buffer = mosflm_empty_output_buffer;
    dest.pub.term_destination = mosflm_term_destination;
  } else {
    dest.pub.empty_output_buffer = mosflm_empty_output_bin;
    dest.pub.term_destination = mosflm_term_destination_bin;
  }
  cinfo.image_width = nx;
  cinfo.image_height = ny;
  /* 
     grey scale only wanted 
  */
  if (colours == 1) {
    image1 = (greylevel *) malloc (sizeof(greylevel) * nx * ny);
    memset(image1,0,sizeof(greylevel) * nx * ny);
    tmpimage1 = image1;
    for (i = 0; i < nx * ny; i++){
      *tmpimage1++ = image++->r;
    }
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;
  } else {
    
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;
  }
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, jpeg_quality, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  if (colours == 1) {
    while(cinfo.next_scanline < cinfo.image_height) {
      row_pointer[0] = (unsigned char *) (&image1[cinfo.next_scanline * nx]);
      jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }
    free(image1);
  } else {
    while(cinfo.next_scanline < cinfo.image_height) {
      row_pointer[0] = (unsigned char *) (&image[cinfo.next_scanline * nx]);
      jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }
  }

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  if(newline != 0) {
#if (defined __MINGW_H) || (defined __WIN32)
    if (file_or_socket == 0) {
      write(sock, "\n", 1);
    } else {
      send(sock, "\n", 1,0);
    }
#else
    write(sock, "\n", 1);
#endif  

  }
  return 0;
}


void jpeg_set_quality_(
		       int *newvalue
		       )
{
  jpeg_quality = *newvalue;
}

int jpeg_set_filename_(
		       char *filename,
		       int *length
		       )
{
  if(*length < 80) {
    memset(imagename, 0, 80);
    strncpy(imagename, filename, *length);
  }
  return 0;
}

/* now includes write_jpeg_spots and write_jpeg_predict */

int write_colour_jpeg_(
		       int *_fd,
		       int *_colours,
		       int *_nrec,
		       int *_iylen,
		       short *_image,
		       int *_pixmin,
		       int *_pixmax,
		       int *_imgmax,
		       int *_imgmin,
		       int *_zoom,
		       int *_smooth,
		       int *_tlc,
		       int *_brc,
		       int *_newline,
		       int *_theta,
		       int *_invert,
		       int *_reverse,
		       int *_binary,
		       int *_type,
		       short int *_spots,
		       int *_nspot,
		       short int *_predictions,
		       short int *_types,
		       int *_predict_size_x,
		       int *_predict_size_y,
		       int pixel_histogram[514]

		       )
{
  int zoom, smooth, pixels_used,xmin, xmax, ymin, ymax, dx, dy, _i, _j,
    start_pred_x, end_pred_x, start_pred_y, end_pred_y,invert_x, invert_y,
     min, max, nx, ny, nxy, i, j, k, l, temp, temp5,extra,ix,jx,knx,lny,
    this_pixel,ii,jj,jjj,jjjj,kkk,kkkk,divisor,here,factor,factor_sq,
    invert, theta, pix, file, binary, nspot,firstpixel,x_step,y_step, 
    wide, tall,imgmax, jump, talljump, widejump,jumpx,jumpy,
    temp_count,i2,i3,i4,pixelr,pixelg,pixelb,istart,ifinish;
  double total,p_x, p_y,range;
  double mostimea,mostimeb,mostime_start,mostime_end;
  triple *image, *lookup255;  
  triple pixel, _overlap, _too_wide,_full, _partial;
  short temppixel;
  int *tempimage;
  unsigned char *ucptr;
  char tempnam[20];
  //    ctimer_(&mostime_start);
  smooth = *_smooth/2;
  pixels_used = *_smooth * *_smooth;
  file = 0;
  if(*_fd == -1) {
    /* create an output file */
    if(strlen(imagename) > 0) {
#if (defined __MINGW_H) || (defined __WIN32)
      *_fd = open(imagename, O_TRUNC | O_WRONLY | O_CREAT | O_BINARY, S_IRWXU);
#else
      *_fd = open(imagename, O_TRUNC | O_WRONLY | O_CREAT , S_IRWXU);
#endif
      file = 1;
    }
  }


  /* this may be called if the file we tried to open belonged to someone 
     else, so open would return -1 */

  if(*_fd == -1) {
    sprintf(tempnam, "tempfile_XXXXXX");
#if (defined __WIN32)
    *_fd = (int) mktemp(tempnam);
#else
    *_fd = mkstemp(tempnam);
#endif
    file = 1;
  }
  /* first, see if the limits are the right way around, and if they're not,
     swap them */

  xmin = _tlc[0];
  ymin = _tlc[1];
  xmax = _brc[0];
  ymax = _brc[1];

  if(xmax < xmin) {
    temp = xmin;
    xmin = xmax;
    xmax = temp;
    invert_x = 1;
  } else {
    invert_x = 0;
  }

  if(ymax < ymin) {
    temp = ymin;
    ymin = ymax;
    ymax = temp;
    invert_y = 1;
  } else {
    invert_y = 0;
  }

  binary = *_binary;
  nspot = *_nspot;
  imgmax = *_imgmax > 0 ? *_imgmax : 32767;
  
  lookup255 = malloc(sizeof (triple) * (imgmax + 1));
  memset(lookup255,0,sizeof (triple) * (imgmax + 1));
  invert = *_invert;
  theta = *_theta;
  while (theta < 0) theta += 360;
  while (theta > 270) theta -= 360;
  
  zoom = *_zoom;
  /* 
     hrp 18.11.2008 only allow zoom = 1 or zoom out (< 0)
  */
  zoom = zoom >= 0 ? 1 : zoom;
  nx = (xmax - xmin);
  ny = (ymax - ymin);
  nxy = (*_iylen * *_nrec);
  extra = nxy%8;
  image = (triple *) malloc (sizeof(triple) * (nxy + extra));
  memset(image,0,sizeof(triple) * (nxy + extra));
  tempimage = (int *) malloc (sizeof(int) * nxy * 2);
  memset(tempimage,0,sizeof(int) * nxy * 2);
  /*  fprintf(stderr, "Writing a jpeg - zoom %d\n", zoom);
      fprintf(stderr, "Width %d Height %d\n", nx, ny);
  */
  //  ctimer_(&mostimeb);
  if((*_pixmin == *_pixmax)) {
    //    printf("manually determining min and max\n");
    /* we need to manually determine the limits on the pixel values -
       this should be done over the image as a whole if not done in Mosflm already*/
    total = 0;
    min = 1000000;
    /* Initialize array used to store pixel intensity histogram */
    for (i= 0 ; i < 511 ; i++) {
      pixel_histogram[i]= 0 ;
    }
    for(i = 0; i < nxy; i++) {
      temppixel = _image[i];
      temppixel = 0 > temppixel ? intpxl_(&temppixel) : temppixel;
      total += 0 > temppixel ? 0 : temppixel;
      if(temppixel != 0){
	min = min < temppixel ? min : temppixel;
      }
      temp5 = temp>>5;
      if(temp5 > 511)temp5 = 511;
      pixel_histogram[temp5]++ ;
    }
    max = min + 4 * ((total/nxy) - min); 
    min = 0 > min ? 0 : min;
  } else {
    max = *_pixmin > *_pixmax ? *_pixmin : *_pixmax;
    max = max >= 1 ? max : 1 ;
    min = *_pixmin < *_pixmax ? *_pixmin : *_pixmax;
    min = min >= 1 ? min : 1 ;
    *_pixmin = min;
    *_pixmax = max;
  }
  //  printf("max = %d, min = %d, imgmax = %d\n",max,min,imgmax);
  /* Store the contrast limits in the histogram array for return to GUI */
  pixel_histogram[512]= min ;
  pixel_histogram[513]= max ;

  /* Lookup table of x/255 greylevel for intensity values */ 
  range = 255.0 / (max - min);
  //    printf("max = %d, min = %d, range = %f\n",max,min,range);
  if (*_colours == 5) {
    i2 = min + (max-min)/2;
    i3 = max - (max-min)/3;
    i4 = max - (max-min)/4;
    for (i= 0 ; i <= imgmax ; i++) {
      if (i <= min) {
	lookup255[i].r = lookup255[i].g = lookup255[i].b = (greylevel) 0;
      } else if (i < i2) {
	lookup255[i].r = (greylevel) 0;
	lookup255[i].g = (greylevel) 0;
	lookup255[i].b = (greylevel) (i * 255/(i2-min));
      } else if (i < i3) {
	lookup255[i].r = (greylevel) 0;
	lookup255[i].g = (greylevel) ((i -i2) * 255/(i3-i2));
	lookup255[i].b = (greylevel) (255 - lookup255[i].g);
      } else if (i < i4) {
	lookup255[i].r = (greylevel) ((i-i3) * 255/(i4-i3));
	lookup255[i].g = (greylevel) (255 - lookup255[i].r);
	lookup255[i].b = (greylevel) 0;
      } else if (i < max) {
	lookup255[i].r = (greylevel) 255;
	lookup255[i].g = (greylevel) ((i-i4) * 255/(max-i4));
	lookup255[i].b = (greylevel) 0;
      } else if (i <= imgmax) {
	lookup255[i].r = (greylevel) 255;
	lookup255[i].g = (greylevel) 255;
	lookup255[i].b = (greylevel) ((i-max) * 255/(imgmax-max));
      } 
	  
	
    }
  } else {
    for (i= 0 ; i <= imgmax ; i++) {
      if (i <= min) {
	lookup255[i].r = lookup255[i].g = lookup255[i].b = (greylevel) 255 ;
      } else if (i >= max) {
	lookup255[i].r = lookup255[i].g = lookup255[i].b = (greylevel) 0 ;
      } else {
	if (range <= 1.0) {
	  lookup255[i].r = lookup255[i].g = lookup255[i].b = (greylevel) ((max-i) * range);
	} else {
	  lookup255[i].r = lookup255[i].g = lookup255[i].b = (greylevel) (255 - (min+i)*sqrt(range));
	}
      }
  } 

    /*    INVERT VIDEO */
    if (*_reverse == 1) {
      for (i= 0 ; i <= imgmax ; i++) {
	lookup255[i].r = (greylevel) 255 - lookup255[i].r;
	lookup255[i].g = (greylevel) 255 - lookup255[i].g;
	lookup255[i].b = (greylevel) 255 - lookup255[i].b;
      }
    }
  }
    /* next make the actual image from the shortmap */

  dx = *_iylen;
  dy = *_nrec;

  // only used for testing
  //    invert = 1;
  //   theta = 0;

  if (theta == 90 || theta == 270 ){
    wide = dx;
    tall = dy;

  } else if (theta == 0 || theta == 180){
    wide = dy;
    tall = dx;
  }
  //  printf("invert = %d, theta = %d\n",invert,theta);

  /* 
     hrp 18.11.2008 only allow zoom = 1 or zoom out 
  */
  //  ctimer_(&mostimeb);
  // cache line jump - how do we find out about this - 64 is best for Intel, 127 for PowerPC?
#ifdef __POWERPC__
  jump = 127;
#else
  jump = 64;
#endif
  talljump = dy;
  jumpx = jumpy = jump;
  /*    printf("invert = %d, theta = %d, wide = %d\n",invert,theta,wide); */
  /*
    these next constructs are almost identical - the only differences are 
    that the mapping to image[...] is different, and the assignment of jjj. 
    By putting the tests for 
    the rotation of the image outside the loops, we gain ~20% in speed
  */
  /* 0 degrees */
  //	printf("Orientation is %d %d\n",invert,theta);

  /* read (short) diffraction image into (int) array as interpreted values
   */
  for(j = 0; j < dy; j+=jump) {
    talljump -=jump;
    widejump = dx;
    jumpx = 0 > talljump ? jump + talljump : jump;
    for(i = 0; i < dx; i+=jump) {
      widejump -=jump;
      jumpy = 0 > widejump ? jump + widejump : jump;
      for (jj = 0; jj < jumpx; jj++) {
	jjj = jj + j;
	jjjj = jjj + 1;
	for (ii = 0; ii < jumpy; ii++) {
	  temppixel = _image[(jjj)*dx+i+ii];
	  /* this code will reverse the colour of the gaps between modules on PADS
	     if (temppixel < 0 ) {
	     if (*_reverse == 1) {
	     temppixel = temppixel == -1 ? 255 : intpxl_(&temppixel) + 2;
	     } else {
	     temppixel = temppixel == -1 ? 0 : intpxl_(&temppixel) + 2;
	     }
	     } 
	     tempimage[(jjj)*dx+i+ii] = temppixel; 
	  */
	  
	  tempimage[(jjj)*dx+i+ii] = temppixel > 0 ? temppixel : intpxl_(&temppixel) + 2;
	}
      }
    }
  }
  /* if we are averaging over a NxN square do it now 
   */
  for(j = smooth; j < dy-smooth; j++) {
    for(i = smooth; i < dx-smooth; i++) {
      tempimage[nxy+(j*dx)+i] = 0;
      for (kkk = -smooth; kkk <= smooth; kkk++ ) {
	for (kkkk = -smooth; kkkk <= smooth; kkkk++ ) {
	  tempimage[nxy+(j*dx)+i] += tempimage[((j+kkk)*dx)+i+kkkk];
	  //	  tempimage[nxy+(j*dx)+i] = tempimage[nxy+(j*dx)+i] > tempimage[((j+kkk)*dx)+i+kkkk] ? tempimage[nxy+(j*dx)+i] : tempimage[((j+kkk)*dx)+i+kkkk] ;
	}
      }
      tempimage[nxy+(j*dx)+i] /= pixels_used;
    }
  }
  if(invert == 0){
    if (theta == 0) {
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[(i*dy)+(wide - j - 1)] = pixel;
	}
      }
    } else if (theta == 90) {
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[nxy-((j*dx)+i) - 1] = pixel;
	}
      }
    } else if (theta == 180) {
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[nxy-((i*dy)+(wide - j -1))] = pixel;
	}
      }
    } else if(theta == 270){
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[(j*dx)+i] = pixel;
	}
      }
    }
  } else {
    /* 
       invertx, i.e. inverted left to right
    */
    if (theta == 0) {
      /* 0 degrees inverted */
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[(i*dy)+j] = pixel;
	}
      }
    } else if (theta == 90) {
      /* 90 degrees inverted */
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[nxy-((j*dx)+(wide -i)) - 1] = pixel;
	}
      }
    } else if (theta == 180) {
      /* 180 degrees inverted */
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[nxy-((i*dy)+j) - 1] = pixel;
	}
      }
    } else if (theta == 270) {
      /* 270 degrees inverted */
      for(j = 0; j < dy; j++) {
	for(i = 0; i < dx; i++) {
	  temp = tempimage[nxy + (j*dx)+i];
	  pixel = temp <= imgmax ? lookup255[temp] : lookup255[imgmax];
	  image[(j*dx)+(wide -i) - 1] = pixel;
	}
      }
    }
  }
  /*
    image[(j+jj)*dx+i+ii] = pixel;
    image[(i+ii)*dy+(wide - j -jj)] = pixel;
    image[nxy-((j+jj)*dx+i+ii)] = pixel;
    image[nxy-((i+ii)*dy+(wide - j - jj))] = pixel;
    ------ invertx ------
    image[(j+jj)*dx+(wide -i -ii)] = pixel;
    image[(i+ii)*dy+j+jj] = pixel;
    image[nxy-((j+jj)*dx+(wide -i -ii))] = pixel;
    image[nxy-((i+ii)*dy+j+jj)] = pixel;
  */
   
	
  // ctimer_(&mostimea);
  //    printf("copying to image took %6.3f seconds\n",mostimea-mostimeb);

  /* now draw the predictions onto the image */
  if(nspot > 0){
    /* set up prediction boxes */

    start_pred_x = -(*_predict_size_x-5)/4 ;
    start_pred_y = -(*_predict_size_y-5)/4;
    end_pred_x = (*_predict_size_x+5)/4;
    end_pred_y = (*_predict_size_y+5)/4;

    pixel.r = 0;
    pixel.g = 0;
    pixel.b = 255;
    
    _overlap.r = 255;
    _overlap.g = 0;
    _overlap.b = 0;
	  
    _too_wide.r = 0;
    _too_wide.g = 255;
    _too_wide.b = 0;
	  
    _full.r = 0;
    _full.g = 0;
    _full.b = 255;
	  
    _partial.r = 255;
    _partial.g = 127;
    _partial.b = 0;
	  
    /* reassign nxy to the values for the output image */
	  
    nxy = nx * ny;
    for(i = 0; i < nspot; i++) {
      l = *_nrec - _predictions[2 * i];
      k = *_iylen -_predictions[2 * i + 1];
      /* 0.5 here to draw the spot in the middle of the prediction */
      p_x = ((float) (nx * (k - xmin - 0.5)) / ((float) xmax - xmin));
      p_y = ((float) (ny * (l - ymin - 0.5)) / ((float) ymax - ymin));
      if(_types[i] == 0) {
	pixel = _full;
      } else if(_types[i] == 2) {
	pixel = _overlap;
      } else if(_types[i] == 3) {
	pixel = _too_wide;
      } else {
	pixel = _partial;
      }
      for(_i = start_pred_x; _i <= end_pred_x; _i++) {
	k = (ny - p_x + _i);
	knx = k * nx;
	for(_j = start_pred_y; _j <= end_pred_y; _j++) {
	  if (_i <= start_pred_x || _i >= end_pred_x || _j <= start_pred_y || _j >= end_pred_y ){
	    l = p_y + _j;
	    lny = l * ny;
	    if((l >= 0) && (l < ny) && (k >= 0) && (knx < nxy)) {
	      if(invert == 0) {
		if(theta == 0){
		  image[knx + l] = pixel;
		} else if(theta == 90){
		  image[lny + ny - k] = pixel;
		} else if(theta == 180){
		  image[nxy -(knx + l)] = pixel;
		} else if(theta == 270){
		  image[nxy - (lny + ny - k)] = pixel;
		} 
	      } else {
		if(theta == 0){
		  image[knx + nx - l] = pixel;
		} else if(theta == 90){
		  image[lny + k] = pixel;
		} else if(theta == 180){
		  image[nxy + l] = pixel;
		} else if(theta == 270){
		  image[nxy - (lny + k)] = pixel;
		}
	      }
	    }
	  }
	}
      }
    }	
  }
/*
	  take full-sized image and shrink it for thumbnail
	*/
  if (zoom < 0) {
    factor = zoom * -1;
    factor_sq = factor * factor;
    nx = (tall /factor);
    ny = (wide /factor);
    for(i = 0; i < nx; i++) {
      for(j = 0; j < ny; j++) {
	pixelr = pixelg = pixelb = 0;
	ii = i * factor;
	jj = j * factor;
	for(k = 0; k < factor; k++) {
	  for(l = 0; l < factor; l++) {
	    here = (ii*wide+jj) + (k*wide+l);
	    pixelr += image[here].r;
	    pixelg += image[here].g;
	    pixelb += image[here].b;
	  }
	}
	pixel.r = pixelr/factor_sq;
	pixel.g = pixelg/factor_sq;
	pixel.b = pixelb/factor_sq;
	image[(i*ny+j)] = pixel;
      }
    }
    wide = wide / factor;
    tall = tall / factor;
    //    ctimer_(&mostimeb);
    // printf("zooming out took %6.3f seconds\n",mostimeb-mostimea);
  }
	
  /* finally, all we need to do is send this image out */
  ctimer_(&mostimea);
  // printf("writing the new image array took: %6.3f\n",mostimea-mostime_start);
  // printf("*_fd is %d\n",*_fd);
  if (*_type == 2){
    mosflm_write_ppm(*_fd, wide,tall, image, *_newline, *_colours, binary);
    //  ctimer_(&mostimeb);
    //  printf("write_ppm took: %6.3f\n",mostimeb-mostimea);
  } else {
    mosflm_write_jpeg(*_fd, wide,tall,image, *_newline, *_colours, binary);
    //  ctimer_(&mostimeb);
    //  printf("write_jpeg took: %6.3f\n",mostimeb-mostimea);
  }
  ctimer_(&mostimeb);
  //  printf("write_jpeg took: %6.3f\n",mostimeb-mostimea);
  fflush(stdout);
	
  free(image);
  free(tempimage);
  free(lookup255);
  fflush(stdout);
	
  /* if we wrote to a file, close it */
  if(file == 1) {
    close(*_fd);
    *_fd = -1;
    file = 0;
  }
  //  ctimer_(&mostime_end);
  //  printf("total time spent in mosflm_colour_jpeg was %7.3f seconds\n",mostime_end-mostime_start);
  return 0;
}
