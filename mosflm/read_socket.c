/* $Id: read_socket.c,v 1.14.14.1 2015/02/16 17:45:48 harry Exp $ */
/* 1 Make the necessary includes and set up the variables: */

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if (defined _MSC_VER) || (defined __WIN32)
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

int read_a_line(int sock,
		char line[])
{
  int i = 0, r;
  char l;

  while(i < 400) {
#if (defined _MSC_VER) || (defined __WIN32)
    r = recv(sock,&l,1,(int) NULL);
#else
    r = read(sock, &l, 1);
#endif
    if (r == -1) {
      /* error reading from the socket */
      return -2;
    }

    if((r == 0) && (i == 0)) {
      /* we have an eof condition */
      return -1;
    }
    if(l < ' ') {
      int j;
      for(j = i; j < 400; j++) {
      	line[j] = ' ';
      }
      return i;
    } else {
      line[i] = l;
    }
    i++;
  }
  return i;
}

int read_socket_(int *client_sockfd, 
		 char *tmpstring, 
		 int *status)
{
  int i,j;
  char *linein;
  int retval;

  linein = malloc(sizeof (char) * 400); 
  retval = read_a_line(*client_sockfd, linein);

  *status = retval <= 0 ? retval : 0;

  /* error reading from socket, exit now */
  if(*status == -2) {
    return *status;
  }
  if(retval == -1) {
    *status = -1;
  } else {
    *status = 0;
  }

  if(linein[0]=='&')
    /* close socket if string starts with '&' */
    {
      printf("C Server will close connection\n");
      write(*client_sockfd, "C Server closing connection\n", 400);
      close_socket_(*client_sockfd);
      return 1;
    }

  memset(tmpstring, ' ', 400);
  strncpy(tmpstring, linein,400);
  free(linein);
  
  return 0;
}
