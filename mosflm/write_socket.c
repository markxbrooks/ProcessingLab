/* $Id: write_socket.c,v 1.12.10.3 2014/12/02 15:06:48 harry Exp $
   1 Make the necessary includes and set up the variables: */

#include <sys/types.h>
#include <stdio.h>
#if (defined _MSC_VER) || (defined __WIN32)
#include <winsock.h>
#include "filenotsocket.h"
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

int write_socket_length_(
  int *fd,
  int *length,
  char *line
)
{
#if (defined _MSC_VER) || (defined __WIN32)
  if (filenotsocket == 1){
  send(*fd, line, *length,0);
  send(*fd, "\n", 1,0);
  } else {
  write(*fd, line, *length,0);
  write(*fd, "\n", 1,0);
  }
#else
  write(*fd, line, *length);
  write(*fd, "\n", 1);
#endif

  return 0;
}

int write_socket_section_(
  int *fd,
  int *length,
  char *line
)
{
#if (defined _MSC_VER) || (defined __WIN32)
  if (filenotsocket == 1){
  send(*fd, line, *length,0);
  } else {
  write(*fd, line, *length,0);
  }
#else
  write(*fd, line, *length);
#endif
  return 0;
}
