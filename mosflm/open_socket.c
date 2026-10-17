/* $Id: open_socket.c,v 1.11.10.2 2014/11/13 18:19:31 harry Exp $
   This is the Mosflm socket CLIENT; it is started up from within mosflm.f
   if the appropriate command line argument(s) have been used.
   These are: 
              REMOTEMOSFLMHOST (optional, defaults to local host)
	      MOSFLMSOCKET     (required)
*/

/* 1 Make the necessary includes and set up the variables: */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#if (defined _MSC_VER) || (defined __WIN32)
#include <winsock.h>
#define FILEORSOCKET
#include "filenotsocket.h"
#define bzero(a,b) memset(a,0,b)
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#endif
#if (defined __WIN32)
#include <winsock2.h>
#endif

int open_socket_(char *remotemosflmhost, int *portnumber, int *serverfd){
    int sockfd;
    int len;
    struct sockaddr_in address;
    int result;
    char *host, **names, **addrs;
    struct hostent *hostinfo;
    char forhostnum[100];
#if (defined _MSC_VER) || (defined __WIN32)
    WORD wVersionRequested;
    WSADATA wsaData;
    int err;
    filenotsocket = 1;

#endif
    /*********************************************************************
                GETNAME stuff 
    **********************************************************************/
    /* (a) associate an IP address with the Host name.... */
#if (defined _MSC_VER) || (defined __WIN32)
    wVersionRequested = MAKEWORD( 1, 1 );
    err = WSAStartup( wVersionRequested, &wsaData );
#endif
    /* we're going to be using the loopback address 
       host = remotemosflmhost;
       gethostname doesn't seem to work on Linux for localhost!
       gethostname(host, 255);

    */
    host = "127.0.0.1";
    remotemosflmhost = host;
	hostinfo = gethostbyname(host);
	
	if(!hostinfo) {
	  fprintf(stderr, "OPEN_SOCKET: cannot get info for host: %s\n", host);
	  exit(1);
	}
	/* (b) Display the hostname and any aliases it may have: */

    printf("OPEN_SOCKET: Connecting to host %s:\n", host);
    printf("OPEN_SOCKET: Full name: %s\n", hostinfo -> h_name);
    names = hostinfo -> h_aliases;
    printf("\n");



    /* (c) Warn and exit if the host in question is not an IP host: */

    if(hostinfo -> h_addrtype != AF_INET) {
        fprintf(stderr, "OPEN_SOCKET: not an IP host!\n");
        exit(1);
    }



    /* (d) Otherwise, display the IP address(es): */

    addrs = hostinfo -> h_addr_list;
    while(*addrs) {
	strcpy(forhostnum,inet_ntoa(*(struct in_addr *)*addrs));
        addrs++;
    }
    printf("OPEN_SOCKET: \n");

    /************************************************************************/
    /* 2 Create a socket for the client: */

    sockfd = socket(AF_INET, SOCK_STREAM, 0);

    /* 3 Name the socket, as agreed with the server: */

    bzero(&address, sizeof(struct sockaddr_in));
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = inet_addr(forhostnum);
    address.sin_port = htons(*portnumber);
    /* address.sin_port = *portnumber; */
    len = sizeof(address);

    /* 4 Now connect our socket to the server's socket: */

    result = connect(sockfd, (struct sockaddr *)&address, len);
    if(result == -1) {
      perror("OPEN_SOCKET: oops: client1");
      exit(1);
    }
    *serverfd = sockfd;
    return(0);
}

void open_file_(int *_portnumber, int *serverfd, char *xmlfile, int *length){
  char *localfile;
#if (defined _MSC_VER) || (defined __WIN32)
  filenotsocket = 0;
#endif
  localfile = malloc(sizeof(char)* (*length + 1));
  strncpy(localfile, xmlfile,*length);
  localfile[*length] = '\0';
  if(*_portnumber == -1) {
    /* create an output file */
    if(strlen(localfile) > 0) {
#if (defined __MINGW_H)
      *_portnumber = open(localfile, O_TRUNC | O_WRONLY | O_CREAT | O_BINARY, S_IRWXU);
#else
      *_portnumber = open(localfile, O_TRUNC | O_WRONLY | O_CREAT , S_IRWXU);
#endif
    }
  }
  *serverfd = *_portnumber;
  free (localfile);
}
