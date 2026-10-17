/* 1 Make the necessary includes and set up the variables: */

#include <sys/types.h>
#include <stdio.h>

#if (defined _MSC_VER) || (defined __WIN32)
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

int  close_socket_(int *client_sockfd)
{
    int err;
    const unsigned int one = 1;
    const unsigned char charone = '1';
    char closeline[26] = "Server closing connection\0";
    write(*client_sockfd, closeline, 26);
    close(*client_sockfd);
#if (defined _MSC_VER) || (defined __WIN32)
    setsockopt(*client_sockfd, SOL_SOCKET, SO_REUSEADDR, &charone, sizeof(one));
    err = WSACleanup();
#else
    setsockopt(*client_sockfd, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
#endif
    return 0;
}
