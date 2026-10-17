/* $Id: filenotsocket.h,v 1.1 2014/11/13 17:51:16 harry Exp $
  filenotsocket - neede for Windows so we can use a file or a socket for 
  XML output. UNIX type OS's allow a "write" to either, but Windows 
  demands a "send" to a socket and a "write" to a file */
#ifndef FILENOTSOCKET
#define FILENOTSOCKET
#ifdef FILEORSOCKET
int filenotsocket;
#else
extern int filenotsocket;
#endif
#endif
