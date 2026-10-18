/* MD5.H - header file for MD5C.C
 */

/* Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
rights reserved.

License to copy and use this software is granted provided that it
is identified as the "RSA Data Security, Inc. MD5 Message-Digest
Algorithm" in all material mentioning or referencing this software
or this function.

License is also granted to make and use derivative works provided
that such works are identified as "derived from the RSA Data
Security, Inc. MD5 Message-Digest Algorithm" in all material
mentioning or referencing the derived work.

RSA Data Security, Inc. makes no representations concerning either
the merchantability of this software or the suitability of this
software for any particular purpose. It is provided "as is"
without express or implied warranty of any kind.

These notices must be retained in any copies of any part of this
documentation and/or software.
 */

#ifndef _MD5_H_
#define _MD5_H_

#include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/* MD5 context. */
typedef struct {
  UINT4 state[4];                                   /* state (ABCD) */
  UINT4 count[2];        /* number of bits, modulo 2^64 (lsb first) */
  unsigned char buffer[64];                         /* input buffer */
} MD5_CTX;

/* ANSI C prototypes (match md5c.c). Do not use PROTO_LIST here:
 * CBFlib's global.h historically defaults PROTOTYPES to 0, which made
 * PROTO_LIST expand to () and triggered -Wstrict-prototypes / C23
 * -Wdeprecated-non-prototype on call sites such as cbf_codes.c. */
void MD5Init(MD5_CTX *context);
/* ``input`` is treated as an opaque byte buffer (``const void *`` so
 * ``char *`` / ``unsigned char *`` call sites type-check cleanly). */
void MD5Update(MD5_CTX *context, const void *input, unsigned int inputLen);
void MD5Final(unsigned char digest[16], MD5_CTX *context);

#ifdef __cplusplus
}
#endif

#endif
