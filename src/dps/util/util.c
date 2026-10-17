#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

/*
 *      makestring()  --  allocates space and duplicates the given string,
 *                        printing an error message in case of failure.
 */
char *makestring( char *s )

{
	char  *t;

	if (s == NULL)
		return NULL;
	t = (char *) malloc( strlen(s) + 1 );
	if (t == NULL) {
		fprintf(stderr,"makestring: out of memory\n");
		exit(1);
	}
	else
		strcpy( t, s );
	return t;
}

/* upper(): Converts a string to upper case */

void upper(char *s)

{
	int i;
	for (i=0; i<strlen(s); i++) s[i]=toupper(s[i]);
}

