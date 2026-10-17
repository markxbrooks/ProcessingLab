#include <stdio.h>
#include "filec.h"

static int filec_debug;

void filec_setdebug (int state)
{
	filec_debug = state;
	if ( filec_debug )
		printf ("Debugging information for filec turned on\n");
}

int filec_getdebug (void)
{
	return filec_debug;
}
