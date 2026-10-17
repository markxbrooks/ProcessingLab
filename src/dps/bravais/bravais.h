/*======================================================================*
 *                   Data Processing Suit                           	*
 * Written by Ingo Steller, Department of Biological Sciences       	*
 *                          Purdue University                       	*
 *                          West Lafayette, Indiana (USA)           	*
 *                                                                      *
 * Wrapper: helpers live next to this file; the bravais() prototype is
 * defined once in include/dps/bravais_api.h (5-arg Cambridge API).
 *======================================================================*/

#ifndef BRAVAIS_H

#include "bravais_help.h"
#include "bravais_goof.h"
#include "bravais_permutate.h"
#include "bravais_transform.h"

#define BRAVAIS_DEBUG

typedef short BravaisMode;

#define BRAVAIS_COMPARE 0
#define BRAVAIS_DETERMINE 1

#include "../../../include/dps/bravais_api.h"

#define BRAVAIS_H
#endif
