#ifndef READ_PEAKS_H

#include "util_def.h"
#include "vector_math.h"

#undef READ_INFO_DEBUG
#undef READ_RS_DEBUG
#undef READXYZ_DEBUG

Boolean readinput(char *, detector_setting *);
int read_rs(char *, struct rs_coord **);
int readxyz(char *, struct vector *);

#define READ_PEAKS_H
#endif
