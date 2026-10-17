#define	acosf	acos
#define	cosf	cos
#define	sinf	sin
#define	logf	log
/* #define	truncf	trunc */
#define	powf	pow
#define	sqrtf	sqrt
#define	fabsf   fabs
#if defined (__WIN32)
#define _rffti_ RFFTI
#define _rfftf_ RFFTF
#else
#define rffti_ RFFTI
#define rfftf_ RFFTF
#endif
#define invert_ INVERT
#define GEN_DIRECT_LATTICE_VECTORS_DEBUG

