/*======================================================================*
 * Canonical DPS / Cambridge Bravais C API (5-argument).
 *
 * Include this header (or a local bravais.h that pulls it in) so that
 * util/, index/, and src/dps/bravais/ cannot drift to incompatible
 * prototypes. See critique.md / mosflm-cpp port notes.
 *======================================================================*/

#ifndef DPS_BRAVAIS_API_H
#define DPS_BRAVAIS_API_H

#ifdef __cplusplus
extern "C" {
#endif

struct matrix_3x3;

/**
 * Classify reciprocal orientation into lattice characters 1–44.
 *
 * \param om           Input reciprocal orientation matrix (modified: reduced)
 * \param cl_om        Output orientation per character (index 1..44)
 * \param crit         Output goof / FOM per character
 * \param sys_id       Output two-char system id per character (e.g. "cP")
 * \param output_type  0 = print-oriented; 1 = fill cl_om/crit/sys_id
 */
void bravais(
    struct matrix_3x3* om,
    struct matrix_3x3 cl_om[45],
    float crit[45],
    char sys_id[45][3],
    int output_type
);

#ifdef __cplusplus
}
#endif

#endif /* DPS_BRAVAIS_API_H */
