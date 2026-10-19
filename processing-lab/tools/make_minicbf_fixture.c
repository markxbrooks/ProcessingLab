/* Generate a small synthetic miniCBF fixture for processing-lab tests.
 *
 * Build (from repo root):
 *   cc -O2 -std=gnu17 -Icbf/include -o processing-lab/tools/make_minicbf_fixture \
 *      processing-lab/tools/make_minicbf_fixture.c cbf/lib/libcbf.a -lm
 * Run:
 *   ./processing-lab/tools/make_minicbf_fixture processing-lab/testdata/cbf/synthetic_minicbf.cbf
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cbf.h"
#include "cbf_simple.h"

#define CHECK(expr)                                                          \
    do {                                                                     \
        int code_ = (expr);                                                  \
        if (code_ != 0) {                                                    \
            fprintf(stderr, "%s failed (%d) at line %d\n", #expr, code_,     \
                    __LINE__);                                               \
            return 1;                                                        \
        }                                                                    \
    } while (0)

/* Add one row with (column, value) pairs; repositions before each set. */
static int add_row_pairs(cbf_handle cbf, const char* const cols[],
                         const char* const vals[], int n) {
    int i;
    CHECK(cbf_new_row(cbf));
    for (i = 0; i < n; ++i) {
        CHECK(cbf_find_column(cbf, cols[i]));
        CHECK(cbf_set_value(cbf, vals[i]));
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <out.cbf>\n", argv[0]);
        return 2;
    }

    const int ndimslow = 128; /* rows (slow) */
    const int ndimfast = 128; /* cols (fast) */
    const size_t nelem = (size_t)ndimslow * (size_t)ndimfast;
    int* pixels = (int*)calloc(nelem, sizeof(int));
    if (!pixels) {
        return 3;
    }

    /* Flat background with a diamond-shaped peak near the centre. */
    const int crow = 64, ccol = 64;
    for (int r = 0; r < ndimslow; ++r) {
        for (int c = 0; c < ndimfast; ++c) {
            const int dr = r - crow, dc = c - ccol;
            const double rr = (double)(dr * dr + dc * dc);
            const double v = 50.0 + 500.0 * (rr < 36.0 ? (1.0 - rr / 36.0) : 0.0);
            pixels[(size_t)r * (size_t)ndimfast + (size_t)c] = (int)(v + 0.5);
        }
    }

    cbf_handle cbf = NULL;
    CHECK(cbf_make_handle(&cbf));
    CHECK(cbf_new_datablock(cbf, "MINICBF"));

    CHECK(cbf_require_category(cbf, "diffrn"));
    CHECK(cbf_new_column(cbf, "id"));
    {
        const char* cols[] = {"id"};
        const char* vals[] = {"SYNTHETIC"};
        if (add_row_pairs(cbf, cols, vals, 1) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }
    CHECK(cbf_set_diffrn_id(cbf, "SYNTHETIC"));

    CHECK(cbf_require_category(cbf, "diffrn_radiation"));
    CHECK(cbf_new_column(cbf, "diffrn_id"));
    CHECK(cbf_new_column(cbf, "wavelength_id"));
    {
        const char* cols[] = {"diffrn_id", "wavelength_id"};
        const char* vals[] = {"SYNTHETIC", "L1"};
        if (add_row_pairs(cbf, cols, vals, 2) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    CHECK(cbf_require_category(cbf, "diffrn_radiation_wavelength"));
    CHECK(cbf_new_column(cbf, "id"));
    CHECK(cbf_new_column(cbf, "wavelength"));
    CHECK(cbf_new_column(cbf, "wt"));
    {
        const char* cols[] = {"id", "wavelength", "wt"};
        const char* vals[] = {"L1", "0.9790", "1.0"};
        if (add_row_pairs(cbf, cols, vals, 3) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    /* Detector / element / frame skeleton for cbf_set_image. */
    CHECK(cbf_require_category(cbf, "diffrn_detector"));
    CHECK(cbf_new_column(cbf, "diffrn_id"));
    CHECK(cbf_new_column(cbf, "id"));
    {
        const char* cols[] = {"diffrn_id", "id"};
        const char* vals[] = {"SYNTHETIC", "DET0"};
        if (add_row_pairs(cbf, cols, vals, 2) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    CHECK(cbf_require_category(cbf, "diffrn_detector_element"));
    CHECK(cbf_new_column(cbf, "detector_id"));
    CHECK(cbf_new_column(cbf, "id"));
    {
        const char* cols[] = {"detector_id", "id"};
        const char* vals[] = {"DET0", "ELEM0"};
        if (add_row_pairs(cbf, cols, vals, 2) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    CHECK(cbf_require_category(cbf, "diffrn_data_frame"));
    CHECK(cbf_new_column(cbf, "detector_element_id"));
    CHECK(cbf_new_column(cbf, "array_id"));
    CHECK(cbf_new_column(cbf, "id"));
    {
        const char* cols[] = {"detector_element_id", "array_id", "id"};
        const char* vals[] = {"ELEM0", "ARR0", "FRAME0"};
        if (add_row_pairs(cbf, cols, vals, 3) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    /* array_structure_list skeleton required by cbf_set_image. */
    CHECK(cbf_require_category(cbf, "array_structure_list"));
    CHECK(cbf_new_column(cbf, "array_id"));
    CHECK(cbf_new_column(cbf, "precedence"));
    CHECK(cbf_new_column(cbf, "dimension"));
    CHECK(cbf_new_column(cbf, "type"));
    {
        const char* cols[] = {"array_id", "precedence", "dimension", "type"};
        const char* v_fast[] = {"ARR0", "1", "128", "int"};
        const char* v_slow[] = {"ARR0", "2", "128", "int"};
        if (add_row_pairs(cbf, cols, v_fast, 4) != 0 ||
            add_row_pairs(cbf, cols, v_slow, 4) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    /* array_data skeleton required by cbf_set_image (incl. `data` column). */
    CHECK(cbf_require_category(cbf, "array_data"));
    CHECK(cbf_new_column(cbf, "array_id"));
    CHECK(cbf_new_column(cbf, "binary_id"));
    CHECK(cbf_new_column(cbf, "data"));
    {
        const char* cols[] = {"array_id", "binary_id", "data"};
        const char* vals[] = {"ARR0", "1", "."};
        if (add_row_pairs(cbf, cols, vals, 3) != 0) {
            cbf_free_handle(cbf);
            free(pixels);
            return 1;
        }
    }

    /* MiniCBF-style `#` comment header (parsed by CbfImageReader fallback).
     * Real Dectris/Eiger CBFs carry these lines as the `header_contents` text
     * item inside the `array_data` category, i.e. tag
     * `_array_data.header_contents`. */
    CHECK(cbf_require_column(cbf, "header_contents"));
    CHECK(cbf_find_column(cbf, "header_contents"));
    CHECK(cbf_set_value(
        cbf,
        "# Wavelength 0.9790 A\n"
        "# Detector_distance 0.1600 m\n"
        "# Pixel_size 0.172e-3 m x 0.172e-3 m\n"
        "# Beam_xy (64.00, 64.00) pixels\n"
        "# Start_angle 0.0000 deg.\n"
        "# Angle_increment 1.0000 deg.\n"
        "# Detector: SYNTHETIC-MINI\n"
    ));

    CHECK(cbf_rewind_row(cbf));
    if (cbf_set_image(cbf, 0, 0, CBF_BYTE_OFFSET, pixels, sizeof(int), 1, ndimslow,
                      ndimfast) != 0) {
        fprintf(stderr, "cbf_set_image failed\n");
        cbf_free_handle(cbf);
        free(pixels);
        return 1;
    }

    FILE* out = fopen(argv[1], "wb");
    if (!out) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        cbf_free_handle(cbf);
        free(pixels);
        return 6;
    }
    const int write_err =
        cbf_write_file(cbf, out, 0, CBF, MIME_HEADERS, CBF_BYTE_OFFSET);
    fclose(out);
    /* CBFlib 0.9.5 double-frees binary fileconnections in cbf_free_handle
     * after cbf_write_file; the process exits immediately, so skip the free. */
    (void)cbf;
    free(pixels);
    if (write_err != 0) {
        fprintf(stderr, "cbf_write_file failed (%d)\n", write_err);
        return 7;
    }

    /* Roundtrip sanity. */
    FILE* in = fopen(argv[1], "rb");
    cbf_handle cb2 = NULL;
    cbf_make_handle(&cb2);
    if (cbf_read_file(cb2, in, MSG_DIGEST) != 0) {
        fprintf(stderr, "roundtrip read failed\n");
        fclose(in);
        cbf_free_handle(cb2);
        return 8;
    }
    double wl = 0.0;
    size_t n1 = 0, n2 = 0;
    cbf_get_wavelength(cb2, &wl);
    cbf_get_image_size(cb2, 0, 0, &n1, &n2);
    fclose(in);
    /* Same CBFlib 0.9.5 fileconnection quirk on read handles with binaries. */
    (void)cb2;
    printf("wrote %s (%dx%d, roundtrip wl=%g dims=%zux%zu)\n", argv[1],
           ndimfast, ndimslow, wl, n2, n1);
    return 0;
}
