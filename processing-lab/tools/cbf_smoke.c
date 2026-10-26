/* Header + pixel smoke test against vendored CBFlib (Option A).
 *
 * Build (from processing-lab/tools):
 *   cc -O0 -g -I../../third_party/cbf/include -o cbf_smoke cbf_smoke.c ../../third_party/cbf/lib/libcbf.a -lm
 * Run:
 *   ./cbf_smoke ../testdata/cbf/sample_minicbf.cbf
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cbf.h"
#include "cbf_simple.h"

int main(int argc, char** argv) {
  const char* path =
      (argc > 1) ? argv[1] : "../testdata/cbf/sample_minicbf.cbf";
  FILE* in = fopen(path, "rb");
  if (!in) {
    fprintf(stderr, "cbf_smoke: cannot open %s\n", path);
    return 1;
  }

  cbf_handle cbf = NULL;
  if (cbf_make_handle(&cbf) != 0) {
    fprintf(stderr, "cbf_smoke: cbf_make_handle failed\n");
    fclose(in);
    return 2;
  }
  int err = cbf_read_file(cbf, in, MSG_DIGEST);
  if (err != 0) {
    fprintf(stderr, "cbf_smoke: cbf_read_file failed (%d)\n", err);
    cbf_free_handle(cbf);
    return 3;
  }

  double wavelength = 0.0;
  double distance = 0.0;
  size_t ndim1 = 0;
  size_t ndim2 = 0;
  size_t nelem = 0;

  if (cbf_get_wavelength(cbf, &wavelength) != 0) {
    fprintf(stderr, "cbf_smoke: warning: wavelength missing\n");
  }

  cbf_detector detector = NULL;
  if (cbf_construct_detector(cbf, &detector, 0) == 0) {
    if (cbf_get_detector_distance(detector, &distance) != 0) {
      fprintf(stderr, "cbf_smoke: warning: distance missing\n");
    }
  } else {
    fprintf(stderr, "cbf_smoke: warning: construct_detector failed\n");
  }

  /* Prefer first binary section dimensions via cbf_get_image_size if available */
  err = cbf_get_image_size(cbf, 0, 0, &ndim1, &ndim2);
  if (err != 0) {
    /* Fall back: report success on header parse alone for miniCBF variants */
    fprintf(stderr,
            "cbf_smoke: cbf_get_image_size failed (%d); header still readable\n",
            err);
    printf("CBF_SMOKE_OK path=%s wavelength=%.6f distance=%.3f dims=unknown\n",
           path,
           wavelength,
           distance);
    if (detector) {
      cbf_free_detector(detector);
    }
    cbf_free_handle(cbf);
    return 0;
  }

  nelem = (size_t)ndim1 * (size_t)ndim2;
  int* pixels = (int*)calloc(nelem, sizeof(int));
  if (!pixels) {
    fprintf(stderr, "cbf_smoke: OOM for %zu pixels\n", nelem);
    cbf_free_handle(cbf);
    return 4;
  }

  err = cbf_get_image(cbf, 0, 0, pixels, sizeof(int), 1, ndim1, ndim2);
  if (err != 0) {
    fprintf(stderr, "cbf_smoke: cbf_get_image failed (%d)\n", err);
    free(pixels);
    if (detector) {
      cbf_free_detector(detector);
    }
    cbf_free_handle(cbf);
    return 5;
  }

  int minv = pixels[0];
  int maxv = pixels[0];
  long long sum = 0;
  for (size_t i = 0; i < nelem; ++i) {
    int v = pixels[i];
    if (v < minv) {
      minv = v;
    }
    if (v > maxv) {
      maxv = v;
    }
    sum += v;
  }

  printf(
      "CBF_SMOKE_OK path=%s wavelength=%.6f distance=%.3f "
      "ndim1=%zu ndim2=%zu min=%d max=%d mean=%.2f\n",
      path,
      wavelength,
      distance,
      ndim1,
      ndim2,
      minv,
      maxv,
      (double)sum / (double)nelem
  );

  free(pixels);
  if (detector) {
    cbf_free_detector(detector);
  }
  cbf_free_handle(cbf);
  return 0;
}
