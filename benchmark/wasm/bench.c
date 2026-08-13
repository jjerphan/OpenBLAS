/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, THE PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/*
 * Node/Emscripten microbenchmark driver for WASM128_GENERIC.
 *
 * Usage: bench <op> <n>
 * Prints one CSV line: op,n,inner,median_s,mflops,min_s,max_s
 *
 * Environment:
 *   WASM_BENCH_WARMUP  (default 5)
 *   WASM_BENCH_TIMED   (default 10)
 *   WASM_BENCH_TARGET_S (default 0.05 seconds per timed sample)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cblas.h"

static volatile double g_sink;

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int env_int(const char *name, int fallback) {
  const char *p = getenv(name);
  return p ? atoi(p) : fallback;
}

static double env_double(const char *name, double fallback) {
  const char *p = getenv(name);
  return p ? atof(p) : fallback;
}

static int cmp_double(const void *a, const void *b) {
  double da = *(const double *)a;
  double db = *(const double *)b;
  return (da > db) - (da < db);
}

static void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p) {
    fprintf(stderr, "out of memory (%zu bytes)\n", n);
    exit(1);
  }
  return p;
}

static void fill_f32(float *a, size_t n, unsigned seed) {
  unsigned s = seed ? seed : 1u;
  for (size_t i = 0; i < n; i++) {
    s = s * 1664525u + 1013904223u;
    a[i] = (float)((int)(s >> 16) % 2001 - 1000) * 0.001f;
  }
}

static void fill_f64(double *a, size_t n, unsigned seed) {
  unsigned s = seed ? seed : 1u;
  for (size_t i = 0; i < n; i++) {
    s = s * 1664525u + 1013904223u;
    a[i] = (double)((int)(s >> 16) % 2001 - 1000) * 0.001;
  }
}

static void tri_f32(float *a, int n) {
  memset(a, 0, (size_t)n * (size_t)n * sizeof(float));
  for (int j = 0; j < n; j++) {
    a[j + j * n] = 1.0f;
    for (int i = j + 1; i < n; i++)
      a[i + j * n] = 0.1f;
  }
}

static void tri_f64(double *a, int n) {
  memset(a, 0, (size_t)n * (size_t)n * sizeof(double));
  for (int j = 0; j < n; j++) {
    a[j + j * n] = 1.0;
    for (int i = j + 1; i < n; i++)
      a[i + j * n] = 0.1;
  }
}

typedef void (*kernel_fn)(void *ctx, int inner);

struct ctx_s {
  int n;
  float *a, *b, *c, *x, *y;
  float alpha, beta;
};

struct ctx_d {
  int n;
  double *a, *b, *c, *x, *y;
  double alpha, beta;
};

#define DEF_SGEMM(name, transa, transb)                                        \
  static void name(void *v, int inner) {                                       \
    struct ctx_s *c = v;                                                       \
    int n = c->n;                                                              \
    for (int i = 0; i < inner; i++)                                            \
      cblas_sgemm(CblasColMajor, transa, transb, n, n, n, c->alpha, c->a, n,   \
                  c->b, n, c->beta, c->c, n);                                  \
    g_sink = c->c[0] + c->c[(size_t)n * n - 1];                                \
  }

DEF_SGEMM(k_sgemm, CblasNoTrans, CblasNoTrans)

static void k_dgemm(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, c->alpha,
                c->a, n, c->b, n, c->beta, c->c, n);
  g_sink = c->c[0] + c->c[(size_t)n * n - 1];
}

static void k_sgemv_n(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_sgemv(CblasColMajor, CblasNoTrans, n, n, c->alpha, c->a, n, c->x, 1,
                c->beta, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_sgemv_t(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_sgemv(CblasColMajor, CblasTrans, n, n, c->alpha, c->a, n, c->x, 1,
                c->beta, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_dgemv_n(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, c->alpha, c->a, n, c->x, 1,
                c->beta, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_dgemv_t(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, c->alpha, c->a, n, c->x, 1,
                c->beta, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_saxpy(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_saxpy(n, c->alpha, c->x, 1, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_daxpy(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_daxpy(n, c->alpha, c->x, 1, c->y, 1);
  g_sink = c->y[0] + c->y[n - 1];
}

static void k_sdot(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  float acc = 0.f;
  for (int i = 0; i < inner; i++)
    acc += cblas_sdot(n, c->x, 1, c->y, 1);
  g_sink = acc;
}

static void k_ddot(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  double acc = 0.0;
  for (int i = 0; i < inner; i++)
    acc += cblas_ddot(n, c->x, 1, c->y, 1);
  g_sink = acc;
}

static void k_sscal(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_sscal(n, c->alpha, c->x, 1);
  g_sink = c->x[0] + c->x[n - 1];
}

static void k_dscal(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_dscal(n, c->alpha, c->x, 1);
  g_sink = c->x[0] + c->x[n - 1];
}

static void k_strsm(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++) {
    memcpy(c->c, c->b, (size_t)n * (size_t)n * sizeof(float));
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                n, n, c->alpha, c->a, n, c->c, n);
  }
  g_sink = c->c[0] + c->c[(size_t)n * n - 1];
}

static void k_dtrsm(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++) {
    memcpy(c->c, c->b, (size_t)n * (size_t)n * sizeof(double));
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                n, n, c->alpha, c->a, n, c->c, n);
  }
  g_sink = c->c[0] + c->c[(size_t)n * n - 1];
}

static void k_ssyrk(void *v, int inner) {
  struct ctx_s *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, n, n, c->alpha, c->a, n,
                c->beta, c->c, n);
  g_sink = c->c[0] + c->c[(size_t)n * n - 1];
}

static void k_dsyrk(void *v, int inner) {
  struct ctx_d *c = v;
  int n = c->n;
  for (int i = 0; i < inner; i++)
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, n, n, c->alpha, c->a, n,
                c->beta, c->c, n);
  g_sink = c->c[0] + c->c[(size_t)n * n - 1];
}

static double flops_for(const char *op, int n) {
  double N = (double)n;
  if (!strcmp(op, "sgemm") || !strcmp(op, "dgemm"))
    return 2.0 * N * N * N;
  if (!strcmp(op, "sgemv") || !strcmp(op, "dgemv") || !strcmp(op, "sgemv_t") ||
      !strcmp(op, "dgemv_t"))
    return 2.0 * N * N;
  if (!strcmp(op, "saxpy") || !strcmp(op, "daxpy") || !strcmp(op, "sdot") ||
      !strcmp(op, "ddot"))
    return 2.0 * N;
  if (!strcmp(op, "sscal") || !strcmp(op, "dscal"))
    return N;
  if (!strcmp(op, "strsm") || !strcmp(op, "dtrsm"))
    return N * N * N;
  if (!strcmp(op, "ssyrk") || !strcmp(op, "dsyrk"))
    return N * N * N;
  return 0.0;
}

static void run_op(const char *op, int n) {
  int warmup = env_int("WASM_BENCH_WARMUP", 5);
  int timed = env_int("WASM_BENCH_TIMED", 10);
  double target = env_double("WASM_BENCH_TARGET_S", 0.05);

  kernel_fn fn = NULL;
  void *ctx = NULL;
  struct ctx_s cs;
  struct ctx_d cd;
  memset(&cs, 0, sizeof(cs));
  memset(&cd, 0, sizeof(cd));

  int is_d = (op[0] == 'd');
  int need_mat = strstr(op, "gemm") || strstr(op, "gemv") || strstr(op, "trsm") ||
                 strstr(op, "syrk");
  int need_vec = strstr(op, "axpy") || strstr(op, "dot") || strstr(op, "scal") ||
                 strstr(op, "gemv");

  if (is_d) {
    cd.n = n;
    cd.alpha = 1.0;
    cd.beta = 0.0;
    if (need_mat) {
      cd.a = xmalloc((size_t)n * n * sizeof(double));
      cd.b = xmalloc((size_t)n * n * sizeof(double));
      cd.c = xmalloc((size_t)n * n * sizeof(double));
      if (!strcmp(op, "dtrsm")) {
        tri_f64(cd.a, n);
        fill_f64(cd.b, (size_t)n * n, 2);
        memcpy(cd.c, cd.b, (size_t)n * n * sizeof(double));
      } else {
        fill_f64(cd.a, (size_t)n * n, 1);
        fill_f64(cd.b, (size_t)n * n, 2);
        fill_f64(cd.c, (size_t)n * n, 3);
      }
    }
    if (need_vec) {
      cd.x = xmalloc((size_t)n * sizeof(double));
      cd.y = xmalloc((size_t)n * sizeof(double));
      fill_f64(cd.x, n, 4);
      fill_f64(cd.y, n, 5);
    }
    ctx = &cd;
    if (!strcmp(op, "dgemm"))
      fn = k_dgemm;
    else if (!strcmp(op, "dgemv"))
      fn = k_dgemv_n;
    else if (!strcmp(op, "dgemv_t"))
      fn = k_dgemv_t;
    else if (!strcmp(op, "daxpy"))
      fn = k_daxpy;
    else if (!strcmp(op, "ddot"))
      fn = k_ddot;
    else if (!strcmp(op, "dscal")) {
      cd.alpha = 0.999;
      fn = k_dscal;
    } else if (!strcmp(op, "dtrsm"))
      fn = k_dtrsm;
    else if (!strcmp(op, "dsyrk"))
      fn = k_dsyrk;
  } else {
    cs.n = n;
    cs.alpha = 1.0f;
    cs.beta = 0.0f;
    if (need_mat) {
      cs.a = xmalloc((size_t)n * n * sizeof(float));
      cs.b = xmalloc((size_t)n * n * sizeof(float));
      cs.c = xmalloc((size_t)n * n * sizeof(float));
      if (!strcmp(op, "strsm")) {
        tri_f32(cs.a, n);
        fill_f32(cs.b, (size_t)n * n, 2);
        memcpy(cs.c, cs.b, (size_t)n * n * sizeof(float));
      } else {
        fill_f32(cs.a, (size_t)n * n, 1);
        fill_f32(cs.b, (size_t)n * n, 2);
        fill_f32(cs.c, (size_t)n * n, 3);
      }
    }
    if (need_vec) {
      cs.x = xmalloc((size_t)n * sizeof(float));
      cs.y = xmalloc((size_t)n * sizeof(float));
      fill_f32(cs.x, n, 4);
      fill_f32(cs.y, n, 5);
    }
    ctx = &cs;
    if (!strcmp(op, "sgemm"))
      fn = k_sgemm;
    else if (!strcmp(op, "sgemv"))
      fn = k_sgemv_n;
    else if (!strcmp(op, "sgemv_t"))
      fn = k_sgemv_t;
    else if (!strcmp(op, "saxpy"))
      fn = k_saxpy;
    else if (!strcmp(op, "sdot"))
      fn = k_sdot;
    else if (!strcmp(op, "sscal")) {
      cs.alpha = 0.999f;
      fn = k_sscal;
    } else if (!strcmp(op, "strsm"))
      fn = k_strsm;
    else if (!strcmp(op, "ssyrk"))
      fn = k_ssyrk;
  }

  if (!fn) {
    fprintf(stderr, "unknown op: %s\n", op);
    exit(2);
  }

  fn(ctx, 1);
  double t0 = now_s();
  fn(ctx, 1);
  double one = now_s() - t0;
  int inner = 1;
  if (one > 0.0 && one < target)
    inner = (int)(target / one);
  if (inner < 1)
    inner = 1;
  if (inner > 100000)
    inner = 100000;

  for (int i = 0; i < warmup; i++)
    fn(ctx, inner);

  double *samples = xmalloc((size_t)timed * sizeof(double));
  for (int i = 0; i < timed; i++) {
    double a = now_s();
    fn(ctx, inner);
    samples[i] = (now_s() - a) / (double)inner;
  }

  qsort(samples, (size_t)timed, sizeof(double), cmp_double);
  double median = samples[timed / 2];
  double minv = samples[0];
  double maxv = samples[timed - 1];
  double mflops = flops_for(op, n) / median * 1e-6;

  printf("%s,%d,%d,%.9f,%.4f,%.9f,%.9f\n", op, n, inner, median, mflops, minv,
         maxv);

  free(samples);
  free(cs.a);
  free(cs.b);
  free(cs.c);
  free(cs.x);
  free(cs.y);
  free(cd.a);
  free(cd.b);
  free(cd.c);
  free(cd.x);
  free(cd.y);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <op> <n>\n", argv[0]);
    return 2;
  }
  run_op(argv[1], atoi(argv[2]));
  return 0;
}
