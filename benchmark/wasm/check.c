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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* CBLAS correctness checks for odd sizes, transposes, and TRSM sides. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cblas.h"

static int g_fail;

static void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p) {
    fprintf(stderr, "oom\n");
    exit(1);
  }
  return p;
}

static void fill_f32(float *a, int n, int seed) {
  for (int i = 0; i < n; i++)
    a[i] = (float)((i * 17 + seed * 13) % 19) / 19.0f - 0.5f;
}

static void fill_f64(double *a, int n, int seed) {
  for (int i = 0; i < n; i++)
    a[i] = (double)((i * 17 + seed * 13) % 19) / 19.0 - 0.5;
}

static void lower_unit_f32(float *a, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (i == j)
        a[i + j * n] = 1.0f;
      else if (i > j)
        a[i + j * n] = 0.1f * (float)((i + j) % 5 + 1);
      else
        a[i + j * n] = 0.0f;
    }
}

static void lower_unit_f64(double *a, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (i == j)
        a[i + j * n] = 1.0;
      else if (i > j)
        a[i + j * n] = 0.1 * (double)((i + j) % 5 + 1);
      else
        a[i + j * n] = 0.0;
    }
}

static void upper_unit_f32(float *a, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (i == j)
        a[i + j * n] = 1.0f;
      else if (i < j)
        a[i + j * n] = 0.1f * (float)((i + j) % 5 + 1);
      else
        a[i + j * n] = 0.0f;
    }
}

static void upper_unit_f64(double *a, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (i == j)
        a[i + j * n] = 1.0;
      else if (i < j)
        a[i + j * n] = 0.1 * (double)((i + j) % 5 + 1);
      else
        a[i + j * n] = 0.0;
    }
}

static void gemm_ref_f32(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int n,
                         float alpha, const float *A, const float *B, float *C) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      float s = 0.0f;
      for (int k = 0; k < n; k++) {
        float av = (ta == CblasNoTrans) ? A[i + k * n] : A[k + i * n];
        float bv = (tb == CblasNoTrans) ? B[k + j * n] : B[j + k * n];
        s += av * bv;
      }
      C[i + j * n] = alpha * s;
    }
}

static void gemm_ref_f64(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int n,
                         double alpha, const double *A, const double *B,
                         double *C) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double s = 0.0;
      for (int k = 0; k < n; k++) {
        double av = (ta == CblasNoTrans) ? A[i + k * n] : A[k + i * n];
        double bv = (tb == CblasNoTrans) ? B[k + j * n] : B[j + k * n];
        s += av * bv;
      }
      C[i + j * n] = alpha * s;
    }
}

static int close_f32(const float *got, const float *ref, int n, float tol) {
  float maxe = 0.0f, maxv = 0.0f;
  for (int i = 0; i < n; i++) {
    float e = fabsf(got[i] - ref[i]);
    float v = fabsf(ref[i]);
    if (e > maxe)
      maxe = e;
    if (v > maxv)
      maxv = v;
  }
  return maxe <= tol * (1.0f + maxv);
}

static int close_f64(const double *got, const double *ref, int n, double tol) {
  double maxe = 0.0, maxv = 0.0;
  for (int i = 0; i < n; i++) {
    double e = fabs(got[i] - ref[i]);
    double v = fabs(ref[i]);
    if (e > maxe)
      maxe = e;
    if (v > maxv)
      maxv = v;
  }
  return maxe <= tol * (1.0 + maxv);
}

static void check_sgemm(int n, enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb) {
  size_t nn = (size_t)n * (size_t)n;
  float *A = xmalloc(nn * sizeof(float));
  float *B = xmalloc(nn * sizeof(float));
  float *C = xmalloc(nn * sizeof(float));
  float *R = xmalloc(nn * sizeof(float));
  fill_f32(A, (int)nn, 1);
  fill_f32(B, (int)nn, 2);
  memset(C, 0, nn * sizeof(float));
  gemm_ref_f32(ta, tb, n, 1.1f, A, B, R);
  cblas_sgemm(CblasColMajor, ta, tb, n, n, n, 1.1f, A, n, B, n, 0.0f, C, n);
  if (!close_f32(C, R, (int)nn, 2e-4f * (float)n)) {
    fprintf(stderr, "FAIL sgemm n=%d ta=%d tb=%d\n", n, (int)ta, (int)tb);
    g_fail++;
  }
  free(A);
  free(B);
  free(C);
  free(R);
}

static void check_dgemm(int n, enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb) {
  size_t nn = (size_t)n * (size_t)n;
  double *A = xmalloc(nn * sizeof(double));
  double *B = xmalloc(nn * sizeof(double));
  double *C = xmalloc(nn * sizeof(double));
  double *R = xmalloc(nn * sizeof(double));
  fill_f64(A, (int)nn, 1);
  fill_f64(B, (int)nn, 2);
  memset(C, 0, nn * sizeof(double));
  gemm_ref_f64(ta, tb, n, 1.1, A, B, R);
  cblas_dgemm(CblasColMajor, ta, tb, n, n, n, 1.1, A, n, B, n, 0.0, C, n);
  if (!close_f64(C, R, (int)nn, 2e-12 * (double)n)) {
    fprintf(stderr, "FAIL dgemm n=%d ta=%d tb=%d\n", n, (int)ta, (int)tb);
    g_fail++;
  }
  free(A);
  free(B);
  free(C);
  free(R);
}

/* Left, lower, non-unit: solve L X = B (col-major). */
static void trsm_lln_ref_f32(int n, const float *L, float *X) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      float s = X[i + j * n];
      for (int k = 0; k < i; k++)
        s -= L[i + k * n] * X[k + j * n];
      X[i + j * n] = s / L[i + i * n];
    }
}

static void trsm_lln_ref_f64(int n, const double *L, double *X) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double s = X[i + j * n];
      for (int k = 0; k < i; k++)
        s -= L[i + k * n] * X[k + j * n];
      X[i + j * n] = s / L[i + i * n];
    }
}

static void check_strsm(int n, enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE t) {
  size_t nn = (size_t)n * (size_t)n;
  float *A = xmalloc(nn * sizeof(float));
  float *B = xmalloc(nn * sizeof(float));
  float *X = xmalloc(nn * sizeof(float));
  lower_unit_f32(A, n);
  fill_f32(B, (int)nn, 3);
  memcpy(X, B, nn * sizeof(float));
  cblas_strsm(CblasColMajor, side, CblasLower, t, CblasNonUnit, n, n, 1.0f, A, n,
              X, n);
  if (side == CblasLeft && t == CblasNoTrans) {
    float *R = xmalloc(nn * sizeof(float));
    memcpy(R, B, nn * sizeof(float));
    trsm_lln_ref_f32(n, A, R);
    if (!close_f32(X, R, (int)nn, 3e-4f * (float)n)) {
      fprintf(stderr, "FAIL strsm LN n=%d\n", n);
      g_fail++;
    }
    free(R);
  } else {
    /* Reconstruct B_hat = op(A) X or X op(A) and compare to original B. */
    float *Bhat = xmalloc(nn * sizeof(float));
    memset(Bhat, 0, nn * sizeof(float));
    if (side == CblasLeft)
      gemm_ref_f32(t, CblasNoTrans, n, 1.0f, A, X, Bhat);
    else
      gemm_ref_f32(CblasNoTrans, t, n, 1.0f, X, A, Bhat);
    if (!close_f32(Bhat, B, (int)nn, 3e-4f * (float)n)) {
      fprintf(stderr, "FAIL strsm side=%d t=%d n=%d\n", (int)side, (int)t, n);
      g_fail++;
    }
    free(Bhat);
  }
  free(A);
  free(B);
  free(X);
}

static void check_dtrsm(int n, enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE t) {
  size_t nn = (size_t)n * (size_t)n;
  double *A = xmalloc(nn * sizeof(double));
  double *B = xmalloc(nn * sizeof(double));
  double *X = xmalloc(nn * sizeof(double));
  lower_unit_f64(A, n);
  fill_f64(B, (int)nn, 3);
  memcpy(X, B, nn * sizeof(double));
  cblas_dtrsm(CblasColMajor, side, CblasLower, t, CblasNonUnit, n, n, 1.0, A, n,
              X, n);
  if (side == CblasLeft && t == CblasNoTrans) {
    double *R = xmalloc(nn * sizeof(double));
    memcpy(R, B, nn * sizeof(double));
    trsm_lln_ref_f64(n, A, R);
    if (!close_f64(X, R, (int)nn, 3e-12 * (double)n)) {
      fprintf(stderr, "FAIL dtrsm LN n=%d\n", n);
      g_fail++;
    }
    free(R);
  } else {
    double *Bhat = xmalloc(nn * sizeof(double));
    memset(Bhat, 0, nn * sizeof(double));
    if (side == CblasLeft)
      gemm_ref_f64(t, CblasNoTrans, n, 1.0, A, X, Bhat);
    else
      gemm_ref_f64(CblasNoTrans, t, n, 1.0, X, A, Bhat);
    if (!close_f64(Bhat, B, (int)nn, 3e-12 * (double)n)) {
      fprintf(stderr, "FAIL dtrsm side=%d t=%d n=%d\n", (int)side, (int)t, n);
      g_fail++;
    }
    free(Bhat);
  }
  free(A);
  free(B);
  free(X);
}

static void check_strmm(int n, enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                        enum CBLAS_TRANSPOSE t) {
  size_t nn = (size_t)n * (size_t)n;
  float *A = xmalloc(nn * sizeof(float));
  float *B = xmalloc(nn * sizeof(float));
  float *X = xmalloc(nn * sizeof(float));
  float *R = xmalloc(nn * sizeof(float));
  if (uplo == CblasLower)
    lower_unit_f32(A, n);
  else
    upper_unit_f32(A, n);
  fill_f32(B, (int)nn, 4);
  memcpy(X, B, nn * sizeof(float));
  cblas_strmm(CblasColMajor, side, uplo, t, CblasNonUnit, n, n, 1.1f, A, n, X,
              n);
  if (side == CblasLeft)
    gemm_ref_f32(t, CblasNoTrans, n, 1.1f, A, B, R);
  else
    gemm_ref_f32(CblasNoTrans, t, n, 1.1f, B, A, R);
  if (!close_f32(X, R, (int)nn, 3e-4f * (float)n)) {
    fprintf(stderr, "FAIL strmm side=%d uplo=%d t=%d n=%d\n", (int)side,
            (int)uplo, (int)t, n);
    g_fail++;
  }
  free(A);
  free(B);
  free(X);
  free(R);
}

static void check_dtrmm(int n, enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                        enum CBLAS_TRANSPOSE t) {
  size_t nn = (size_t)n * (size_t)n;
  double *A = xmalloc(nn * sizeof(double));
  double *B = xmalloc(nn * sizeof(double));
  double *X = xmalloc(nn * sizeof(double));
  double *R = xmalloc(nn * sizeof(double));
  if (uplo == CblasLower)
    lower_unit_f64(A, n);
  else
    upper_unit_f64(A, n);
  fill_f64(B, (int)nn, 4);
  memcpy(X, B, nn * sizeof(double));
  cblas_dtrmm(CblasColMajor, side, uplo, t, CblasNonUnit, n, n, 1.1, A, n, X, n);
  if (side == CblasLeft)
    gemm_ref_f64(t, CblasNoTrans, n, 1.1, A, B, R);
  else
    gemm_ref_f64(CblasNoTrans, t, n, 1.1, B, A, R);
  if (!close_f64(X, R, (int)nn, 3e-12 * (double)n)) {
    fprintf(stderr, "FAIL dtrmm side=%d uplo=%d t=%d n=%d\n", (int)side,
            (int)uplo, (int)t, n);
    g_fail++;
  }
  free(A);
  free(B);
  free(X);
  free(R);
}

int main(void) {
  const int sizes[] = {1, 2, 3, 4, 5, 7, 8, 15, 16, 17, 31, 32, 33};
  const enum CBLAS_TRANSPOSE tr[4] = {CblasNoTrans, CblasTrans, CblasNoTrans,
                                      CblasTrans};
  int ns = (int)(sizeof(sizes) / sizeof(sizes[0]));
  for (int s = 0; s < ns; s++) {
    int n = sizes[s];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
        check_sgemm(n, tr[i], tr[j]);
        check_dgemm(n, tr[i], tr[j]);
      }
    check_strsm(n, CblasLeft, CblasNoTrans);
    check_strsm(n, CblasLeft, CblasTrans);
    check_strsm(n, CblasRight, CblasNoTrans);
    check_strsm(n, CblasRight, CblasTrans);
    check_dtrsm(n, CblasLeft, CblasNoTrans);
    check_dtrsm(n, CblasLeft, CblasTrans);
    check_dtrsm(n, CblasRight, CblasNoTrans);
    check_dtrsm(n, CblasRight, CblasTrans);
    check_strmm(n, CblasLeft, CblasLower, CblasNoTrans);
    check_strmm(n, CblasLeft, CblasLower, CblasTrans);
    check_strmm(n, CblasRight, CblasLower, CblasNoTrans);
    check_strmm(n, CblasRight, CblasLower, CblasTrans);
    check_strmm(n, CblasLeft, CblasUpper, CblasNoTrans);
    check_strmm(n, CblasLeft, CblasUpper, CblasTrans);
    check_strmm(n, CblasRight, CblasUpper, CblasNoTrans);
    check_strmm(n, CblasRight, CblasUpper, CblasTrans);
    check_dtrmm(n, CblasLeft, CblasLower, CblasNoTrans);
    check_dtrmm(n, CblasLeft, CblasLower, CblasTrans);
    check_dtrmm(n, CblasRight, CblasLower, CblasNoTrans);
    check_dtrmm(n, CblasRight, CblasLower, CblasTrans);
    check_dtrmm(n, CblasLeft, CblasUpper, CblasNoTrans);
    check_dtrmm(n, CblasLeft, CblasUpper, CblasTrans);
    check_dtrmm(n, CblasRight, CblasUpper, CblasNoTrans);
    check_dtrmm(n, CblasRight, CblasUpper, CblasTrans);
  }
  if (g_fail) {
    fprintf(stderr, "%d checks failed\n", g_fail);
    return 1;
  }
  printf("check: GEMM/TRSM/TRMM odd sizes + transposes + LN/LT/RN/RT + U OK\n");
  return 0;
}
