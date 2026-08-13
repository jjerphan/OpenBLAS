#!/usr/bin/env bash
# Correctness gate: CBLAS odd-size GEMM/TRSM check, utest (no Fortran-ABI
# gemv tests), and CBLAS ctest under node.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/benchmark/wasm/out}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
COMMON_OPT="${COMMON_OPT:--O2}"

activate_emscripten() {
  if command -v emcc >/dev/null 2>&1; then
    return 0
  fi
  local prefix="${OPENBLAS_EM_PREFIX:-}"
  if [[ -z "$prefix" ]]; then
    local cand
    for cand in \
      "$ROOT/.em-prefix" \
      /home/jjerphan/dev/recipes/output/bld/rattler-build_openblas_*/build_env
    do
      if [[ -f "$cand/bin/activate_emscripten.sh" ]]; then
        prefix="$cand"
        break
      fi
    done
  fi
  if [[ -z "$prefix" || ! -f "$prefix/bin/activate_emscripten.sh" ]]; then
    echo "emcc not found. Set OPENBLAS_EM_PREFIX to an emscripten-forge env." >&2
    exit 1
  fi
  export CONDA_PREFIX="$prefix"
  export PREFIX="${PREFIX:-$prefix}"
  export PATH="$prefix/bin:$prefix/opt/emsdk/upstream/emscripten:$prefix/opt/emsdk/upstream/bin:$PATH"
  # shellcheck disable=SC1091
  source "$prefix/bin/activate_emscripten.sh"
}

activate_emscripten
cd "$ROOT"

LDFLAGS="${LDFLAGS:--O2 -msimd128 -mrelaxed-simd -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=8MB -sINITIAL_MEMORY=512MB -sEXIT_RUNTIME=1}"

MAKE_COMMON=(
  HOSTCC=gcc
  CC=emcc
  AR=emar
  RANLIB=emranlib
  TARGET=WASM128_GENERIC
  USE_THREAD=0
  NOFORTRAN=1
  NO_LAPACK=1
  COMMON_OPT="$COMMON_OPT"
  LDFLAGS="$LDFLAGS"
  CEXTRALIB="$LDFLAGS"
  CROSS=0
  COLORCODE=0
)

LIB=$(ls -1 "$ROOT"/libopenblas_wasm128*.a | head -n1)
if [[ -z "$LIB" ]]; then
  echo "libopenblas_wasm128*.a not found; run build.sh first" >&2
  exit 1
fi

mkdir -p "$OUT"

echo "==> CBLAS GEMM/TRSM odd-size check"
emcc $LDFLAGS -I"$ROOT" -o "$OUT/check.js" "$ROOT/benchmark/wasm/check.c" "$LIB"
node "$OUT/check.js"

echo "==> utest core (skip test_gemv + extensions: Fortran hidden char-length ABI on this branch)"
# $(UTESTBIN) is .PHONY, so build via run_test's first binary only.
UTEST_OBJS="utest_main.o test_min.o test_amax.o test_ismin.o test_rotmg.o test_axpy.o test_dotu.o test_dsdot.o test_swap.o test_rot.o test_dnrm2.o test_zscal.o test_amin.o test_axpby.o"
emmake make -C utest -j"$JOBS" "${MAKE_COMMON[@]}" $UTEST_OBJS
emcc $LDFLAGS -o "$ROOT/utest/openblas_utest" \
  $(printf " $ROOT/utest/%s" $UTEST_OBJS) \
  "$LIB"
node "$ROOT/utest/openblas_utest"

echo "==> ctest"
CTEST_LOG="$OUT/ctest.log"
# CBLAS testers print FAILED/FATAL and continue; make does not fail the target.
set +e
emmake make -C ctest all -j"$JOBS" "${MAKE_COMMON[@]}" 2>&1 | tee "$CTEST_LOG"
ctest_make=$?
set -e
if [[ "$ctest_make" -ne 0 ]]; then
  echo "ctest make failed (exit $ctest_make)" >&2
  exit 1
fi
if grep -E 'FAILED ON CALL|FAILED  -|COMPUTED RESULT IS LESS THAN HALF ACCURATE' "$CTEST_LOG" >/dev/null; then
  echo "ctest reported computational failures (see $CTEST_LOG)" >&2
  grep -E 'FAILED ON CALL|FAILED  -|COMPUTED RESULT IS LESS THAN HALF ACCURATE|cblas_s|cblas_d|cblas_c|cblas_z' "$CTEST_LOG" | head -n 40 >&2
  exit 1
fi

echo "correctness gate passed"
