#!/usr/bin/env bash
# Build OpenBLAS for wasm32 + the Node benchmark driver.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/benchmark/wasm/out}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
COMMON_OPT="${COMMON_OPT:--O2}"
EXTRA_CFLAGS="${EXTRA_CFLAGS:-}"
# Keep EXTRA_CFLAGS on COMMON_OPT so we do not override Makefile.wasm's CCOMMON_OPT += -msimd128.
COMMON_OPT="$COMMON_OPT $EXTRA_CFLAGS"

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
  # shellcheck disable=SC1091
  export CONDA_PREFIX="$prefix"
  export PREFIX="${PREFIX:-$prefix}"
  export PATH="$prefix/bin:$prefix/opt/emsdk/upstream/emscripten:$prefix/opt/emsdk/upstream/bin:$PATH"
  # shellcheck disable=SC1091
  source "$prefix/bin/activate_emscripten.sh"
}

activate_emscripten

mkdir -p "$OUT"
cd "$ROOT"

# Drop leftover host objects (a previous native build in this tree will
# otherwise get archived into libopenblas_wasm128*.a).
if [[ "${WASM_BENCH_CLEAN:-1}" != "0" ]]; then
  echo "==> make clean (drop host objects leftover from native builds)"
  emmake make clean COLORCODE=0 >/dev/null || true
  find "$ROOT" -name '*.o' -not -path '*/.git/*' -delete
  rm -f "$ROOT"/libopenblas_wasm128*.a "$ROOT"/libopenblas.a
fi

echo "==> OpenBLAS libs TARGET=WASM128_GENERIC USE_THREAD=0 COMMON_OPT=$COMMON_OPT"
emmake make libs \
  -j"$JOBS" \
  HOSTCC=gcc \
  CC=emcc \
  AR=emar \
  RANLIB=emranlib \
  TARGET=WASM128_GENERIC \
  USE_THREAD=0 \
  NOFORTRAN=1 \
  NO_LAPACK=1 \
  COMMON_OPT="$COMMON_OPT" \
  COLORCODE=0

LIB=$(ls -1 "$ROOT"/libopenblas_wasm128*.a | head -n1)
if [[ -z "$LIB" ]]; then
  echo "libopenblas_wasm128*.a not found" >&2
  exit 1
fi

LINKFLAGS=(
  -O3
  -msimd128
  -mrelaxed-simd
  -sALLOW_MEMORY_GROWTH=1
  -sSTACK_SIZE=8MB
  -sINITIAL_MEMORY=256MB
  -sEXIT_RUNTIME=1
  -sASSERTIONS=0
)

echo "==> linking $OUT/bench.js against $LIB"
emcc "${LINKFLAGS[@]}" \
  -I"$ROOT" \
  -o "$OUT/bench.js" \
  "$ROOT/benchmark/wasm/bench.c" \
  "$LIB"

echo "built $OUT/bench.js"
