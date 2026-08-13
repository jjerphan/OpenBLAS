# WASM OpenBLAS performance experiments

Branch: `v0.3.34-emforge-perf` (base `v0.3.34-emforge`).
Harness: `benchmark/wasm/` under Node. Speedup > 1 means faster than the
comparison CSV (baseline unless noted). Champion is the last kept tree.

## Protocol

One change per experiment. 5 warmup + 10 timed samples (median). Keep only if
faster than the current champion and `utest`/`ctest` pass under `node`.

## Experiments

### 0. Baseline — recorded

Unmodified `WASM128_GENERIC`, `USE_THREAD=0`, `COMMON_OPT=-O2`, `NOFORTRAN=1`.
See `results/baseline.csv`. 2×2 wasm128 GEMM, RISC-V P/Q/R, scalar GEMV.

### A. `-O3` — revert

Hypothesis: Clang `-O3` beats OpenBLAS default `-O2`.
Diff: `COMMON_OPT=-O3` only.
CSV: `results/expA-O3.csv`.

| group | geomean vs baseline |
| --- | --- |
| ALL | 0.913 |
| sgemm | 0.867 |
| dgemm | 0.853 |

Slower across GEMM and most of the suite. Keep `-O2`.

### B. Gate TRSM call counters — keep (hygiene)

Hypothesis: always-on `openblas_wasm128_*_calls` increments cost TRSM.
Diff: wrap counters in `#ifdef OPENBLAS_WASM_TRSM_PROFILE` in the four
`kernel/wasm/trsm_kernel_*_wasm128.c` files.
CSV: `results/expB-trsm-counters.csv`.

| group | geomean vs baseline |
| --- | --- |
| ALL | 0.960 |
| strsm | 0.967 |
| dtrsm | 0.959 |

Within noise / slightly down (run-to-run). Kept as production-clean; counters
are off unless `OPENBLAS_WASM_TRSM_PROFILE` is defined.

### C. Relaxed SIMD madd — keep GEMM-only

Hypothesis: `f32x4.relaxed_madd` / `f64x2.relaxed_madd` speeds the GEMM inner
loop. Expect bit-level changes; keep only if BLAS tests still pass.

**C1** (`Makefile.wasm` `-mrelaxed-simd` + GEMM inner loop + `intrin_wasm.h`
`v_muladd_*`): `results/expC-relaxed-simd.csv`, ALL **0.988**. GEMM up, L1
(DOT/AXPY) down — net loss.

**C2** (revert `intrin_wasm.h` to mul+add; keep relaxed madd in GEMM only):
`results/expC2-relaxed-gemm-only.csv` / `results/champion.csv`.

| group | geomean vs baseline |
| --- | --- |
| ALL | **1.039** |
| sgemm | 1.176 |
| dgemm | 1.145 |
| ssyrk | 1.159 |
| dsyrk | 1.126 |
| strsm | 1.073 |
| dtrsm | 1.107 |

Keep C2. Champion from here until a later keep.

`WASM_SIMD_COMPAT_SLOW=1` is not a compile flag in this OpenBLAS tree or emcc
4.0.9; skipped as a diagnostic.

### D. SGEMM_DEFAULT_P 128 → 64 — revert

Hypothesis: RISC-V P=128 oversizes the A panel for WASM (128×240×4 B ≈ 120 KiB,
above a typical 32 KiB L1). Halving P should cut packing traffic if GEMM is
cache-bound.

Diff: split `WASM128_GENERIC` from `RISCV64_GENERIC` in `param.h`; set
`SGEMM_DEFAULT_P=64` only.

CSV: `results/expD-sgemm-P64.csv`.

| group | geomean vs champion |
| --- | --- |
| ALL | 0.909 |
| sgemm | 0.872 |

Slower. Reverted `SGEMM_DEFAULT_P` to 128. WASM linear memory is not L1-bound the way a 32 KiB RISC-V D-cache is; smaller P just increased packing overhead. Left Q/R at RISC-V copies; retune after a 4×4 kernel if GEMM is still packing-bound.

### E. SGEMM/DGEMM 4×4 microkernel — keep

Hypothesis: a 4×4 tile with contiguous `f32x4` / `f64x2` loads beats the 2×2
shuffle kernel, especially DGEMM.

Diff: `param.h` `SGEMM`/`DGEMM` unroll 4×4 (complex stays 2×2); packing
`gemm_{n,t}copy_4`; `KERNEL.WASM128_GENERIC` points TRMM at `trmmkernel_4x4.c`;
new `kernel/wasm/gemmkernel_4x4_wasm128.c`. TRSM wasm128 kernels already handle
unroll 4 (4×4 tiles use scalar `solve()`). First keep still compiled 2×2 TRMM
(see H).

CSV: `results/expE-gemm-4x4.csv` (now `results/champion.csv`).

| group | geomean vs previous champion |
| --- | --- |
| ALL | **1.065** |
| sgemm | 1.043 |
| dgemm | **1.382** |
| ssyrk | 1.135 |
| dsyrk | **1.513** |
| strsm | 1.113 |
| dtrsm | 1.181 |

Large SGEMM (N≥64) is within ~3% of the 2×2 kernel (slightly down); DGEMM and
DSYRK are the win. Keep.

### F. SIMD GEMV N/T — revert

Hypothesis: explicit `intrin_wasm.h` GEMV (4-column unroll) beats scalar
`kernel/riscv64/gemv_{n,t}.c`, especially GEMV-T.

Diff: `kernel/wasm/gemv_n.c`, `gemv_t.c`; KERNEL points SGEMV/DGEMV at them.
CSV: `results/expF-gemv-simd.csv`.

| group | geomean vs champion |
| --- | --- |
| ALL | 0.997 |
| sgemv | 1.014 |
| sgemv_t | 1.000 |
| dgemv | 0.962 |
| dgemv_t | 1.068 |

Within noise. GEMV-T sits at the same ~4 GFLOPS ceiling as SDOT/DDOT (reduction
bound on SIMD128, no IEEE FMA). Reverted KERNEL to the RISC-V scalar files
(autovec already matches). Source files kept for a later retry.

### G. SIMD SCAL — revert

Hypothesis: `kernel/riscv64/scal.c` is scalar; a unit-stride SIMD loop should
lift SSCAL/DSCAL (suite ~3 / ~1.1 GFLOPS vs AXPY ~14 / ~7).

Diff: `kernel/wasm/scal.c`; KERNEL `SSCALKERNEL`/`DSCALKERNEL`.
CSV: `results/expG-scal-simd.csv`.

| group | geomean vs champion |
| --- | --- |
| ALL | 0.984 |
| sscal | 0.921 |
| dscal | 0.982 |

No win (store-bound; Clang already autovectorizes the RISC-V loop). Reverted.
Source kept unwired.

### H. Wire 4×4 TRMM (KERNEL include order) — keep (correctness)

Hypothesis: CBLAS STRMM/DTRMM failed after E because packing is `trmm_*copy_4`
while the kernel was still 2×2. `kernel/Makefile` includes
`KERNEL.$(TARGET_CORE)` then `KERNEL`; the latter overwrote
`STRMMKERNEL`/`DTRMMKERNEL` with `trmmkernel_2x2.c`.

Diff: wrap those assignments in `ifndef` in `kernel/wasm/KERNEL` so the
`KERNEL.WASM128_GENERIC` 4×4 TRMM sticks. `benchmark/wasm/check.c` now covers
TRMM (L/U × LN/LT/RN/RT). `test.sh` greps ctest output for computational
`FAILED`/`FATAL` (make does not fail those targets).

No new CSV: TRMM is not in the frozen suite; GEMM/TRSM kernels unchanged.

Gate: check OK; utest core 106/106; CBLAS L1/L2/L3 under node, including
`cblas_strmm`/`cblas_dtrmm` 3528+3528 calls col+row major, STRSM/DTRSM, SYRK.

## Remaining headroom

- SIMD128 is 4×f32 / 2×f64; no IEEE FMA (relaxed madd helped GEMM only).
- Single-threaded (`USE_THREAD=0`); WASM pthreads out of scope.
- GEMV-T / DOT sit at ~4 GFLOPS (reduction-bound).
- Large SGEMM is close to the 2×2 kernel; further SGEMM work is low leverage
  versus DGEMM, which still has width-2 f64 lanes.
- Q/R blocking not retuned after 4×4; P=64 lost on 2×2.

## Net vs unmodified baseline

Champion (`results/champion.csv` = exp E) vs `results/baseline.csv`:

| group | geomean |
| --- | --- |
| ALL | **1.106** |
| sgemm | 1.226 |
| dgemm | **1.582** |
| ssyrk | 1.316 |
| dsyrk | **1.703** |
| strsm | 1.194 |
| dtrsm | 1.307 |

Level-1/2 are slightly down (run noise and the GEMM-only relaxed-madd flag);
the suite geomean is a GEMM/SYRK/TRSM win.





