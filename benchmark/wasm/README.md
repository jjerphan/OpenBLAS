# WASM / Emscripten OpenBLAS microbenchmarks

Node-based harness for `TARGET=WASM128_GENERIC`. Warmup discards JIT, then
reports the median of timed samples.

```bash
# OPENBLAS_EM_PREFIX points at an emscripten-forge env that provides emcc.
./benchmark/wasm/build.sh
python3 benchmark/wasm/run.py -o benchmark/wasm/results/latest.csv
python3 benchmark/wasm/compare.py benchmark/wasm/results/baseline.csv \
    benchmark/wasm/results/latest.csv
```

`COMMON_OPT` defaults to `-O2` (OpenBLAS default). Override for experiments:

```bash
COMMON_OPT=-O3 ./benchmark/wasm/build.sh
```

After a keep, run correctness under Node:

```bash
./benchmark/wasm/test.sh
```

`test.sh` runs the odd-size CBLAS check, utest core (no Fortran-ABI `test_gemv`),
and CBLAS `ctest`. It fails if ctest prints computational `FAILED`/`FATAL`
(those testers do not fail the make target).
