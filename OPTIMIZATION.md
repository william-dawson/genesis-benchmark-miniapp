# Experimental RIKYU optimization snapshot

This branch contains the optimized benchmark mini-app measured on RIKYU
B200 with NVHPC/HPC-X 26.3, one MPI rank and 32 OpenMP threads, in double
precision. It is based on upstream commit
`fc545c73b543bcdde9141bb2cfaa25ba1497c183` and exports production build
345258 from the FugakuNEXT benchmarking campaign.

The snapshot includes tuned CUDA nonbonded kernels; resident order-4
NOOPT_2DALLTOALL PME using cuFFT; sparse GPU exclusion-mask construction
and pair-list compaction; GPU bonded-force production/assembly; the resident
VRES timestep, thermostat and coordinate producers; stable GPU outgoing
water/constrained-group migration; and parallel CPU incoming migration.
CPU fallbacks and opt-in diagnostic helpers remain. The next exclusion
reset experiment is not included.

## Build and run

On RIKYU:

```bash
module load nvhpc-hpcx/26.3
bash bootstrap
cuda_root="$(dirname "$(dirname "$(readlink -f "$(command -v nvcc)")")")"
bash configure FC=mpif90 CC=mpicc --enable-gpu \
  --with-cuda="$cuda_root" --with-gpuarch=sm_100
make -j16 PME_CUFFT_LIBS=-cudalib=cufft
```

The added CUDA Fortran kernels explicitly target `cc100`. Other GPU
architectures, precision modes, force fields and multi-rank residency paths
have not been accepted by these results. The original nonresident paths
remain available by leaving the optimization flags disabled.

The benchmark data are external, from input revision
`498091f5ce62291e0c7e10828da31d86545b5174`. Use ApoA1's
`npt/genesis2.0beta_3.5fs/apoa1/p8.inp` with the original directory structure,
changing `domain_x`, `domain_y`, `domain_z` to 1 and adding
`pme_scheme = NOOPT_2DALLTOALL` in ENERGY. Keep the full 6,000 steps. For
controlled accuracy comparisons set DYNAMICS `iseed = 734381`: the default
clock-based seed and incompatible eight-rank restart RNG state otherwise
give different random streams at one rank.

Run within a one-GPU allocation with 32 CPU cores. From the input directory:

```bash
bash /path/to/checkout/benchmarks/run-rikyu-optimized.sh \
  /path/to/checkout/src/spdyn/spdyn p1-long.inp
```

The preset exports the exact production settings used in the campaign.
Use `GENESIS_OPT_VERIFY=1` with the same command for the independent CPU
force/energy/virial, exact migration, thermostat and coordinate-pack checks.
Verification is much slower and must be excluded from timing.

## Measured results and numerical scope

Fixed-seed 6,000-step ApoA1 comparison, two runs/arm on one exclusive B200
(job 345319; average application timer column):

| Version | Dynamics median | Total including setup |
|---|---:|---:|
| Untuned upstream campaign control | 40.817 s | 42.7755 s |
| Previous outgoing-water version | 19.209 s | 21.1660 s |
| This snapshot | **17.606 s** | **19.5475 s** |

The complete campaign achieves 2.32× faster dynamics on this case. These
are benchmark-mini-app results, not full GENESIS application predictions.

Full-length verification job 345320 passes 6,002 same-input force and
energy/virial comparisons, 1,000 exact events for each outgoing/incoming
water/constrained-group operation, 999 thermostat preparations, 3,996
scales, 599 COM checks and 5,999 exact periodic-coordinate checks.
All 100 printed observable rows are finite; no invalid/divide-by-zero/
overflow flag is reported. Frozen initial/120-step force gates pass with
repeated controls (345282). Migration checks include payload bytes, global
IDs, stable record order, counts, topology lists and reverse maps.

The GPU integrator's separate strict state gate still has small investigated
SHAKE/SETTLE differences near the stopping tolerance. Default maximum solute
differences observed in that study were 7.43e-11 Å in position and 1.038e-9
internal velocity (about 2.12e-8 Å/ps). A 100× tighter diagnostic solve
shrinks those differences about 100×. These findings explain the observed
cases but do not establish complete trajectory equivalence. The strict gate
has not been relaxed or reclassified as passed; `GENESIS_STEP_VERIFY=2`
is diagnostic collection, not acceptance. The run preset leaves it disabled.

The inherited mini-app omits essential physics modules; the warning in
the main README continues to apply. This experimental optimization branch
is for benchmarking and port development, not scientific simulation.

## Provenance

`OPTIMIZED_SOURCE.sha256` identifies the source files changed from the
upstream base; generated/preprocessed files, objects, binaries and profiler
reports are excluded. The production source was recovered from the complete
build-tree archive, SHA256
`d3685bd2ff22ce43b1c8920a958f8850654cae9bf79e33f60f633f2b46718627`,
with only the separately applied frozen-validation driver removed.
Opt-in migration oracles and maintenance timers remain in the source.

Historical binary SHA256:

```text
performance 9f8dca59e07c9a89e60d26e8a24b2ed3b8c1dd517b17c9f171a8dcb9ee38b161
validation  710b20038b6a91e22b86f0eb8bfa1004d16c5ff78e22e40dca0fc8122e6a8beb
```

Detailed equations, data and validation methodology live in the companion
benchmarking repository's `report/genesis-domain-tuning/` reports,
especially `LONG_RUN.md`, `INCOMING_AND_GROUPS.md`, and the PME, pair-list
and timestep reports they reference. This snapshot preserves all 77
compiled PTX entries of the latest production source.

Clean-checkout validation of this export (same compiler/GPU configuration):
build **345579** verifies all 34 changed source hashes and rebuilds from
scratch. All 77 normalized PTX entries match build 345258. Two full-length
production runs **345607** measure 17.625/17.519 s dynamics and
19.598/19.479 s total, seed 734381. Independent full-length check **345608**
passes all 6,002 force/energy/virial comparisons, 1,000 events for each
migration oracle, thermostat/COM checks and 5,999 exact periodic-coordinate
checks. Maximum normalized force and energy/virial errors are
5.979e-5 / 5.211e-4 (pass ≤1). No invalid/divide-by-zero/overflow is reported.
The clean rebuilt binary SHA256 is
`d626af9ee2d3e4e769d2bb5c806c3fbe3d8ad57abc992452df196dafa6a52a81`.
