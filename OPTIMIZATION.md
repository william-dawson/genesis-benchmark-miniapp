# Experimental RIKYU optimization snapshot

This branch contains the optimized benchmark mini-app measured on RIKYU
B200 with NVHPC/HPC-X 26.3, one MPI rank and 32 OpenMP threads, in double
precision. It is based on upstream commit
`fc545c73b543bcdde9141bb2cfaa25ba1497c183` and exports production build
345258 from the FugakuNEXT benchmarking campaign, followed by the sparse
CPU exclusion-reset optimization in build 345704.

The snapshot includes tuned CUDA nonbonded kernels; resident order-4
NOOPT_2DALLTOALL PME using cuFFT; sparse GPU exclusion-mask construction
and pair-list compaction; GPU bonded-force production/assembly; the resident
VRES timestep, thermostat and coordinate producers; stable GPU outgoing
water/constrained-group migration; and parallel CPU incoming migration.
CPU fallbacks and opt-in diagnostic helpers remain. CPU exclusion rebuilding
restores the previous sparse zero records instead of repeatedly filling
dense masks. Set `GENESIS_EXCLUSION_RESET=1` on the existing one-rank CHARMM
GPU-mask path with domain mode 6 and no halos. The optimized run preset
enables it; setting it to 0 retains the original full initialization.

Pair-list work estimates and stable descending sorting also stay on the
GPU with `GENESIS_PAIR_SORT=1` (enabled by the optimized preset). The old
CPU algorithm and transfers remain available with flag 0. This requires
one rank and no halos on the existing GPU pair-list/nonbond path.

Build **345739**, matched unprofiled **345756**, two runs per arm on one
exclusive B200, seed 734381, 6,000 steps: parent 345704/new-disabled/
GPU-sort dynamics **15.7965/15.7815/15.6635 s**; candidate total 17.6055 s.
This saves **0.133 s / 0.84%** against the matched parent, removing about
117 MB of host/device count/index traffic per 600 steps. GPU sorting takes
~4.4 ms versus ~23 ms CPU calculation; GPU pair-list construction remains
~51 ms per 600 steps. These scopes are not additive whole-run savings.

Full **345755** passes 1,001 exact sort keys/count/stable-order comparisons,
all 1,001 exact exclusion/GPU-consumer checks, 6,002 force and energy/virial
checks and all existing migration/thermostat/COM/coordinate oracles. Max
normalized force/energy-virial errors: 5.917e-5/5.852e-4 (pass ≤1), finite
observables, no reported invalid/divide-by-zero/overflow. Separate 345744
passes with CPU counts/keys/order poisoned at all 101 short-run rebuilds.
Frozen initial/evolved 345757 passes the unchanged investigated-startup
policy and controls. All 77 previous GPU kernels retain identical PTX;
ten CUDA C++ objects match between performance/check builds. The existing
strict integration-state caveat remains unresolved.

`GENESIS_PAIR_SORT_VERIFY=1` runs the original CPU algorithm from downloaded
counts and compares all active results byte for byte. The optimized
verification preset enables it. `GENESIS_PAIR_SORT_POISON=1` separately
checks the production path without valid CPU mirrors. Details and theory:
companion `report/genesis-domain-tuning/PAIR_SORT_RESIDENCY.md`.

## Build and run

Build **346414** adds cached GPU outgoing bonded migration with
`GENESIS_BONDED_MIGRATION_GPU=1` (enabled by the optimized preset; flag0
retains the original CPU routines). It batches bonds, angles and dihedral
terms, caches fixed parameters by unique interaction ID, and performs
classification, stable compaction/sorting and payload gathering on GPU.
CPU scatter and incoming payload updates remain; the private identity grid
follows the original incoming permutation. Scope: one-rank, no active halos,
fixed-topology CHARMM double precision on the existing GPU path.

Final matched **346423**, four/arm, 6,000 steps, seed 734381: disabled/GPU
dynamics **15.7615→15.5760 s**, saving **0.1855 s / 1.18%**; total
17.7125→17.5430 s. GPU calculation is ~2.7 ms/600 versus ~44 ms for the
old outgoing CPU scan/packing, excluding transfers, scatter and first-use
cache costs. The comparable complete outgoing+incoming diagnostic falls
58.356→48.753 ms/600. No large application speedup is claimed.

Full correctness **346421** passes all 1,000 exact outgoing payload/order
and incoming identity/conservation events for each of the three kinds,
6,002 force and energy/virial checks, and all existing migration,
thermostat, COM, coordinate, exclusion and pair-sort gates. Max normalized
force/energy-virial errors are6.589e-5/4.641e-4 (pass≤1); observables finite,
seed 734381, no reported invalid/divide-by-zero/overflow. Frozen initial/
evolved **346422** passes the unchanged investigated-startup policy and
controls. Existing strict integration-state acceptance caveat remains.
All84 prior PTX entries and shared CUB bodies stay unchanged; all 11 CUDA
C++ objects match between production/check builds. New classification/
gather kernels have no stack/local traffic after explicit descriptor
selection. Detailed equations, failed prototypes, phase attribution,
PTX/counters and reproduction: companion
`report/genesis-domain-tuning/BONDED_MIGRATION.md`.

`GENESIS_BONDED_MIGRATION_VERIFY=1` enables exact outgoing and incoming-ID
checks; `GENESIS_OPT_VERIFY=1` includes them in the preset. Verification
runtime is excluded from throughput claims.

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

## Sparse-history exclusion reset

Build **345704** preserves all 77 GPU PTX entries. Diagnostic job 345703,
three repeats per arm, measures **146.934→35.729 ms** per 600 steps for
exclusion rebuilding (4.11×). Initialization/restoration itself falls
107.827→5.267 ms. This changes CPU work only, adding no GPU transfers.
The remaining 1–2/1–3 and 1–4 construction code retains its original order.

Unprofiled job **345722**, two runs per arm, compares the published 003bbed
binary with this build on one exclusive B200, 6,000 steps, seed 734381:

| Path | Dynamics median (range), seconds | Total median, seconds |
|---|---:|---:|
| Published 003bbed | 17.423 (17.406–17.440) | 19.549 |
| New build, reset off | 17.685 (17.472–17.898) | 19.6255 |
| Sparse reset | 15.738 (15.738–15.738) | 17.6595 |

The gain against the published checkpoint is **1.685 s / 9.67% dynamics**.
The slower disabled-control repeat is retained. The full-run saving exceeds
the directly timed exclusion leaf saving; other CPU phases also improve,
consistent with reduced cache pollution. Do not attribute it all to reset.

Independent full-length job **345727** passes **1,001 exact mask and ordered
list/scale comparisons**, plus 1,001 GPU mask/packing and atom-set checks,
6,002 force and energy/virial comparisons, all 1,000 events of each migration
oracle, 999 thermostat events, 3,996 scales, 599 COM events and 5,999 exact
coordinate checks. Maximum normalized force/energy-virial errors are
6.070e-5 / 4.907e-4 (pass ≤1). All 100 observable rows are finite, both
printed seeds are 734381, and no invalid/divide-by-zero/overflow is reported.
Frozen initial/evolved gates **345723** pass the unchanged investigated
MPI-startup policy and controls. The older strict integration-state
acceptance caveat above remains unresolved.

`GENESIS_EXCLUSION_VERIFY=1` saves candidate results, repeats the original
dense construction, and compares active masks, counts, ordered exclusion/
1–4 lists and scales byte for byte. `GENESIS_OPT_VERIFY=1` enables it and
the downstream GPU-mask checks. Verification time is excluded from timing
claims. Detailed equations, raw-result locations and figure reproduction:
companion report `report/genesis-domain-tuning/EXCLUSION_REBUILD.md`.
