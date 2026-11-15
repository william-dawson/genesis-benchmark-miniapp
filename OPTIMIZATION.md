# Experimental RIKYU optimization snapshot

Latest unprofiled timing for build351026 (one B200, one MPI rank x32 CPU
threads, double ApoA1,92,224atoms,6,000steps): four-run medians are13.553s
total and11.645s dynamics; a separate four-run confirmation gives13.5325s
and11.6295s. These application timers overlap and must not be summed.

An exclusive partition of trace351040's540warm steps (1,158.279ms) gives
65.5% GPU kernels active,3.6% copies/clears without kernels, and30.9% GPU
idle. Idle time is9.3% in copy APIs,4.2% in synchronization APIs,5.1% in
launch/other CUDA APIs and12.3% outside CUDA APIs; this is not a measurement
of CPU arithmetic or fully recoverable speedup. Explicit warm-window traffic
is76.359MB host-to-device,41.171MB device-to-host and14.287GB device-to-device.
These byte totals exclude coherent GPU stores to host metadata arrays.
The nonbonded GPU kernels are the largest exclusive kernel family (32.0%
of wall time); another6.1% has multiple kernel families active concurrently.
GPU kernel tuning and the slow device-copy population are now primary
investigations alongside measured removal of remaining host boundaries.

The [accepted report, additive ledger and reproducible analysis](https://github.com/RIKEN-RCCS/block1-eeas/tree/b185d81/report/genesis-observables)
record correctness evidence, full timing sets and the100% wall breakdown.

Build **351026** keeps periodic INFO reporting resident: reduce squared
thermostat-snapshot forces and mass-weighted pre-update reference velocities
on GPU, return two doubles, and retain the existing scalar report formulas.
The preset selects `GENESIS_OBSERVABLES_GPU=1`;0 keeps full report sync,
2 compares the original CPU sums on identical inputs,3 poisons retired CPU
reporting mirrors. Output files and startup/final/unsupported paths retain
their old synchronization. Full351037 passes99reports/reference arm and all
previous numerical/ownership gates; maximum relative difference2.11e-14.
Normal/tight351052/351053 and strict frozen351051 pass; I/O351054 produces
complete finite coordinate/velocity frames and restart files. All159old PTX
bodies match350496; two reduction kernels are the only additions.
Traces351039/351040 remove all81large periodic report downloads per540warm
steps: D2H605.924→41.171MB. Full matched351038/351084 improve dynamics
2.34%/2.29% (~0.27–0.28s per6,000steps). Evidence:
`report/genesis-observables/README.md` in the companion knowhow repository.

Build **350496** closes the maintenance coordinate/velocity round trip and
translation-pack host handoff for the one-GPU ApoA1 path. Resident migration
already updates the live integration allocations: keep those values, retire
their maintenance downloads/reloads and publish only other required metadata.
Configure coordinate packing from resident translations/charge/counts/offsets;
retire host translation publication. Real output and fallbacks retain sync.
The preset selects `GENESIS_STATE_CYCLE=1` /`GENESIS_PACK_RESIDENT=1`;0 keeps
old paths,2 retains independent CPU mirrors and checks exact active values,
3 poisons retired CPU mirrors. Production retirement excludes CPU oracles.
Full350560 passes1,000rebuilds/arm and every original physical/ownership gate.
Normal/tight350565/350566 and strict frozen350564 pass. All157old PTX bodies
match350411; two new metadata publishers are the only additions.
Traces350562/350563 remove96.48% H2D and67.44% D2H bytes per540warm steps;
large remaining state downloads occur only at output cadence. Two interleaved
full timing sets350561/350591 reduce dynamics4.37%/4.47%, total0.601/0.580s.
This closes the two boundaries, not all remaining host work. Evidence:
`report/genesis-state-translation/README.md` in the companion knowhow repository.

Build **350411** replaces nine maintenance metadata uploads with direct
HBM copies from the current incoming/boundary views. The optimized preset
selects `GENESIS_METADATA_IMPORT=1`;0 retains uploads,2 checks every active
field byte against CPU metadata,3 poisons device destinations before overwrite.
The path requires the completed resident maintenance handoff; other paths
retain their uploads. Full350418 passes1,000events in each of four arms,
all existing force/energy/virial and ownership gates, and unchanged budgets.
Frozen350430 passes strict initial/advanced frames. All157retained kernel
PTX bodies match350329. Interleaved full sets350419/350434 improve dynamics
3.15%/2.61%, saving0.430/0.345s total. Ordinary traces350427/350428 remove
2.425GB (52.8%) of H2D traffic per540warm steps; D2H is unchanged.
Full residency remains incomplete: next are state and translation transfers.
Companion evidence: `report/genesis-maintenance-import/README.md`.

Build **349465** retires the CPU raw nonbonded-force mirror on CHARMM
short steps whose GPU force assembly reads the device result directly.
`GENESIS_NONBOND_MIRROR=1` is enabled in the optimized preset; 0 retains
downloads; 2 poisons retired mirrors and checks direct consumption.
Actual CPU virial/reference readers keep downloads. Both full matched
sets 349470/349478 show **0.52–0.60% shorter dynamics** against same-build
controls (13.364→13.2945 s, 13.351→13.2705 s). First-set ranges overlap;
net gain against the previous executable is 0.24–0.78%, so do not claim a
large arithmetic speedup. Trace 349469 proves 300 copies / 664,012,800 B
removed per 600 steps, H2D/D2D unchanged. Full 349472 passes 3,000 poisoned
mirrors/direct consumers and all previous numerical/ownership gates;
frozen 349471 passes both policies. All 97 device bodies stay identical;
17 performance/check CUDA C++ objects match. Full-step CPU virial remains.
Companion evidence: `report/genesis-nonbond-virial/README.md`.

For integrator/constraint changes or unexplained state/constraint-virial
differences, compare normal and 100× tighter SHAKE tolerances (1e-10 and
1e-12 Å here), using identical inputs at each CPU/GPU operator boundary,
fixed seed and unchanged error budgets. Keep diagnostic runs separate
from normal-tolerance performance. Analytical SETTLE is independent of
SHAKE tolerance; a diagnostic run finishing does not make recorded budget
failures pass. Current check 349456 passes all position/velocity-state
comparisons in both arms but retains 16/500 and 11/500 constraint-virial
budget failures. Details: `report/genesis-timestep-port/current-shake-check/README.md`.

Build **349387** prepares remapped GPU state before force evaluation.
Coordinates, velocities and maps upload once, then the existing periodic
pack publishes with unit scale. Bonded raw/periodic inputs and PME reuse
that state; late VV2 prepare avoids a second import. Actual CPU virial/
reference readers retain current host data. `GENESIS_REMAP_PREPARE=1`
enables the path; 0 keeps late prepare; 2 checks active coordinates,
velocity, translations and both packed layouts bit-for-bit, poisons packed
GPU outputs before overwrite, and requires GPU consumers before VV2.

Matched full **349392/349409**, four primary runs per arm/set, gives
late/early dynamics **14.6040/13.3875 s** and **14.4305/13.2790 s**:
**7.98–8.33%** dynamics and **7.28–7.61%** total gain against same-build
controls. Previous-executable controls give net **7.3–8.2%** dynamics gain;
confirmation's disabled arm is ~0.75% slower than its parent. Traces move
310 copies / 1.310 GB per 600 steps from H2D to D2D, with 100 additional
invocations of unchanged periodic packing. Full **349393** passes 1,000
instances of all new exact/poison/consumer gates and every previous
numerical gate. Frozen **349389** passes both policies. All 97 old PTX
entries and 16 performance/check CUDA C++ objects match. The separate
strict integrator-state caveat remains unresolved. Presets select 1/2.
Companion evidence, allocation-heavy trace anomaly and full controls:
`report/genesis-remap-prepare/README.md`.

Build **349289** trims maintenance synchronization to coordinates and
velocities. Migration does not consume the five omitted work/reference/
force fields; GPU VV2 regenerates scratch and uses freshly published
forces. Velocity-reference bytes remain on GPU instead of making an
unchanged round trip. `GENESIS_MAINT_STATE=1` enables it, 0 keeps full sync,
2 privately verifies retained reference bytes, poisons five CPU mirrors,
checks that CPU reference fields remain unwritten and poisons GPU work
before regeneration. Actual CPU output/fallbacks keep full state sync;
independent integrator-state reference mode keeps the old handoff.

Full matched **349331/349336**, four primary runs per arm in each set:
full/trimmed dynamics **15.0090/14.5575 s** and **14.9695/14.5500 s**:
**2.80–3.01%** dynamics and **2.47–2.65%** total-time reduction, saving
**0.42–0.45 s** per 6,000 steps. Trace removes 500 D2H +100 H2D slab copies /
4.183 GB per 600; full validated event counts imply **41.83 GB** less
traffic. Full **349332** passes 1,000 instances of all four new ownership
checks plus all existing numerical gates. Frozen **349327** passes both
policies; all 97 prior PTX entries and 15 performance/check CUDA C++ objects
match. The earlier strict integrator-state caveat remains unresolved.
The optimized/verification presets select modes 1/2. Companion evidence:
`report/genesis-maintenance-state/README.md`.

Build **349243** retires GPU-owned translations before maintenance
recomputes every active value, and skips translation refreshes for ordinary
output that does not read them. `GENESIS_TRANSLATION_REBUILD=1` enables it;
0 keeps the prior refresh path; 2 poisons old/output mirrors and compares
all rebuilt active values bit-for-bit with the original boundary routine.
Rebuild intervals forbid real consumers or cache configuration until
completion. CPU fallback/reference consumers retain mirror refreshes.
The optimized/verification presets select modes 1/2, respectively.

Matched full sets **349248/349257**, four primary runs per arm per set,
give refresh/retire dynamics medians **14.9555/14.9030 s** and
**14.9575/14.8870 s**: **0.35–0.47%** shorter dynamics. Total median savings
are 14/54 ms. Production traces eliminate all 109 translation downloads per
600 steps; the full production path loses 1,099 / **7.66 GB** of downloads.
CPU maintenance still rebuilds translations and H2D cache refill remains.
Full **349249** passes 1,000 poisoned exact overwrite events, 99 output
poisons, all existing numerical gates; frozen **349245** passes both policies.
All 97 prior PTX entries remain unchanged. The separate strict integration-
state caveat remains unresolved. Companion evidence:
`report/genesis-translation-rebuild/README.md`.

Build **349165** makes the existing GPU periodic pack the translation
owner. `GENESIS_TRANSLATION_GPU=1` removes the duplicate CPU scaling loop;
0 retains it; 2 compares every active update bit-for-bit with the original
Fortran CPU scaling sequence and poisons the live host mirror. Actual
coordinate/state downloads refresh translations, and CPU consumer/cache
reconfiguration/invalidation boundaries reject stale mirrors. The optimized
preset enables mode 1 and verification selects mode 2. No new kernel or
launch is introduced; all 97 prior normalized PTX entries are unchanged.

Two matched 6,000-step sets **349175/349183**, four primary runs per arm
in each, give CPU/GPU-owned dynamics medians **15.0535/15.0015 s** and
**15.0620/14.9905 s**: a small **0.35–0.47%** gain. Total wall-time change
is not consistent. Refreshes add ~7.66 GB of padded D2H per full run,
limiting the benefit. Full **349176** passes 5,999 exact scaling/host-poison
checks, all 6,002 force/energy/virial and ordered totals, and every existing
component gate. Frozen **349174** passes both policies. The old separate
strict integrator-state caveat remains unresolved. Companion evidence:
`report/genesis-translation-residency/README.md`.

The separate `GENESIS_FORCE_CLEAR` experiment remains disabled (default 0):
although correct, it did not demonstrate a whole-run gain over build 348347.

Build **348347** keeps PME long forces on GPU through integration,
removing their download, CPU scatter/clears and reupload. After one
bootstrap handoff, PME writes the registered integrator destination
directly. `GENESIS_LONG_FORCE=1 GENESIS_LONG_DIRECT=1` enables the path;
long-force mode 2 checks every active value bit-for-bit against CPU scatter
and poisons the CPU mirror. The optimized/verification presets select
these modes. Startup, output, maintenance and state-reference fallbacks
retain CPU mirrors; scope is the existing one-rank packed GPU PME path.

Matched full **348351**, four primary runs/arm, gives same-build CPU/direct
dynamics **15.3190→14.9935 s** (2.12%), total 17.3080→16.9550 s. Production
trace 348356 removes 3.449 GB of host/device traffic per 600 steps. Full 348352
passes 3,000 exact publication/poison checks and 6,002 force/energy/virial
and ordered-total checks; frozen348354 passes strict and unchanged startup
policy locally. All 97 prior PTX entries match, as do 13 CUDA C++ objects
between production/check builds. The separate strict integration-state
caveat remains. Details: companion `report/genesis-long-force-residency/README.md`.

Build **348320** moves final bonded energy/virial accumulation to GPU,
with cooperative staging preserving the CPU term/cell summation order.
`GENESIS_BONDED_TOTALS=1` enables it; 0 retains CPU accumulation, and 2
also checks every per-worker result bit-for-bit. The optimized preset
enables mode1; its verification preset uses mode 2.
`GENESIS_BONDED_TOTALS_STAGED=0` selects the first direct-load GPU helper.

Matched full **348324**, four primary runs/arm, same build disabled/GPU:
dynamics **15.5110→15.2765 s** (1.51%), total 17.4495→17.2260 s. The helper
is 33.408→6.976 µs after staging, and removes about 438 MB of D2H traffic
per 600 steps. Full **348325** passes 6,002 exact ordered accumulations and
force/energy/virial shadows, plus all existing component gates. Frozen
**348329** passes strict initial/evolved checks and the unchanged startup
policy. All 94 pre-existing PTX entries stay unchanged. The separate strict
integration-state caveat remains. Evidence and reproduction are in the
companion `report/genesis-force-finalization/README.md`.

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
