# One-time hash table for Circuit ORAM batch stash access

Status: implementation plan

## Summary

`CircuitORAM::ORAM::BatchReadAndRemove` currently scans all `stashSize`
blocks for every request. The existing default is 33, but the intended batch
eviction scheme has a heavier stash tail and needs a larger bound.

The proposed optimization is sound if the stash is converted, for the
duration of one `BatchReadAndRemove` call, into a freshly keyed one-choice
bucketed hash table. A request always probes:

1. its four-entry main bucket, and
2. a shared four-entry overflow bucket.

This reduces the online stash work from the new stash bound to 8 entries. The
table must be built obliviously, use a fresh key for every batch, handle
duplicate requests without revealing equality, and be obliviously converted
back to the normal flat stash before `BatchReadAndRemove` returns.

The intended batch scheme performs no eviction on the read path and performs
two deterministic evictions on two different paths. The matching checked-in
trace is `circuit_oram_stash_original_no_read_two_paths.log`. Its fitted tail
gives a provisional 46-slot stash bound at `2^-64` per query, rather than 33.

The hash table is built only once per batch. For a batch of `q` queries, a
`2^-64` per-query goal permits this batch event probability to be as large as
`q * 2^-64`. Thus `2^-58` is sufficient when `q >= 64`, before any split among
other failure sources. With `S=46`, bucket size 4, overflow size 4, and one
candidate, the empirical convolution reaches approximately `2^-58.04` at 72
main buckets. A fitted-tail upper-envelope calculation reaches `2^-58.04` at
74 buckets; 80 buckets gives approximately `2^-58.94` and is a better
engineering starting point. This corresponds to a worst-case nominal main
table load factor of about 0.144. These figures remain experimental estimates,
not proved cryptographic bounds.

## Current behavior and integration point

The relevant code is in `omap/odsl/circuit_oram.hpp`.

- `path[0:stashSize]` is the persistent flat stash.
- Both in-memory and `DISK_IO` branches of `BatchReadAndRemove` call
  `ReadElementAndRemoveFromPath` on the entire stash once per request.
- Tree buckets are then searched on the requested path.
- `BatchReadAndRemove` returns with requested blocks removed. A later
  `BatchWriteBack` inserts updated blocks and runs Circuit ORAM eviction.
- `evictPath` and `WriteNewBlockToPath` assume that the first `stashSize`
  entries of `path` contain the authoritative flat stash.
- UIDs are sorted. `deDuplicatePoses` gives later occurrences of a duplicate
  UID random tree positions, and `duplicateVal` propagates the first result.

Consequently, a temporary table cannot merely contain copies while the old
stash remains authoritative: removing a block only from a copy would leave a
duplicate for writeback. Nor can the bucketed layout be left active after the
method returns, because the current writeback and eviction code scans and
mutates a flat `stashSize`-entry region.

The required lifecycle is therefore:

```text
flat stash snapshot
        |
        v
fixed-work oblivious table build
        |
        v
8-entry stash probe for every request + existing tree-path probe
        |
        v
oblivious compaction of remaining table entries
        |
        v
restored flat stash -> existing BatchWriteBack and eviction
```

The hash-table optimization initially applies only to `BatchReadAndRemove`.
Single-item `Read` and `Update` remain unchanged, and the hash layout itself
must not alter eviction semantics. The separate batch-writeback mismatch noted
below still needs to be resolved so the implementation matches the eviction
scheme being sized.

## Eviction configuration and stash bound

Size the batch optimization for this access schedule:

- do not evict the requested read path;
- perform two deterministic evictions per logical query;
- use two different deterministic paths (`evict_freq=2`, `evict_group=1`,
  `evict_on_read=false`).

In `tests/oram.cpp`, this is the `OriginalORAM` configuration. Its checked-in
log has 142,400,000,000 samples, a maximum observed load of 23, and the fitted
tail

```text
log2 Pr[L >= k] = -1.3311364218*k - 3.4137756687
R^2 = 0.9967359.
```

Representative empirical tail points are:

| load `k` | samples with `L >= k` | empirical `log2 Pr[L >= k]` |
| ---: | ---: | ---: |
| 8 | 7,732,338 | -14.17 |
| 12 | 242,102 | -19.17 |
| 16 | 6,358 | -24.42 |
| 20 | 187 | -29.50 |
| 22 | 17 | -32.96 |
| 23 | 2 | -36.05 |

Solving the fitted line for the per-query stash target gives:

| stash tail target | fitted crossing | first integer at or below target |
| ---: | ---: | ---: |
| `2^-58` | 41.01 | 42 |
| `2^-60` | 42.51 | 43 |
| `2^-64` | 45.51 | 46 |
| `2^-80` | 57.53 | 58 |

Use 46 as the provisional `2^-64` stash capacity for the calculations below.
If the two-slot engineering margin implicit in choosing 33 for the current
scheme is to be preserved, use 48 instead and rerun the hash analysis. The
physical Circuit ORAM stash is exercised throughout read/writeback processing,
so its capacity should still be selected from the per-query overflow goal. The
more relaxed batch-event budget applies specifically to the one hash-table
construction shared by the batch.

Before implementation, verify that the actual batch writeback performs exactly
this two-distinct-path schedule and add stash-load instrumentation at real
batch boundaries. The existing log is the correct experimental proxy for the
stated schedule, but it is not a trace of the full batched recursive-ORAM call
path.

The batch writeback path now propagates `evict_on_read` consistently. When it
is enabled, `BatchReadAndRemove` retains the de-duplicated requested paths and
`BatchWriteBack` uses each one for the first partial eviction, followed by the
configured deterministic evictions. When it is disabled, the insertion path
is still read and written but is not evicted; only the subsequent deterministic
paths are evicted. The in-memory and `DISK_IO` branches implement the same
schedule.

## Proposed table

Use compile-time parameters so the probability and performance choices remain
explicit:

| Symbol | Meaning | Initial value to evaluate |
| --- | --- | --- |
| `S` | Circuit ORAM stash capacity | provisional 46 for the batch scheme |
| `b` | main bucket size | 4 |
| `m` | number of main buckets | 46, 64, 74/80, 92, and 128 candidates |
| `g` | shared overflow bucket size | 4 |
| `r` | independently keyed candidates built in fixed time | start with 1 |
| `C` | slots in one candidate | `m*b + g` |

The main area occupies `[0, m*b)`. Bucket `j` occupies
`[j*b, (j+1)*b)`, and the global overflow bucket occupies
`[m*b, m*b+g)`. All unused slots contain dummy blocks.

The table stores actual `Block_` objects, not indexes into the old stash. This
keeps a query to exactly `b + g` block reads and avoids a second
data-dependent access through an index. The original flat stash remains an
untouched recovery copy until a candidate has been selected. Once online
probing begins, the selected table is the logical authoritative stash.

The table is a local workspace whose lifetime ends inside
`BatchReadAndRemove`. This avoids permanently allocating one padded table for
every recursive ORAM level. Peak EPC/heap usage still needs to be measured and
reported.

## Hash and randomness requirements

For a real block with UID `u`, compute

```text
tag = PRF(batch_key, STASH_MAIN_DOMAIN || canonical_encode(u))
bucket = reduce(tag, m)
```

The implementation should use AES-NI, with a separate local AES context. It
must not reuse the global external-memory encryption key or the long-lived
salts used by the oblivious map. Generate a fresh 128-bit or 256-bit key for
each candidate after the stash snapshot is fixed. Wipe candidate keys and tag
metadata after restoration.

For power-of-two `m`, reduction is a mask. For other `m`, use enough PRF bits
that reduction bias is included in the analysis and is far below the failure
budget. A power-of-two table simplifies both the implementation and the
probability model, but may cost more build work.

Freshness has two purposes:

- the table layout is independent of Circuit ORAM positions and stash load;
- bucket traces cannot be correlated across batches to recognize repeated
  UIDs.

Hash by UID, never by the block's assigned tree position or its current tree
location. Position prefixes can be highly skewed, for example when many
blocks are assigned to a crowded leaf. That skew invalidates the balls-into-
bins analysis and can make overflow adversarially correlated with the ORAM
state.

## Oblivious construction

Introduce a small helper, tentatively
`omap/odsl/one_time_stash_hash.hpp`, parameterized by `Block_`, `UidType`,
`S`, `b`, `m`, `g`, and `r`. A hash function injection point should be
available in tests so collision cases are deterministic.

### 1. Tag and oblivious histogram pass

Copy the `S` stash blocks into candidate records. Each record contains its
block and small metadata: `isReal`, `bucket`, `rankInBucket`, `isResidual`,
and `destination` or an equivalent sort key.

Maintain `m` 32-bit counters initialized to zero. For every one of the `S`
input slots:

1. Compute the AES tag and bucket for the slot. Dummy blocks receive dummy
   metadata and never increment a counter.
2. Scan every counter in a fixed pattern. Compare its public lane number with
   the secret bucket number and conditionally increment exactly the matching
   counter for a real block.
3. Obtain the old counter value through register operations and set
   `rankInBucket = old + 1`.
4. Mark a real block residual when `rankInBucket > b`.
5. Increment a register-resident residual count without branching.

The AVX-512 implementation can process sixteen 32-bit counters at a time with
lane IDs, equality masks, and masked adds. It must scan all chunks, including
a masked final chunk. Do not use a gather, scatter, or direct
`counters[bucket]` access. Add AVX2 and scalar fixed-scan fallbacks under the
same interface; the build already uses `-march=native`, but host and enclave
feature configurations still need to compile consistently.

For each non-residual block, define the unique destination

```text
bucket*b + (rankInBucket - 1).
```

For residuals, assign consecutive destinations starting at `m*b`. If the
residual rank exceeds `g`, set a secret candidate-failure bit and turn that
record into a dummy in the candidate. The original flat stash remains intact,
so a failed candidate never loses data.

The pass count and memory addresses must be identical for every stash load.
Do not choose `m`, stop early, or allocate a differently sized buffer based on
the live stash count.

### 2. Oblivious sort

Sort the `S` records by destination with dummies last, using the existing
bitonic sorting network in `omap/algorithm/bitonic.hpp`. Sorting by the
destination is equivalent to sorting main entries by bucket and rank and
placing residuals at the end. Use separate metadata and block payloads where
that reduces movement of large `Block_` values.

Only the `S` source records need sorting. Initialize the full `C`-slot output
to dummies and copy the sorted records into its prefix.

### 3. Reverse compaction into bucket slots

Build a `C+1` prefix array for occupied destination slots by scanning public
bucket and slot numbers:

- main slot `(j, k)` is occupied iff `counter[j] > k`, capped at `b`;
- overflow slot `k` is occupied iff `residualCount > k`, capped at `g`.

This is a sequential scan; no secret destination is used as an address. Feed
the prefix array and the compact sorted records to the existing
`Algorithm::OrDistributeSeparateMark` primitive in
`omap/algorithm/or_compact_shuffle.hpp`. The result has every real block at
its unique destination and dummies in all gaps.

Validate the exact prefix semantics with small exhaustive tests before
integrating the helper. The primitive is already used in `par_omap.hpp` and
`linear_oram.hpp` for the same compact-then-distribute pattern.

### 4. Fixed candidate selection

If `r > 1`, always construct all `r` candidates with independent keys,
regardless of whether the first succeeds. Select the first successful table
and its key into a canonical buffer with conditional moves. A secret pointer
or a branch choosing a candidate would reveal which candidate overflowed.

This fixed selection hides the build-time branch, but it does not make the
selected layout an unconditional random hash table. Whenever candidate 0
fails, the query trace comes from a later key conditioned on successful
placement. A hybrid argument bounds the resulting trace deviation by roughly
the probability that candidate 0 fails, not by the probability that all `r`
candidates fail. Multiple candidates are therefore an availability and
correctness fallback; they are not a substitute for sizing an individual
candidate to the privacy budget. A sharper distributional proof may improve
this bound, but must cover joint bucket traces of queried UIDs that are
members of the stash.

If all candidates fail, increment a dedicated performance counter and use the
legacy full-stash scan for that batch, or throw before changing positions or
tree buckets. The event is the explicitly budgeted failure/leakage event. A
variable number of rekeys is not acceptable because its timing reveals
information correlated with the stash load.

## Online lookup and duplicate requests

For every request, including misses and duplicates, read and compare all `b`
entries in one main bucket and all `g` entries in the global overflow bucket.
Use the existing block-level conditional removal logic. Never skip the global
bucket based on the secret residual count.

UIDs are sorted, so the first occurrence flag is

```text
first = (i == 0) || (uid[i] != uid[i-1]).
```

For a first occurrence, probe `PRF(batch_key, STASH_MAIN_DOMAIN || uid)`. For
every later duplicate, probe an independently pseudorandom dummy bucket, for
example

```text
PRF(batch_key, STASH_DUP_DOMAIN || request_index).
```

Select between the two bucket numbers without a branch. Compute both if
needed to keep instruction traces uniform. The global bucket is still probed.
After all path reads, the existing `duplicateVal` propagation supplies the
first result to duplicate requests.

This masking is necessary. If all occurrences probed `hash(uid)`, repeated
main-bucket addresses would reveal duplicate UIDs and defeat the purpose of
`deDuplicatePoses` for the stash portion of the access.

The intended leakage argument is computational rather than a claim that every
query touches the same address: for each distinct UID, the observed main
bucket is pseudorandom under a one-time key and independent of the UID and
Circuit ORAM state. The duplicate masking makes all later occurrences look
like fresh random probes as well.

## Restoring the Circuit ORAM stash

Before either branch of `BatchReadAndRemove` returns:

1. Obliviously compact all non-dummy blocks in the canonical `C`-slot table to
   the front with `Algorithm::OrCompactSeparateMark`.
2. Copy exactly the first `S` compacted slots over `path[0:S]`. Since the
   selected table began with at most `S` real blocks and lookups only remove
   blocks, this prefix contains every remaining real block followed by
   dummies.
3. Destroy or wipe the temporary table, tags, counters, and candidate keys.

The restore runs even when no request hit the stash. It must happen before
`BatchWriteBack`, `evictPath`, `WriteNewBlockToPath`, `GetStash`, or any other
code observes the stash again. A small RAII guard is advisable to keep the
flat-stash invariant on ordinary C++ exceptions; integrity-failure handling
may still terminate according to the existing policy.

## Correctness invariants

The implementation should assert these invariants in debug/test builds:

1. Before construction and after restoration, `path[0:S]` is the sole
   authoritative stash.
2. During online probing, the canonical table represents exactly the multiset
   that was in the flat stash, minus successful reads.
3. Every real block has exactly one copy in the canonical table.
4. Every main-table block is in the bucket derived from its UID and selected
   candidate key.
5. Every residual block is in the global overflow bucket.
6. At most `g` residuals exist for a successful candidate.
7. Restoration preserves every unqueried block and never produces more than
   `S` real blocks.
8. Table construction and restoration do not inspect or change
   `Block_::position`; subsequent Circuit ORAM eviction sees the same blocks
   and positions it would have seen under a flat stash scan.

## Failure probability

Let `L` be the number of real stash blocks at the beginning of the batch. For
a fixed load `L = l`, let `(X_1, ..., X_m)` be the occupancy vector produced by
fresh UID hashing. Define the number of residual blocks as

```text
R_l = sum_j max(0, X_j - b).
```

A candidate fails when `R_l > g`. Under the independent uniform-bucket model,
the exact conditional probability is

```text
p_hash(l) =
  sum over x_1+...+x_m=l and sum_j max(0,x_j-b)>g
  l! / (m^l * product_j x_j!).
```

Compute this with a dynamic program over `(processed buckets, assigned balls,
residual count capped at g+1)`. Do not use a union bound over individual
buckets when selecting final parameters; the exact DP is inexpensive for
`S <= 50` and captures multiple overfull buckets.

Circuit ORAM's real-block UIDs are unique. The conditional calculation may
therefore fix an arbitrary, even adversarial, set of `l` distinct UIDs: a
fresh PRF key supplies the bucket randomness. If raw AES is modeled as a PRP,
its outputs for distinct inputs are sampled without replacement from the
128-bit block space. Include the resulting PRP/PRF switching term (at most on
the order of `l^2 / 2^128`) in the final bound, along with any range-reduction
bias.

Because the hash key is fresh and independent of the Circuit ORAM state, one
candidate has total failure probability

```text
P_fail = sum_l Pr[L=l] * p_hash(l).
```

For `r` fixed, independently keyed candidates, all of which are always built,
the probability that no usable table exists is

```text
P_fail(r) = sum_l Pr[L=l] * p_hash(l)^r.
```

It is important not to compute `(sum_l Pr[L=l] * p_hash(l))^r`: the candidates
share the same stash load.

If the Circuit ORAM analysis supplies upper tail bounds
`Q_l >= Pr[L >= l]` instead of exact point probabilities, use monotonicity and
summation by parts. For `f(l) = p_hash(l)^r` and `f(0) = 0`,

```text
P_fail(r) <= sum_{l=1}^S (f(l) - f(l-1)) * Q_l.
```

This permits a conservative convolution without pretending that an empirical
histogram is a proof.

### Sanity check from the checked-in trace

Use `logs/circuit_oram_stash_original_no_read_two_paths.log`, which matches no
read-path eviction and two deterministic evictions on different paths. The
following table fixes `S=46`, `b=4`, `g=4`, and `r=1`.

Three estimates are shown:

- **empirical** convolves the observed point counts and assigns zero mass
  above the maximum observed load 23;
- **fitted tail** uses the exponential tail regression given above through
  load 46;
- **upper envelope** raises the fitted intercept by 0.5611 bits, enough to put
  the line above every observed empirical tail point. This is a useful
  sensitivity calculation, not a statistical or cryptographic confidence
  bound.

| main buckets `m` | nominal load factor `46/(4m)` | empirical `log2 P_fail` | fitted-tail `log2 P_fail` | upper-envelope `log2 P_fail` | minimum batch `q` allowed by envelope |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 46 | 0.2500 | -52.93 | -53.18 | -52.62 | 2,662 |
| 64 | 0.1797 | -56.69 | -56.95 | -56.39 | 196 |
| 72 | 0.1597 | -58.04 | -58.29 | -57.73 | 78 |
| 74 | 0.1554 | -58.35 | -58.61 | -58.04 | 63 |
| 80 | 0.1438 | -59.24 | -59.50 | -58.94 | 34 |
| 92 | 0.1250 | -60.84 | -61.10 | -60.54 | 12 |
| 128 | 0.0898 | -64.63 | -64.88 | -64.32 | 1 |

The final column is

```text
ceil(P_fail / 2^-64) = ceil(2^(64 + log2(P_fail))).
```

It is the smallest batch size for which the one-build-per-batch failure meets
the `2^-64` per-query accounting, before dividing the budget among multiple
hash-table builds or other failure sources.

For a minimum batch size of 64, the raw empirical calculation says 72 buckets
is just sufficient (`2^-58.04`). The upper-envelope sensitivity calculation
needs 74 buckets (`2^-58.04`). Starting with 80 buckets gives approximately
0.94 bits of margin against a `2^-58` batch target while preserving the same
eight-entry online lookup. If the implementation only enables hashing at a
larger public batch threshold, fewer buckets may suffice; if it enables it for
smaller batches or splits the budget across many recursive ORAM instances,
use more.

If the physical stash is rounded up from 46 to 48 for engineering margin, 80
main buckets have nominal load factor 0.15. The empirical and upper-envelope
convolutions remain `-59.24` and `-58.94` respectively to the shown precision,
because the probability mass at loads 47 and 48 is already extremely small.

The trace was produced by the single-access `CircuitORAM.StashLoad` benchmark
at one ORAM size, not by every real batched recursive-ORAM workload. Its
regression tool, `tools/tail_bound_analysis.py`, performs a curve fit and does
not provide a cryptographic upper confidence bound.

Add `tools/stash_hash_failure.py` to:

- compute `p_hash(l)` exactly for every `0 <= l <= S`;
- parse a stash-load histogram and perform the convolution;
- support fixed candidate count `r`;
- report both all-candidates-fail probability and the selected-layout hybrid
  bound;
- consume proved or conservative tail bounds when available;
- report per-load contributions, sensitivity to unseen tail mass, and the
  requested overall failure budget;
- brute-force small cases to validate the DP.

The stash-load input must be measured at the beginning of real batch reads for
every tree-backed recursive level and relevant ORAM size/configuration.
Single-access samples are useful but insufficient. Empirical confidence
intervals and tail-model extrapolation should be labeled as such. The final
security argument should use a stash tail bound that holds for every allowed
logical access sequence; a random-access workload histogram cannot establish
that property.

### Failure-budget accounting

Treat hash-table failure as a batch event. If the contract is at most
`2^-64` failure per logical request and a batch has `q` requests, the entire
batch may spend at most `q * 2^-64`, divided among all failure sources.

The repository currently mentions Circuit ORAM overflow, cuckoo insertion
stash overflow, and parallel load-balancing failure. The new event must be
added to that accounting. Also union-bound across every one-time stash table
built by all recursive levels, cuckoo tables, and shards during one logical
operation. A component-level estimate just below `2^-64` is therefore not
automatically sufficient. Select final `(m, r)` from the composed budget with
margin; do not select it solely from the empirical table above. Apply the
privacy budget to an individual candidate/selection hybrid as well as applying
the correctness/availability budget to the all-candidates-fail event.

## Why tree-top entries are out of scope

Do not directly distribute cached top-tree entries to lower subtrees according
to their assigned positions as part of this change. Position prefixes are not
an independent load-balancing hash: a crowded leaf or path can send many
entries to the same subtree. This can overflow a target even if average load
looks small, and it couples the failure event to the ORAM state.

The one-time stash table works because its bucket assignment is a fresh PRF of
the UID. If top-tree accesses are optimized later, use the same independent
UID-based snapshot/index/restore principle and give it a separate failure and
correctness analysis. Do not move top-tree blocks into subtrees or alter their
Circuit ORAM placement invariant in this implementation.

## Code changes

### New helper

Add `omap/odsl/one_time_stash_hash.hpp` with an interface along these lines:

```cpp
template <typename Block, typename Uid, size_t S, size_t B,
          size_t M, size_t G, size_t Candidates>
class OneTimeStashHash {
 public:
  bool Build(const Block* flatStash);
  bool ReadAndRemove(const Uid& uid, uint64_t requestIndex,
                     bool firstOccurrence, typename BlockData<Block>& out);
  void Restore(Block* flatStash);
};
```

The concrete interface may use iterators and infer the data type from
`Block_`. Keep construction, candidate selection, lookup, and restoration in
one helper so the authoritative-storage transition is hard to misuse.

### Circuit ORAM integration

In `omap/odsl/circuit_oram.hpp`:

- add compile-time configuration at the end of the ORAM template parameter
  list, or a nested policy type, without disturbing existing positional
  template arguments;
- build the helper once at the start of `BatchReadAndRemove` when the public
  batch size is above a benchmarked threshold;
- replace both full-stash scans, in-memory and `DISK_IO`, with the same helper
  lookup;
- restore the flat stash before `duplicateVal` and before either branch
  returns;
- retain the legacy scan for small public batches and the all-candidates-fail
  path;
- leave single accesses and writeback/eviction untouched.

Avoid duplicating the lookup logic in the two storage branches by factoring a
small local callable or method.

### AES utility

Add a batch-tag helper in `omap/common/encutils.hpp/.cpp` only if it is useful
beyond the stash helper. It should use the BearSSL AES-NI implementation
already linked for host and enclave builds, create a separate context, support
domain-separated 128-bit inputs, and have a fixed-time fallback. Benchmark a
four- or eight-block AES path rather than invoking a heavyweight hash API for
every UID.

### Counters and diagnostics

Add release performance counters for at least:

- hash-table builds;
- candidate overflows;
- all-candidates-fail fallbacks;
- main and overflow entries probed;
- restore operations.

Never log the live stash load, selected candidate, bucket index, or residual
count in production. Test-only histograms may record stash load before build.

## Performance model and threshold

The legacy stash work is approximately

```text
q * S block comparisons.
```

The new work is approximately

```text
r * (S AES tags + S*m/vector_width histogram lanes
     + oblivious sort of S records + reverse distribution over C slots)
+ q * (one AES tag + b+g block comparisons)
+ oblivious compaction over C slots + S-block restore copy.
```

The build and restore are fixed costs amortized by the batch. A small-batch
crossover is expected, especially for large `Block_` payloads. Branching on
the already-observable/public batch size is acceptable in the current API; if
batch size is intended to be secret, it must be padded before this layer.

Benchmark the full operation, not only the eight-entry query loop. Compare:

- main bucket counts 46, 64, 74/80, 92, and 128;
- one lower-load-factor candidate versus fixed multi-candidate fallback;
- scalar, AVX2, and AVX-512 histograms;
- representative `Block_` sizes for internal and leaf recursive ORAM levels;
- fully cached and `DISK_IO` trees;
- batch sizes from 1 through the expected production range.

Use the measurements to set a conservative public
`stashHashMinBatchSize`. Allocation should be separated in profiles; if it is
material, introduce a reusable per-thread workspace only after ensuring its
memory is counted and cannot be used concurrently.

## Test plan

### Standalone helper tests

- Exhaust every load `0..S` with randomized block order and many keys. On
  successful build, every real UID must be found exactly once and every miss
  must remain a miss.
- Inject hashes producing exactly 4, 5, 8, and 9 blocks in one main bucket to
  test main capacity, global residual placement, and overflow detection.
- Test multiple overfull main buckets whose residual counts sum to `g` and
  `g+1`.
- Test arbitrary dummy placement and ensure dummy UIDs do not affect counters.
- Remove an arbitrary subset, restore, and compare the resulting block
  multiset, including positions and data, with a reference flat scan.
- Force candidate 0 to fail and candidate 1 to succeed; verify canonical
  selection without a secret-dependent pointer. Force all candidates to fail
  and verify the legacy fallback preserves all data.
- Verify the duplicate-probe domain produces a dummy lookup for later sorted
  duplicates while `duplicateVal` preserves the public API result.
- Instantiate 32- and 64-bit UID/position types and representative payload
  sizes.

### Circuit ORAM integration tests

- Differentially run legacy and hash-table batch modes with the same logical
  operations and compare outputs and the multiset of all stash/tree blocks.
- Cover hits in a main bucket, hits in the global bucket, tree-only hits,
  misses, duplicate UIDs, removals, and writeback flags.
- Run repeated `BatchReadAndRemove`/`BatchWriteBack` cycles and validate against
  a normal map.
- Exercise fully cached and `DISK_IO` accessor paths with freshness checking
  enabled.
- Test batch size zero and the public crossover boundary.

### Obliviousness/trace tests

Add a test instrumentation mode that records abstract memory-region accesses
and pass counts. For equal public parameters, construction and restoration
must have the same trace shape for every stash load. Every online request must
record one complete main bucket and the complete global bucket. Distinct and
duplicate request sequences should have indistinguishable bucket-index
distributions under deterministic test keys/domains; later duplicates must
not deterministically repeat the real UID bucket.

### Probability and long-run tests

- Cross-check the conditional DP against exhaustive enumeration for small
  `m`, `b`, `g`, and `l`.
- Cross-check it against Monte Carlo simulation at probabilities large enough
  to measure.
- Extend stash-load collection to record load at real batch-build boundaries,
  separated by ORAM level and configuration.
- Re-run failure composition whenever eviction parameters, stash capacity,
  bucket parameters, or recursive layout change.

## Rollout

1. Implement the exact failure-analysis tool, including selected-layout
   conditioning, and collect batch-boundary stash distributions. Choose a
   provisional parameter policy; do not enable the optimization by default
   yet.
2. Implement the standalone helper with scalar fixed-scan histogram,
   deterministic hash injection, fixed candidate construction, and restore.
   Land exhaustive correctness tests.
3. Integrate behind a compile-time policy or feature flag in both
   `BatchReadAndRemove` branches. Add differential and long-running tests.
4. Add AES-NI batching and AVX-512/AVX2 histogram paths, retaining scalar
   fallbacks. Establish the batch-size crossover and memory overhead.
5. Complete composed failure-budget review across recursive levels, tables,
   and shards. Enable by default only after the bound and trace tests meet the
   target with margin.

## Acceptance criteria

The optimization is ready to enable when all of the following hold:

- successful batches read exactly `b+g = 8` stash blocks per request;
- the selected parameter policy meets the composed failure budget using a
  conservative stash-load bound, not only regression or observed maxima;
- an individual candidate meets the privacy hybrid budget, fixed candidate
  work avoids a build-time retry leak, and duplicate dummy probes hide request
  equality;
- all surviving blocks are restored to the ordinary `S`-slot stash before
  writeback;
- legacy and optimized modes pass differential tests in cached, external
  memory, and enclave builds;
- end-to-end batch throughput improves beyond a documented public crossover,
  including construction and restoration costs;
- no position-based load balancing of stash or top-tree entries is introduced.
