#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

default_test_binary="$repo_root/build/tests/test_oram"
test_binary="${CIRCUIT_ORAM_TEST_BINARY:-$default_test_binary}"
warmup_windows="${STASH_LOAD_WARMUP_WINDOWS:-100000}"
window_count="${STASH_LOAD_WINDOWS:-4450000000}"

# The stash-load variants are selected in tests/oram.cpp and most of the ORAM
# implementation is header-only.  An existing executable can therefore be
# stale even when the runner itself has not changed.  Keep custom binaries
# caller-managed, but always bring the repository's default target up to date.
if [[ "$test_binary" == "$default_test_binary" ]]; then
  echo "Building Circuit ORAM test binary"
  cmake --build "$repo_root/build" --target test_oram --parallel
fi

if [[ ! -x "$test_binary" ]]; then
  echo "Test binary not found or not executable: $test_binary" >&2
  exit 1
fi

mkdir -p logs

declare -a variants=(
  "current:circuit_oram_stash_current_partial_read_two_same.log"
  "batched_current:circuit_oram_stash_batched_current.log"
  "original:circuit_oram_stash_original_no_read_two_paths.log"
  "partial_two_paths:circuit_oram_stash_partial_read_two_paths.log"
)

for variant_and_log in "${variants[@]}"; do
  variant="${variant_and_log%%:*}"
  log_file="${variant_and_log#*:}"
  echo "Starting $variant; output: logs/$log_file"
  CIRCUIT_ORAM_STASH_LOAD_VARIANT="$variant" \
    STASH_LOAD_WARMUP_WINDOWS="$warmup_windows" \
    STASH_LOAD_WINDOWS="$window_count" \
    "$test_binary" --gtest_filter=CircuitORAM.StashLoad \
    > "logs/$log_file" 2>&1
  echo "Completed $variant"
done

echo "All CircuitORAM.StashLoad variants completed."
