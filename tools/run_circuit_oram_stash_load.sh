#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

test_binary="${CIRCUIT_ORAM_TEST_BINARY:-$repo_root/build/tests/test_oram}"
warmup_windows="${STASH_LOAD_WARMUP_WINDOWS:-100000}"
window_count="${STASH_LOAD_WINDOWS:-4450000000}"

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
