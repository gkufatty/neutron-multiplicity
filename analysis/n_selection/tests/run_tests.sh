#!/usr/bin/env bash
set -euo pipefail
selection_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/n-selection-tests.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT
"${CXX:-g++}" -std=c++17 -O2 -Wall -Wextra -Werror -pedantic \
    -I"$selection_dir/headers" \
    "$selection_dir/sources/selection_input.cc" \
    "$selection_dir/tests/test_selection_input.cc" \
    -o "$test_dir/test_selection_input"
"$test_dir/test_selection_input" "$test_dir" "$@"
