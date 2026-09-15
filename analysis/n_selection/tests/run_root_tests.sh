#!/usr/bin/env bash
# Requires the ROOT and duneanaobj environment used by macro.sh.
set -euo pipefail
selection_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/n-selection-root-tests.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT
"${CXX:-g++}" -std=c++17 -O2 -DN_SELECTION_NO_MAIN \
    "$selection_dir/sources/main.cc" "$selection_dir/sources/utils.cc" \
    "$selection_dir/sources/cuts.cc" "$selection_dir/sources/selection_input.cc" \
    "$selection_dir/tests/test_caf_selection.cc" -o "$test_dir/test_caf_selection" \
    -I"$selection_dir/headers" $(root-config --cflags --glibs) -lEG \
    -I"${DUNEANAOBJ_INC:?Set up duneanaobj first}" \
    -L"${DUNEANAOBJ_LIB:?Set up duneanaobj first}" \
    -lduneanaobj_StandardRecord -lduneanaobj_StandardRecord_dict
"$test_dir/test_caf_selection" "$test_dir"
