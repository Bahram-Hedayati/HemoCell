#!/usr/bin/env bash
set -euo pipefail
example_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$example_dir/../.." && pwd)"
cmake -S "$repo_dir" -B "$repo_dir/build" -DBUILD_TESTING=OFF
cmake --build "$repo_dir/build" --target vessel --parallel "${JOBS:-4}"
