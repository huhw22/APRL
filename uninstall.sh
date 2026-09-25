#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./uninstall.sh [options]

Remove only files listed in the installation receipt and, by default, the
disposable compilation directory. Project source files are never removed.

Options:
  --prefix DIR          installation prefix (default: /usr/local)
  --manifest FILE       explicit installation receipt
  --build-dir DIR       compilation directory to clear (default: build-install)
  --keep-build          retain compilation products
  --dry-run             print removals without changing files
  -h, --help            show this help
EOF
}

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
prefix="/usr/local"
manifest=""
build_dir="$source_root/build-install"
keep_build=false
dry_run=false

while (($#)); do
  case "$1" in
    --prefix)
      [[ $# -ge 2 ]] || { echo "--prefix needs a value" >&2; exit 2; }
      prefix="$2"; shift 2 ;;
    --manifest)
      [[ $# -ge 2 ]] || { echo "--manifest needs a value" >&2; exit 2; }
      manifest="$2"; shift 2 ;;
    --build-dir)
      [[ $# -ge 2 ]] || { echo "--build-dir needs a value" >&2; exit 2; }
      build_dir="$2"; shift 2 ;;
    --keep-build) keep_build=true; shift ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

absolute_path() {
  local value="$1"
  local candidate=""
  if [[ "$value" == /* ]]; then
    candidate="$value"
  else
    candidate="$(pwd -P)/$value"
  fi
  realpath -m -- "$candidate"
}

prefix="$(absolute_path "$prefix")"
build_dir="$(absolute_path "$build_dir")"
[[ -n "$manifest" ]] ||
  manifest="$prefix/share/unnamed_fel_program/install-manifest.txt"
manifest="$(absolute_path "$manifest")"

if [[ -z "$prefix" || "$prefix" == "/" ]]; then
  echo "Refusing an empty or root install prefix" >&2
  exit 2
fi
if [[ ! -f "$manifest" ]]; then
  echo "Installation receipt not found: $manifest" >&2
  echo "No installed file was removed. Use --manifest only with a receipt created by install.sh." >&2
  exit 1
fi

manifest_copy="$(mktemp "${TMPDIR:-/tmp}/fel-uninstall-manifest.XXXXXX")"
cp -- "$manifest" "$manifest_copy"
trap 'rm -f -- "$manifest_copy"' EXIT

removed=0
while IFS= read -r installed_path || [[ -n "$installed_path" ]]; do
  [[ -n "$installed_path" ]] || continue
  normalized_path="$(realpath -m -- "$installed_path")"
  case "$normalized_path" in
    "$prefix"/*) ;;
    *)
      echo "Unsafe receipt entry outside prefix; refusing uninstall: $installed_path" >&2
      exit 1 ;;
  esac
  if [[ -f "$normalized_path" || -L "$normalized_path" ]]; then
    if $dry_run; then
      echo "would remove: $normalized_path"
    else
      rm -f -- "$normalized_path"
      echo "removed: $normalized_path"
    fi
    removed=$((removed + 1))
  fi
done < "$manifest_copy"

data_dir="$prefix/share/unnamed_fel_program"
if ! $dry_run && [[ -d "$data_dir" ]]; then
  find "$data_dir" -depth -type d -empty -exec rmdir -- {} \; 2>/dev/null || true
fi

if ! $keep_build; then
  if [[ -z "$build_dir" || "$build_dir" == "/" ||
        "$build_dir" == "$source_root" ]]; then
    echo "Refusing unsafe build directory cleanup: $build_dir" >&2
    exit 2
  fi
  case "$source_root/" in
    "$build_dir/"*)
      echo "Refusing a build directory that contains the source tree: $build_dir" >&2
      exit 2 ;;
  esac
  if [[ -d "$build_dir" ]]; then
    if $dry_run; then
      echo "would remove compilation directory: $build_dir"
    else
      cmake -E remove_directory "$build_dir"
      echo "removed compilation directory: $build_dir"
    fi
  fi
fi

if $dry_run; then
  echo "Dry run complete: $removed installed files would be removed."
else
  echo "Uninstall complete: $removed installed files removed."
fi
echo "Source tree preserved: $source_root"
