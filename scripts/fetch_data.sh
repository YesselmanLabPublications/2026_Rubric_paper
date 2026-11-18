#!/usr/bin/env bash
set -euo pipefail

BASE_URL="${RUBRIC_DATA_URL:-https://zenodo.org/records/24214090/files}"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DOWNLOADS="$REPO_ROOT/.downloads"

fetch() {
    local file="rubric_$1.zip"
    local archive="$DOWNLOADS/$file"
    local url="$BASE_URL/$file?download=1"
    if [[ -f "$archive" ]] && unzip -tq "$archive" >/dev/null 2>&1; then
        echo "$file is already downloaded."
    else
        echo "Downloading $file from Zenodo..."
        mkdir -p "$DOWNLOADS"
        if ! curl -L --fail --retry 5 --retry-delay 10 -C - -o "$archive" "$url"; then
            echo "Could not download $file. Run the script again to resume the download." >&2
            exit 1
        fi
    fi
    if ! unzip -tq "$archive" >/dev/null 2>&1; then
        echo "$file is not a valid archive (got $(wc -c <"$archive") bytes)." >&2
        rm -f "$archive"
        exit 1
    fi
    echo "Adding the files of $file that are not yet in $REPO_ROOT/data ..."
    unzip -qn "$archive" -d "$REPO_ROOT"
    rm -f "$archive"
    rmdir "$DOWNLOADS" 2>/dev/null || true
}

case "${1:-data}" in
    data | structures | program_runs) fetch "${1:-data}" ;;
    all) for name in data structures program_runs; do fetch "$name"; done ;;
    *)
        echo "Usage: bash scripts/fetch_data.sh [data | structures | program_runs | all]" >&2
        exit 1
        ;;
esac
echo "Done."
