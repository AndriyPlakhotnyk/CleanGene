#!/usr/bin/env bash
set -Eeuo pipefail

DB_DIR=${1:?Usage: build_kraken2_database.sh DB_DIR THREADS CLEAN_INTERMEDIATES DATABASE_SIZE}
THREADS=${2:-4}
CLEAN_INTERMEDIATES=${3:-true}
DATABASE_SIZE=${4:-${KRAKEN2_DATABASE_SIZE:-standard-8}}

case "$DATABASE_SIZE" in
  standard-8|standard_8|8)
    URL=${KRAKEN2_STANDARD_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz}
    MD5_URL=${KRAKEN2_STANDARD_MD5_URL:-https://genome-idx.s3.amazonaws.com/kraken/standard_08_GB_20260626/standard_08_GB.md5}
    ;;
  standard-16|standard_16|16)
    URL=${KRAKEN2_STANDARD_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260626.tar.gz}
    MD5_URL=${KRAKEN2_STANDARD_MD5_URL:-https://genome-idx.s3.amazonaws.com/kraken/standard_16_GB_20260626/standard_16_GB.md5}
    ;;
  standard|full|full-standard)
    URL=${KRAKEN2_STANDARD_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz}
    MD5_URL=${KRAKEN2_STANDARD_MD5_URL:-https://genome-idx.s3.amazonaws.com/kraken/standard_20260626/standard.md5}
    ;;
  *) echo "KRAKEN2_DATABASE_SIZE must be standard-8, standard-16, or standard" >&2; exit 2 ;;
esac

mkdir -p "$DB_DIR"
ARCHIVE="$DB_DIR/$(basename "$URL")"
CHECKSUMS="$DB_DIR/standard.md5"

echo "Downloading Kraken2 database: $URL"
if command -v aria2c >/dev/null 2>&1; then
  aria2c --continue=true --max-connection-per-server=8 --split=8 --min-split-size=4M \
    --max-tries=8 --retry-wait=10 --dir "$DB_DIR" --out "$(basename "$ARCHIVE")" "$URL"
else
  curl --fail --location --continue-at - --retry 8 --retry-all-errors --retry-delay 10 \
    --connect-timeout 30 --output "$ARCHIVE" "$URL"
fi
curl --fail --silent --show-error --location --retry 5 --output "$CHECKSUMS" "$MD5_URL"
(
  cd "$DB_DIR"
  expected=$(awk -v archive="$(basename "$ARCHIVE")" '$2 == archive {print; exit}' "$CHECKSUMS")
  [[ -n "$expected" ]]
  printf '%s\n' "$expected" | md5sum --check -
)
tar -xzf "$ARCHIVE" -C "$DB_DIR"
for file in hash.k2d opts.k2d taxo.k2d; do
  test -s "$DB_DIR/$file"
done
if [[ "$CLEAN_INTERMEDIATES" =~ ^(1|true|yes|on)$ ]]; then
  rm -f "$ARCHIVE"
fi
printf 'install_date\t%s\ndatabase_size\t%s\nsource_url\t%s\nthreads\t%s\n' \
  "$(date -Iseconds)" "$DATABASE_SIZE" "$URL" "$THREADS" \
  > "$DB_DIR/cleangene_database_metadata.tsv"
