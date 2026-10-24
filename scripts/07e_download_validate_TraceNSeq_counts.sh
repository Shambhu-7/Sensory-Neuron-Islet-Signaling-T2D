#!/bin/bash
set -euo pipefail

OUT="data/raw/E-MTAB-12940"

mkdir -p "$OUT"
mkdir -p reproducibility

BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/940/E-MTAB-12940/Files"

echo "===================================================="
echo "SCRIPT 07e: DOWNLOAD TRACE-N-SEQ COUNT MATRICES"
echo "===================================================="

download_file () {

    NAME="$1"
    URL="$BASE/$NAME"
    DEST="$OUT/$NAME"
    PART="$DEST.part"

    echo
    echo "----------------------------------------------------"
    echo "FILE: $NAME"
    echo "----------------------------------------------------"
    echo "URL: $URL"

    if [ -s "$DEST" ]; then

        echo "Existing completed file found."

    else

        echo "Downloading..."

        curl \
          --fail \
          --location \
          --retry 10 \
          --retry-delay 5 \
          --retry-all-errors \
          --continue-at - \
          --output "$PART" \
          "$URL"

        mv "$PART" "$DEST"

    fi

    echo
    echo "Local file:"
    ls -lh "$DEST"

    echo "Bytes:"
    wc -c "$DEST"

    echo "MD5:"
    md5 "$DEST"

    echo
}

download_file "smartseq_raw_counts.txt"
download_file "smartseq_raw_counts2.txt"

echo
echo "===================================================="
echo "BOTH COUNT MATRICES DOWNLOADED"
echo "===================================================="

ls -lh "$OUT"/smartseq_raw_counts*.txt

echo
echo "SCRIPT 07e COMPLETED SUCCESSFULLY"
