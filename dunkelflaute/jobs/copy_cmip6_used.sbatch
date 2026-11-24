#!/bin/bash
# Copy every CMIP6 file the paper's extraction scripts read (results/cmip6/used_files/ALL.txt, about 18,400 files,
# about 6 TB, built by scripts/cmip6/used_files.py) from the DKRZ pool to /work/ab0246/cmip6/data, before the
# pool cleanup of 2026-11 (CMIP7 storage). Paths keep the pool layout below data/. The pool's data/ is a symlink to
# /work/ik1017/CMIP6/data/CMIP6, which rsync cannot traverse, so the copy reads from the resolved directory.
# Files already present (same size and time) are skipped, so the job can be resubmitted.
#SBATCH --job-name=copy_cmip6_used
#SBATCH --partition=shared
#SBATCH --account=ab0246
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/copy_cmip6_used_%j.log
set -uo pipefail
R=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute
SRC=/work/ik1017/CMIP6/data/CMIP6/
DST=/work/ab0246/cmip6/data/
LIST=$R/results/cmip6/used_files/ALL.txt
TMP=$(mktemp -d)
mkdir -p "$DST"
echo "=== $(date) start, $(wc -l < "$LIST") files"
# 16 parallel streams over chunks of the list
split -n l/16 "$LIST" "$TMP/chunk_"
ls "$TMP"/chunk_* | xargs -P 16 -I{} rsync -a --files-from={} "$SRC" "$DST"
rc=$?
echo "=== $(date) rsync done, exit $rc"
n=$(find "$DST" -type f | wc -l)
echo "files at destination: $n"
# verify: every listed file exists with the source size
python3 - "$LIST" "$SRC" "$DST" <<'EOF'
import os, sys
L, S, D = sys.argv[1:]
bad = [p.strip() for p in open(L) if p.strip() and (not os.path.exists(D + p.strip())
       or os.path.getsize(D + p.strip()) != os.path.getsize(S + p.strip()))]
print("missing or size mismatch:", len(bad))
for p in bad[:20]:
    print("  ", p)
EOF
rm -rf "$TMP"
echo "=== $(date) DONE"
