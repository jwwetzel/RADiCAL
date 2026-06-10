#!/usr/bin/env bash
# ============================================================================
# pull_argon_data.sh — parallel mirror of the RADiCAL data tree from Argon.
#
#   Source : argon.hpc.uiowa.edu:/Shared/lss_yonel/jwwetzel/   (ssh port 40, Duo 2FA)
#   Dest   : $DEST (default ~/RADiCAL_DATA), tree structure preserved
#
# 2FA strategy: SSH connection multiplexing (ControlMaster). You authenticate
# once per master connection (MASTERS Duo prompts total, sequential, at the
# start); every parallel rsync then rides an existing authenticated connection
# with no further prompts. sshd typically caps ~10 sessions per connection,
# so total parallelism is spread across the masters (PARALLEL/MASTERS each).
#
# Usage:
#   tools/pull_argon_data.sh                 # 24-way, 3 Duo prompts
#   PARALLEL=8 MASTERS=1 tools/pull_argon_data.sh   # gentler: 1 Duo prompt
#   DEST=/Volumes/Big/RADiCAL tools/pull_argon_data.sh
#
# Re-running is safe and cheap: rsync skips files already present and
# identical (size+mtime), so an interrupted pull just resumes.
# Run inside caffeinate so the Mac doesn't sleep mid-transfer:
#   caffeinate -i tools/pull_argon_data.sh
# ============================================================================
set -uo pipefail

# ----- config (override via environment) -----
ARGON_USER="${ARGON_USER:-jwwetzel}"
ARGON_HOST="${ARGON_HOST:-argon.hpc.uiowa.edu}"
ARGON_PORT="${ARGON_PORT:-40}"
SRC="${SRC:-/Shared/lss_yonel/jwwetzel}"
DEST="${DEST:-$HOME/RADiCAL_DATA}"
PARALLEL="${PARALLEL:-24}"   # total concurrent rsync streams
MASTERS="${MASTERS:-3}"      # ssh master connections = number of Duo prompts

PER=$(( (PARALLEL + MASTERS - 1) / MASTERS ))
if (( PER > 8 )); then
  echo "WARNING: $PER sessions per master exceeds the usual sshd MaxSessions"
  echo "         headroom (10). If transfers fail to start, raise MASTERS or"
  echo "         lower PARALLEL."
fi

CTL_DIR="$HOME/.ssh/argon-ctl"; mkdir -p "$CTL_DIR" "$DEST"
WORK="$(mktemp -d /tmp/argonpull.XXXXXX)"
FAIL="$WORK/failures.txt"; : > "$FAIL"

sock() { echo "$CTL_DIR/m$1.sock"; }

# ----- 1. open the master connections (one Duo prompt each, sequential) -----
for (( i=0; i<MASTERS; i++ )); do
  S="$(sock $i)"
  if ssh -o ControlPath="$S" -O check "$ARGON_USER@$ARGON_HOST" 2>/dev/null; then
    echo "master $i: already connected (no new 2FA needed)"
  else
    echo "master $i: connecting — expect a Duo prompt..."
    ssh -p "$ARGON_PORT" -fN \
        -o ControlMaster=yes -o ControlPath="$S" -o ControlPersist=600 \
        -o ServerAliveInterval=30 -o Compression=no \
        "$ARGON_USER@$ARGON_HOST" || { echo "master $i FAILED"; exit 1; }
  fi
done

mssh() { ssh -o ControlPath="$(sock 0)" "$ARGON_USER@$ARGON_HOST" "$@"; }

# ----- 2. manifest: every file under SRC, NUL-delimited (safe for spaces) -----
echo "building file manifest from $SRC ..."
mssh "cd '$SRC' && find . -type f -print0" > "$WORK/manifest0" \
  || { echo "manifest failed"; exit 1; }
TOTAL=$(tr -cd '\0' < "$WORK/manifest0" | wc -c | tr -d ' ')
echo "  $TOTAL files; remote size: $(mssh "du -sh '$SRC' 2>/dev/null | cut -f1")"
(( TOTAL > 0 )) || { echo "nothing to copy"; exit 0; }

# ----- 3. split the manifest round-robin across the masters -----
python3 - "$WORK/manifest0" "$WORK/chunk" "$MASTERS" << 'EOF'
import sys
src, pre, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
items = open(src,'rb').read().split(b'\0')[:-1]
outs = [open(f"{pre}{i}", 'wb') for i in range(n)]
for k, it in enumerate(items): outs[k % n].write(it + b'\0')
for o in outs: o.close()
EOF

# ----- 4. parallel pull: PER streams per master, each file with 3 retries -----
echo "pulling with $PARALLEL streams across $MASTERS connections -> $DEST"
export ARGON_USER ARGON_HOST SRC DEST FAIL
transfer_one() {  # $1 = control socket, $2 = ./relative/path
  local rel="${2#./}"
  for try in 1 2 3; do
    rsync -aR --partial -e "ssh -o ControlPath=$1 -o Compression=no" \
      "$ARGON_USER@$ARGON_HOST:$SRC/./$rel" "$DEST/" 2>/dev/null \
      && { echo "done: $rel"; return 0; }
    sleep $(( try * 5 ))
  done
  echo "FAILED: $rel" | tee -a "$FAIL" >&2
}
export -f transfer_one

for (( i=0; i<MASTERS; i++ )); do
  xargs -0 -n1 -P "$PER" -I{} bash -c 'transfer_one "$1" "$2"' _ "$(sock $i)" {} \
    < "$WORK/chunk$i" &
done
wait

# ----- 5. summary -----
NFAIL=$(wc -l < "$FAIL" | tr -d ' ')
echo "----------------------------------------------------------------------"
echo "complete: $(( TOTAL - NFAIL )) / $TOTAL files; local size: $(du -sh "$DEST" | cut -f1)"
if (( NFAIL > 0 )); then
  echo "$NFAIL failures (listed in $FAIL) — just re-run the script; it resumes."
else
  rm -rf "$WORK"
fi
echo "(masters stay alive ~10 min — an immediate re-run needs no new Duo prompt)"
