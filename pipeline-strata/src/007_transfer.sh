#!/bin/bash
#SBATCH --job-name=007_transfer
#SBATCH --output=src/logs/007_transfer.log
#SBATCH --error=src/logs/007_transfer.err
#SBATCH --time=4:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# Script: src/007_transfer.sh
# Purpose: Upload a staged directory to collaborators via SFTP (uses Python/paramiko).
#
# Uploads RECURSIVELY, recreating the staged tree on the remote:
#   <transfer>/<seqid>/<group>.tsv.gz  ->  <remote>/<name>/<seqid>/<group>.tsv.gz
#
# Usage:
#   bash src/007_transfer.sh --dry-run            # show the plan, connect to nothing
#   sbatch src/007_transfer.sh                    # upload
#   bash src/007_transfer.sh --list               # list what is already on the remote
#
#   --dir <path>   local staged directory   (default: GWAS_STAGING_DIR)
#   --name <name>  remote folder to upload into, created under the remote dir
#                  (default: GWAS_TRANSFER_SUBDIR, else the local directory name)
#
# Read from .env:
#   GWAS_STAGING_DIR         Local directory to upload (default for --dir)
# Read from the untracked server.env:
#   GWAS_TRANSFER_SUBDIR     Remote folder name        (default for --name)
#   GWAS_TRANSFER_HOST       SFTP hostname
#   GWAS_TRANSFER_USER       SFTP username
#   GWAS_TRANSFER_PASS       SFTP password
#   GWAS_TRANSFER_REMOTE_DIR Remote destination path

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${GWAS_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}}"

# Shared .env loader (src/lib/env.sh): an exported variable always beats .env.
_gwas_find_lib() {
  local c
  for c in "${SLURM_SUBMIT_DIR:-/nonexistent}/src/lib/env.sh" \
           "${DEFAULT_PROJ_ROOT:-/nonexistent}/src/lib/env.sh" \
           "${PROJ_ROOT:-/nonexistent}/src/lib/env.sh" \
           "${SCRIPT_DIR:-/nonexistent}/lib/env.sh" \
           "${SCRIPT_DIR:-/nonexistent}/../lib/env.sh" \
           "$(dirname -- "${BASH_SOURCE[0]:-$0}")/lib/env.sh" \
           "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../lib/env.sh"; do
    if [[ -f "$c" ]]; then printf '%s\n' "$c"; return 0; fi
  done
  return 1
}
_GWAS_ENV_LIB="$(_gwas_find_lib)" || {
  echo "ERROR: src/lib/env.sh not found. The tree is incompletely synced;" >&2
  echo "       rsync pipeline-strata/ again (see TRANSFER.md)." >&2
  exit 1
}
# shellcheck disable=SC1090
source "$_GWAS_ENV_LIB"

load_env_file "${PROJ_ROOT}/.env" || true
load_env_file "${PROJ_ROOT}/server.env" || true
mkdir -p "${PROJ_ROOT}/src/logs"

LOCAL_DIR="${GWAS_STAGING_DIR:-}"
REMOTE_NAME="${GWAS_TRANSFER_SUBDIR:-}"
LIST_ONLY=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) LOCAL_DIR="$2"; shift 2 ;;
    --name) REMOTE_NAME="$2"; shift 2 ;;
    --list) LIST_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '10,30p' "${BASH_SOURCE[0]:-$0}"
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

SFTP_HOST="${GWAS_TRANSFER_HOST:?GWAS_TRANSFER_HOST not set}"
SFTP_USER="${GWAS_TRANSFER_USER:?GWAS_TRANSFER_USER not set}"
SFTP_PASS="${GWAS_TRANSFER_PASS:?GWAS_TRANSFER_PASS not set}"
REMOTE_DIR="${GWAS_TRANSFER_REMOTE_DIR:?GWAS_TRANSFER_REMOTE_DIR not set}"

if [[ "$LIST_ONLY" -eq 0 ]]; then
  [[ -n "$LOCAL_DIR" ]] || { echo "ERROR: --dir is required (or set GWAS_STAGING_DIR in .env)" >&2; exit 1; }
  [[ -d "$LOCAL_DIR" ]] || { echo "ERROR: Directory not found: $LOCAL_DIR" >&2; exit 1; }
  # Default the remote folder to the local directory's name, but it is usually
  # worth naming explicitly: the local staging directory is named for where it
  # sits on scratch, not for what the collaborator should see.
  [[ -n "$REMOTE_NAME" ]] || REMOTE_NAME="$(basename -- "${LOCAL_DIR%/}")"
  # Recursive: the staged tree is <seqid>/<group>.tsv.gz, so a non-recursive
  # count reports only the manifest at the top level.
  N_FILES="$(find "$LOCAL_DIR" \( -type f -o -type l \) | wc -l | tr -d ' ')"
  echo "Local:  $LOCAL_DIR  (${N_FILES} file(s), recursive)"
  echo "Remote: ${SFTP_USER}@${SFTP_HOST}:${REMOTE_DIR%/}/${REMOTE_NAME}"
  [[ "$DRY_RUN" -eq 1 ]] && echo "DRY RUN: nothing will be uploaded and no connection is made."
fi

# Secrets and paths go to Python through the environment, not interpolated into
# the heredoc: a password containing a quote or backslash would otherwise break
# the script or inject into it.
export XFER_HOST="$SFTP_HOST" XFER_USER="$SFTP_USER" XFER_PASS="$SFTP_PASS"
export XFER_REMOTE="$REMOTE_DIR" XFER_LOCAL="${LOCAL_DIR:-}" XFER_NAME="${REMOTE_NAME:-}"
export XFER_LIST="$LIST_ONLY" XFER_DRY="$DRY_RUN"

python3 - <<'PYEOF'
import os, pathlib, sys

host       = os.environ["XFER_HOST"]
user       = os.environ["XFER_USER"]
password   = os.environ["XFER_PASS"]
remote_dir = os.environ["XFER_REMOTE"].rstrip("/")
local_dir  = pathlib.Path(os.environ["XFER_LOCAL"]) if os.environ.get("XFER_LOCAL") else None
name       = os.environ.get("XFER_NAME", "")
list_only  = os.environ["XFER_LIST"] == "1"
dry_run    = os.environ["XFER_DRY"] == "1"

dest_root = f"{remote_dir}/{name}" if name else remote_dir

# The plan is built from the local tree alone, so --dry-run needs no connection
# and no credentials to be exercised.
if not list_only:
    files = sorted(p for p in local_dir.rglob("*") if p.is_file() or p.is_symlink())
    if not files:
        sys.exit(f"ERROR: no files found under {local_dir}")
    plan = [(f, f"{dest_root}/{f.relative_to(local_dir).as_posix()}") for f in files]

    if dry_run:
        for f, dest in plan:
            print(f"  {f.relative_to(local_dir).as_posix()}  ->  {dest}")
        dirs = sorted({str(pathlib.PurePosixPath(d).parent) for _, d in plan})
        print(f"\n{len(plan)} file(s) into {len(dirs)} remote director(ies) under {dest_root}")
        print("DRY RUN: nothing uploaded.")
        sys.exit(0)

import paramiko

transport = paramiko.Transport((host, 22))
transport.connect(username=user, password=password)
sftp = paramiko.SFTPClient.from_transport(transport)

try:
    if list_only:
        def walk(path, depth=0):
            try:
                entries = sorted(sftp.listdir_attr(path), key=lambda e: e.filename)
            except IOError:
                print(f"  (cannot list {path})")
                return 0
            n = 0
            for e in entries:
                child = f"{path}/{e.filename}"
                if e.st_mode is not None and (e.st_mode & 0o40000):
                    print("  " * depth + f"{e.filename}/")
                    n += walk(child, depth + 1)
                else:
                    print("  " * depth + f"{e.filename:<30} {e.st_size:>14,}")
                    n += 1
            return n
        print(f"Remote: {dest_root}")
        total = walk(dest_root)
        print(f"\n{total} file(s)")
    else:
        # SFTP has no mkdir -p, so build each level and ignore "already exists".
        made = set()
        def ensure_dir(path):
            if path in made or path in ("", "/"):
                return
            parts = [p for p in path.strip("/").split("/") if p]
            cur = ""
            for part in parts:
                cur = f"{cur}/{part}"
                if cur in made:
                    continue
                try:
                    sftp.mkdir(cur)
                except OSError:
                    pass
                made.add(cur)

        for f, dest in plan:
            ensure_dir(str(pathlib.PurePosixPath(dest).parent))
            print(f"  {f.relative_to(local_dir).as_posix()}")
            sftp.put(str(f.resolve()), dest)
        print(f"\nTransfer complete: {len(plan)} file(s) -> {dest_root}")
finally:
    sftp.close()
    transport.close()
PYEOF
