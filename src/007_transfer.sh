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
# Usage:
#   sbatch src/007_transfer.sh [--dir <path>]
#
# Read from .env:
#   GWAS_STAGING_DIR         Local directory to upload (default for --dir)
# Read from the untracked server.env:
#   GWAS_TRANSFER_HOST       SFTP hostname
#   GWAS_TRANSFER_USER       SFTP username
#   GWAS_TRANSFER_PASS       SFTP password
#   GWAS_TRANSFER_REMOTE_DIR Remote destination path

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="${GWAS_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}}"

if [[ -f "${PROJ_ROOT}/.env" ]]; then
  set -a; source "${PROJ_ROOT}/.env"; set +a
fi
if [[ -f "${PROJ_ROOT}/server.env" ]]; then
  set -a; source "${PROJ_ROOT}/server.env"; set +a
fi
mkdir -p "${PROJ_ROOT}/src/logs"

LOCAL_DIR="${GWAS_STAGING_DIR:-}"
LIST_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) LOCAL_DIR="$2"; shift 2 ;;
    --list) LIST_ONLY=1; shift ;;
    -h|--help)
      echo "Usage: bash src/007_transfer.sh [--dir <staged-directory>] [--list]"
      echo "       Default: GWAS_STAGING_DIR from server.env"
      echo "       --list  List files on the remote without uploading"
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

SFTP_HOST="${GWAS_TRANSFER_HOST:?GWAS_TRANSFER_HOST not set}"
SFTP_USER="${GWAS_TRANSFER_USER:?GWAS_TRANSFER_USER not set}"
SFTP_PASS="${GWAS_TRANSFER_PASS:?GWAS_TRANSFER_PASS not set}"
REMOTE_DIR="${GWAS_TRANSFER_REMOTE_DIR:?GWAS_TRANSFER_REMOTE_DIR not set}"

if [[ "$LIST_ONLY" -eq 1 ]]; then
  echo "Listing remote: ${SFTP_USER}@${SFTP_HOST}:${REMOTE_DIR}"
else
  [[ -n "$LOCAL_DIR" ]] || { echo "ERROR: --dir is required (or set GWAS_STAGING_DIR in server.env)" >&2; exit 1; }
  [[ -d "$LOCAL_DIR" ]] || { echo "ERROR: Directory not found: $LOCAL_DIR" >&2; exit 1; }
  N_FILES="$(find "$LOCAL_DIR" -maxdepth 1 \( -type f -o -type l \) | wc -l)"
  echo "Uploading ${N_FILES} file(s) from $LOCAL_DIR"
  echo "  -> ${SFTP_USER}@${SFTP_HOST}:${REMOTE_DIR}"
fi

# Secrets and paths go to Python through the environment, not interpolated into
# the heredoc: a password containing a quote or backslash would otherwise break
# the script or inject into it.
export XFER_HOST="$SFTP_HOST" XFER_USER="$SFTP_USER" XFER_PASS="$SFTP_PASS"
export XFER_REMOTE="$REMOTE_DIR" XFER_LOCAL="${LOCAL_DIR:-}" XFER_LIST="$LIST_ONLY"

python3 - <<'PYEOF'
import os, pathlib, paramiko

host       = os.environ["XFER_HOST"]
user       = os.environ["XFER_USER"]
password   = os.environ["XFER_PASS"]
remote_dir = os.environ["XFER_REMOTE"]
list_only  = os.environ["XFER_LIST"] == "1"

transport = paramiko.Transport((host, 22))
transport.connect(username=user, password=password)
sftp = paramiko.SFTPClient.from_transport(transport)

if list_only:
    entries = sorted(sftp.listdir_attr(remote_dir), key=lambda e: e.filename)
    print(f"{'Size':>12}  {'Name'}")
    print("-" * 50)
    for e in entries:
        print(f"{e.st_size:>12,}  {e.filename}")
    print(f"\n{len(entries)} file(s) in {remote_dir}")
else:
    local_dir = pathlib.Path(os.environ["XFER_LOCAL"])
    # Transfer the directory itself: create remote_dir/<local_dir name>/ and
    # upload the files into it (rather than flattening into remote_dir).
    dest_dir = remote_dir.rstrip("/") + "/" + local_dir.name
    for d in (remote_dir, dest_dir):
        try:
            sftp.mkdir(d)
        except OSError:
            pass  # already exists
    files = sorted(p for p in local_dir.iterdir() if p.is_file() or p.is_symlink())
    print(f"Remote target: {dest_dir}")
    for f in files:
        dest = dest_dir + "/" + f.name
        print(f"  {local_dir.name}/{f.name}")
        sftp.put(str(f.resolve()), dest)
    print(f"Transfer complete: {len(files)} file(s) uploaded to {dest_dir}")

sftp.close()
transport.close()
PYEOF
