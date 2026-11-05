#!/bin/bash
#SBATCH --job-name=001_source
#SBATCH --output=src/logs/001_source.out
#SBATCH --error=src/logs/001_source.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# Purpose: Set up the lightweight Nextflow environment and bundled pipeline tools.

set -euo pipefail

SCRIPT_SOURCE="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(dirname -- "$SCRIPT_SOURCE")" && pwd)"
# Under SLURM the job runs from the slurmd spool dir and the script is a copy
# there, so derive the project root from SLURM_SUBMIT_DIR (the directory sbatch
# was launched from) rather than from the script/PWD location.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR}"
else
  DEFAULT_PROJ_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
ENV_NAME="nf_master"

source_env_file() {
  local env_file="$1"
  if [[ -f "$env_file" ]]; then
    set -a
    source "$env_file"
    set +a
  fi
}

source_env_file "${PWD}/.env"
if [[ "$DEFAULT_PROJ_ROOT" != "$PWD" ]]; then
  source_env_file "${DEFAULT_PROJ_ROOT}/.env"
fi

# Resolve PROJ_ROOT only after .env is sourced, so GWAS_PROJECT_ROOT is honoured.
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$DEFAULT_PROJ_ROOT}"
# SLURM opens the #SBATCH log paths before a job starts.
mkdir -p "${PROJ_ROOT}/src/logs"
PIPELINE_DIR="${PROJ_ROOT}/pipeline"
# Prebuilt pipeline conda environments live here (one prefix env per spec in
# pipeline/envs/). The Nextflow modules point their `conda` directive at
# ${projectDir}/../env/<name>, i.e. this directory, so no env is ever built on
# scratch at run time.
ENV_DIR="${PROJ_ROOT}/env"
TOOLS_DIR="${PIPELINE_DIR}/tools"
METAL_URL="https://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz"
FORCE=0
SETUP_ENV=1
SETUP_TOOLS=1

choose_conda_cmd() {
  if command -v mamba >/dev/null 2>&1; then
    printf '%s\n' "mamba"
    return 0
  fi

  if command -v conda >/dev/null 2>&1; then
    printf '%s\n' "conda"
    return 0
  fi

  echo "ERROR: Neither conda nor mamba is available in your PATH." >&2
  exit 1
}

setup_nextflow_env() {
  local conda_cmd
  conda_cmd="$(choose_conda_cmd)"

  echo "=========================================="
  echo " Setting up Nextflow Master Environment"
  echo "=========================================="

  if "$conda_cmd" env list | awk 'NR > 2 {print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Environment '${ENV_NAME}' already exists. Skipping creation."
    return 0
  fi

  echo "Creating environment '${ENV_NAME}' with ${conda_cmd}..."
  "$conda_cmd" create -n "$ENV_NAME" -c conda-forge -c bioconda nextflow openjdk=21 python=3.9 -y

  echo "=========================================="
  echo "Environment ready. Activate with: 'conda activate ${ENV_NAME}'"
  echo "=========================================="

  echo "Installing Python dependencies (paramiko)..."
  pip install --quiet paramiko
}

# Pre-build every pipeline conda environment as a prefix env under ${ENV_DIR}
# (one per spec in pipeline/envs/). The Nextflow modules reference these via
# ${projectDir}/../env/<name>, and src/002_data-genetics.sh activates env/plink2,
# so no environment is ever built on scratch at pipeline run time (which was slow
# and unreliable on NFS, with builds timing out / failing mid-verify). Built once
# here, persisted in the project root, reused by every run. Use --force to rebuild.
setup_pipeline_envs() {
  local conda_cmd
  conda_cmd="$(choose_conda_cmd)"

  echo "=========================================="
  echo " Building pipeline conda environments"
  echo " Target: ${ENV_DIR}"
  echo "=========================================="

  if [[ ! -d "${PIPELINE_DIR}/envs" ]]; then
    echo "ERROR: env spec directory not found: ${PIPELINE_DIR}/envs" >&2
    exit 1
  fi

  # METAL is not a conda package on these channels; it is installed as a binary
  # by install_metal() into pipeline/tools/metal and referenced via
  # params.metal_binary, so it is never built as a conda env here.
  local skip_envs=" metal "

  mkdir -p "$ENV_DIR"
  local failed=() yml name prefix ready
  for yml in "${PIPELINE_DIR}"/envs/*.yml; do
    [[ -e "$yml" ]] || continue
    name="$(basename "$yml" .yml)"
    prefix="${ENV_DIR}/${name}"
    ready="${prefix}/.gwas_env_ready"

    if [[ "$skip_envs" == *" $name "* ]]; then
      echo "  [skip]    ${name} (provided as a binary, not a conda env)"
      continue
    fi

    # Skip only fully-built envs (marked by the readiness sentinel). A bare
    # directory means a previous build failed part-way; rebuild it from scratch
    # rather than trusting a partial env.
    if [[ -f "$ready" && "$FORCE" -ne 1 ]]; then
      echo "  [skip]    ${name} (already built; use --force to rebuild)"
      continue
    fi
    if [[ -d "$prefix" ]]; then
      echo "  [clean]   removing existing/partial ${prefix}"
      rm -rf "$prefix"
    fi

    echo "  [build]   ${name} -> ${prefix}"
    if "$conda_cmd" env create --prefix "$prefix" -f "$yml"; then
      touch "$ready"
      echo "  [ok]      ${name}"
    else
      echo "  [FAIL]    ${name}" >&2
      failed+=("$name")
    fi
  done

  if [[ ${#failed[@]} -gt 0 ]]; then
    echo "ERROR: failed to build environment(s): ${failed[*]}" >&2
    exit 1
  fi

  echo "All pipeline environments ready under ${ENV_DIR}."
}

install_metal() {
  echo "=========================================="
  echo " Installing Bundled METAL Binary"
  echo "=========================================="

  mkdir -p "$TOOLS_DIR"
  ARCHIVE="${TOOLS_DIR}/Linux-metal.tar.gz"

  if [[ -f "$ARCHIVE" && "$FORCE" -ne 1 ]]; then
    echo "Archive already exists: $ARCHIVE (use --force to re-download)"
  else
    echo "Downloading METAL from: $METAL_URL"
    if command -v curl >/dev/null 2>&1; then
      curl -L -o "$ARCHIVE" "$METAL_URL"
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$ARCHIVE" "$METAL_URL"
    else
      echo "ERROR: curl or wget is required to download METAL."
      exit 1
    fi
  fi

  TMPDIR="$(mktemp -d -p "$TOOLS_DIR")"
  tar -xzf "$ARCHIVE" -C "$TMPDIR"

  METAL_BIN="$(find "$TMPDIR" -type f -name metal -perm -u+x | head -n 1)"
  if [[ -z "$METAL_BIN" ]]; then
    METAL_BIN="$(find "$TMPDIR" -type f -name metal | head -n 1)"
  fi

  if [[ -z "$METAL_BIN" ]]; then
    echo "ERROR: Could not locate METAL binary in extracted archive."
    rm -rf "$TMPDIR"
    exit 1
  fi

  mkdir -p "${TOOLS_DIR}/metal"
  cp "$METAL_BIN" "${TOOLS_DIR}/metal/metal"
  chmod +x "${TOOLS_DIR}/metal/metal"

  rm -rf "$TMPDIR"

  echo "METAL installed to: ${TOOLS_DIR}/metal/metal"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --url) METAL_URL="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --env-only)
      SETUP_TOOLS=0
      shift
      ;;
    --tools-only)
      SETUP_ENV=0
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--env-name NAME] [--tools-dir PATH] [--url URL] [--force] [--env-only|--tools-only]"
      echo "Default env name: ${ENV_NAME}"
      echo "Default install location: ${TOOLS_DIR}"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ "$SETUP_ENV" -eq 1 ]]; then
  setup_nextflow_env
  setup_pipeline_envs
fi

if [[ "$SETUP_TOOLS" -eq 1 ]]; then
  install_metal
fi
