# src/lib/env.sh — shared .env loader.
#
# Sourced by the src/*.sh wrappers. The single rule it enforces:
#
#   an exported variable ALWAYS beats the value in .env
#
# A server.env beside the .env is read after it, so its real values replace the
# tracked placeholders. The same rule applies: an exported variable beats both.
#
# The obvious `set -a; source .env; set +a` gets this backwards. It clobbers the
# caller's environment, so `GWAS_RESULTS_DIR=... sbatch src/006_staging.sh`
# silently ran against the .env path instead and failed on a directory the user
# never asked for. .env is a file of DEFAULTS; anything explicit outranks it.
#
# Only names actually assigned in the .env are considered, and every one of them
# is a GWAS_* setting, so this never shadows PATH or the SLURM_* variables.

# load_env_file <path>
#   Returns 0 if the file was found and loaded, 1 if it does not exist, so
#   callers can search a list of candidates:  load_env_file "$f" && break
load_env_file() {
  local env_file="$1"
  [[ -f "$env_file" ]] || return 1
  local server_env
  server_env="$(dirname -- "$env_file")/server.env"
  [[ -f "$server_env" ]] || server_env=""

  # Snapshot the values the caller already has for names this file assigns.
  local -a _preset=()
  local _name
  while IFS= read -r _name; do
    if [[ -n "${!_name+x}" ]]; then
      _preset+=("${_name}=${!_name}")
    fi
  done < <(
    sed -nE 's/^[[:space:]]*(export[[:space:]]+)?([A-Za-z_][A-Za-z0-9_]*)=.*/\2/p' "$env_file" ${server_env:+"$server_env"} \
      | sort -u
  )

  set -a
  # shellcheck disable=SC1090
  source "$env_file"
  if [[ -n "$server_env" ]]; then
    # shellcheck disable=SC1090
    source "$server_env"
  fi
  set +a

  # Put the caller's values back over the top.
  local _kv
  for _kv in "${_preset[@]+"${_preset[@]}"}"; do
    export "${_kv%%=*}=${_kv#*=}"
  done

  GWAS_ENV_FILE_LOADED="$env_file"
  return 0
}

# Kept as an alias: several scripts already call this name.
source_env_file() { load_env_file "$1" || return 0; }
