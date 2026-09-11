#!/usr/bin/env python3

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


ID_ALIASES = {"#iid", "iid", "sampleid", "sample_id"}
COLUMN_ALIASES = {
    "pheno": ("pheno1", "phenotype", "phenotype1"),
}


def load_dot_env(env_file: str | Path = ".env") -> None:
    env_path = Path(env_file)
    if not env_path.exists():
        return

    for raw_line in env_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        os.environ[key.strip()] = value.strip().strip("'\"")


def load_project_dotenv(project_root: str | Path | None = None) -> Path | None:
    root = project_root or os.getenv("GWAS_PROJECT_ROOT", "")
    search_paths: list[Path] = []

    if root:
        search_paths.append(Path(root))

    search_paths.extend([Path("."), Path(".."), Path("../..")])

    seen: set[Path] = set()
    for candidate in search_paths:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)

        env_path = resolved / ".env"
        if env_path.exists():
            load_dot_env(env_path)
            return env_path

    return None


def get_env(key: str, default: str | None = None) -> str | None:
    value = os.getenv(key, "")
    return default if value == "" else value


def resolve_column_name(columns: Iterable[str], requested: str | None) -> str | None:
    if requested is None:
        return None

    cleaned = requested.strip()
    if cleaned == "":
        return None

    lookup = {str(column).strip().lower(): str(column) for column in columns}

    direct_match = lookup.get(cleaned.lower())
    if direct_match is not None:
        return direct_match

    for alias in COLUMN_ALIASES.get(cleaned.lower(), ()):
        alias_match = lookup.get(alias.lower())
        if alias_match is not None:
            return alias_match

    return None


def _guess_separator(line: str) -> str:
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return r"\s+"


def read_table_auto(path: str | Path, dtype=None, header="infer") -> pd.DataFrame:
    file_path = Path(path)
    first_content_line = None

    with file_path.open() as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if stripped:
                first_content_line = stripped
                break

    if first_content_line is None:
        raise ValueError(f"Input file is empty: {file_path}")

    sep = _guess_separator(first_content_line)
    kwargs = {"sep": sep, "dtype": dtype, "header": header}
    if sep == r"\s+":
        kwargs["engine"] = "python"
    else:
        kwargs["low_memory"] = False

    return pd.read_csv(file_path, **kwargs)


def read_psam_file(path: str | Path) -> pd.DataFrame:
    file_path = Path(path)
    header_row = None

    with file_path.open() as handle:
        for index, raw_line in enumerate(handle):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("##"):
                continue
            header_row = index
            break

    if header_row is None:
        raise ValueError(f"Could not locate PSAM header in: {file_path}")

    return pd.read_csv(
        file_path,
        sep=r"\s+",
        engine="python",
        skiprows=header_row,
        dtype=str,
    )


def _coerce_non_identifier_columns(df: pd.DataFrame) -> pd.DataFrame:
    protected = {"FID", "#FID", "IID", "#IID", "PID", "MID"}

    for column in df.columns:
        if column in protected:
            continue

        series = df[column]
        if not pd.api.types.is_object_dtype(series) and not pd.api.types.is_string_dtype(series):
            continue

        non_empty = series.notna() & (series.astype(str).str.strip() != "")
        if not non_empty.any():
            continue

        converted = pd.to_numeric(series, errors="coerce")
        if converted[non_empty].notna().all():
            df[column] = converted

    return df


def read_sample_file(path: str | Path) -> pd.DataFrame:
    file_path = Path(path)
    suffix = file_path.suffix.lower()

    if suffix == ".psam":
        df = read_psam_file(file_path)
    else:
        df = read_table_auto(file_path, dtype=str)

    df = _coerce_non_identifier_columns(df)
    df = standardize_ids(df)

    # PLINK2 sample files sometimes expose the default phenotype column as
    # PHENO1 rather than PHENO. Normalize it so downstream grouping logic can
    # consistently request PHENO across .psam inputs.
    phenotype_column = resolve_column_name(df.columns, "PHENO")
    if phenotype_column is not None and phenotype_column != "PHENO":
        df = df.rename(columns={phenotype_column: "PHENO"})

    return df


def standardize_ids(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        if column.lower() in ID_ALIASES:
            rename_map[column] = "IID"
            break

    if rename_map:
        df = df.rename(columns=rename_map)

    return df


def read_path_list(raw_value: str | None) -> list[str] | None:
    if raw_value is None:
        return None

    cleaned = raw_value.strip()
    if cleaned == "":
        return None

    path_candidate = Path(cleaned)
    if path_candidate.exists():
        values = [line.strip() for line in path_candidate.read_text().splitlines()]
        return [value for value in values if value]

    return [part.strip() for part in cleaned.split(",") if part.strip()]


def ensure_string_ids(df: pd.DataFrame, columns: Iterable[str] = ("FID", "IID")) -> pd.DataFrame:
    for column in columns:
        if column in df.columns:
            df[column] = df[column].astype("string")
    return df


def make_safe_name(name: str) -> str:
    safe = re.sub(r"[^0-9A-Za-z_.]", ".", str(name))
    safe = re.sub(r"\.+", ".", safe)

    if safe == "":
        safe = "X"

    if re.match(r"^[^A-Za-z.]|^\.[0-9]", safe):
        safe = f"X{safe}"

    return safe


def make_safe_names(names: Iterable[str]) -> list[str]:
    used: set[str] = set()
    safe_names: list[str] = []

    for name in names:
        base = make_safe_name(name)
        candidate = base
        suffix = 1

        while candidate in used:
            candidate = f"{base}.{suffix}"
            suffix += 1

        used.add(candidate)
        safe_names.append(candidate)

    return safe_names


load_project_dotenv()
