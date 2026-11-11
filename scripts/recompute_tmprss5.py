#!/usr/bin/env python3
"""
Recompute Tmprss5 statistics from GSE103799 sciatic-nerve TPM values
(WT_PBS, C22_PBS, C22_ASO). Optionally verify CSV values against NCBI GEO.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import os
import re
import sys
import urllib.error
import urllib.request
from collections import defaultdict

from scipy import stats


CSV_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "gse103799_tmprss5_tpm.csv",
)

MAIN_GROUPS = ("WT_PBS", "C22_PBS", "C22_ASO")
GROUPS = MAIN_GROUPS + ("WT_ASO",)
SERIES_ACC = "GSE103799"
MATCH_TOL = 0.005

# Tmprss5 (Ensembl gene ENSMUSG00000032268) transcript IDs, retrieved 2026-11-10.
# Union of two sources; IDs are unversioned and matched with any version suffix stripped.
#   [r102] Ensembl release 102 (GRCm38) GTF, transcript lines with
#          gene_id "ENSMUSG00000032268" (gene_name "Tmprss5"):
#          https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
#   [rest] Current Ensembl REST lookup:
#          https://rest.ensembl.org/lookup/symbol/mus_musculus/Tmprss5?expand=1
# The GEO records cite "mouse ensembl build 81"; the release 81 GTF
# (https://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz)
# lists the same seven transcripts as release 102.
TMPRSS5_TRANSCRIPTS = (
    "ENSMUST00000070390",  # r102, rest  (protein_coding)
    "ENSMUST00000165088",  # r102, rest  (protein_coding)
    "ENSMUST00000166272",  # r102, rest  (nonsense_mediated_decay)
    "ENSMUST00000167095",  # r102, rest  (protein_coding)
    "ENSMUST00000170246",  # r102, rest  (protein_coding)
    "ENSMUST00000170426",  # r102, rest  (nonsense_mediated_decay)
    "ENSMUST00000171217",  # r102, rest  (retained_intron)
    "ENSMUST00000513301",  # rest only
    "ENSMUST00000513302",  # rest only
    "ENSMUST00000513303",  # rest only
)
USER_AGENT = "cmt1a-tmprss5-reanalysis/1.0 (GEO verification; contact via GitHub)"


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def sample_sd(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    m = mean(values)
    return math.sqrt(sum((x - m) ** 2 for x in values) / (n - 1))


def fmt_p(p: float) -> str:
    if p >= 1e-4:
        return f"{p:.5f}"
    return f"{p:.1e}"


def load_rows(path: str) -> list[dict]:
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["group"].strip() not in GROUPS:
                raise ValueError(
                    f"{row['sample_id']}: unknown group {row['group']!r}; expected one of {GROUPS}"
                )
            rows.append(
                {
                    "sample_id": row["sample_id"].strip(),
                    "group": row["group"].strip(),
                    "tmprss5_tpm": float(row["tmprss5_tpm"]),
                }
            )
    return rows


def load_groups(rows: list[dict]) -> dict[str, list[float]]:
    groups: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        groups[row["group"]].append(row["tmprss5_tpm"])
    return groups


def print_group_line(group: str, values: list[float]) -> None:
    print(
        f"  {group:8s}  n={len(values)}  mean={mean(values):.2f}  SD={sample_sd(values):.2f}  "
        f"min-max={min(values):.2f}-{max(values):.2f}"
    )


def print_tests(a: list[float], b: list[float]) -> None:
    welch_raw = stats.ttest_ind(a, b, equal_var=False)
    log_a = [math.log2(x) for x in a]
    log_b = [math.log2(x) for x in b]
    welch_log = stats.ttest_ind(log_a, log_b, equal_var=False)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    print(f"    Welch t-test on raw TPM:            p={fmt_p(welch_raw.pvalue)}")
    print(f"    Welch t-test on log2(TPM):          p={fmt_p(welch_log.pvalue)}")
    print(f"    Mann-Whitney U:                     p={fmt_p(mw.pvalue)}")


def print_analysis(rows: list[dict]) -> None:
    groups = load_groups(rows)
    print("=" * 64)
    print("GSE103799 Tmprss5 sciatic nerve (15-week subset)")
    print("=" * 64)
    print()

    summary: dict[str, dict] = {}
    print("Group summary (n, mean, sample SD, min-max):")
    for group in MAIN_GROUPS:
        values = groups[group]
        summary[group] = {"mean": mean(values), "values": values}
        print_group_line(group, values)

    wt = summary["WT_PBS"]["mean"]
    c22 = summary["C22_PBS"]["mean"]
    aso = summary["C22_ASO"]["mean"]
    wt_vals = summary["WT_PBS"]["values"]
    c22_vals = summary["C22_PBS"]["values"]
    aso_vals = summary["C22_ASO"]["values"]

    fold_over_wt = c22 / wt
    disease_lfc = math.log2(c22 / wt)
    gap = c22 - wt

    ratio = aso / c22
    treatment_lfc = math.log2(aso / c22)
    shift = c22 - aso
    reduction_from_untreated = shift / c22
    share_of_gap_closed = shift / gap
    treated_over_wt = aso / wt

    print()
    print("Disease contrast (C22_PBS vs WT_PBS):")
    print(f"  Fold over WT:                         {fold_over_wt:.2f}x")
    print(f"  log2FC:                               {disease_lfc:+.2f}")
    print(f"  Gap over WT:                          {gap:.2f} TPM")

    print()
    print("Treatment contrast (C22_ASO vs C22_PBS):")
    print(f"  Ratio (C22_ASO / C22_PBS):             {ratio:.2f}")
    print(f"  log2FC:                               {treatment_lfc:+.2f}")
    print(f"  Shift (C22_PBS - C22_ASO):             {shift:.2f} TPM")
    print(
        "  (a) Reduction from untreated "
        f"(C22_PBS - C22_ASO) / C22_PBS:          {100.0 * reduction_from_untreated:.1f}%"
    )
    print(
        "  (b) Share of the gap over WT closed "
        f"(C22_PBS - C22_ASO) / (C22_PBS - WT_PBS): {100.0 * share_of_gap_closed:.1f}%"
    )
    print(f"  Treated remain above WT:              {treated_over_wt:.2f}x")

    contrasts = [
        ("Disease (C22_PBS vs WT_PBS)", c22_vals, wt_vals),
        ("Treatment (C22_ASO vs C22_PBS)", aso_vals, c22_vals),
        ("Treated vs WT (C22_ASO vs WT_PBS)", aso_vals, wt_vals),
    ]

    print()
    print("Two-sided p-values (one gene; no multiple-testing correction):")
    for label, a, b in contrasts:
        print(f"  {label}:")
        print_tests(a, b)

    aso_lo, aso_hi = min(aso_vals), max(aso_vals)
    over_100 = sorted(v for v in c22_vals if v > 100)
    in_aso_range = sorted(v for v in c22_vals if aso_lo <= v <= aso_hi)

    print()
    print("C22_PBS distribution notes:")
    if over_100:
        print(
            "  Values exceeding 100 TPM: "
            + ", ".join(f"{v:.2f}" for v in over_100)
        )
    else:
        print("  Values exceeding 100 TPM: none")
    print(
        f"  C22_ASO range: {aso_lo:.2f}-{aso_hi:.2f} TPM; "
        "C22_PBS values inside that range: "
        + (", ".join(f"{v:.2f}" for v in in_aso_range) if in_aso_range else "none")
    )
    print()

    wt_aso_vals = groups["WT_ASO"]
    wt_aso = mean(wt_aso_vals)
    print("ASO control in wild-type (WT_ASO vs WT_PBS):")
    print_group_line("WT_ASO", wt_aso_vals)
    print(f"  Ratio (WT_ASO / WT_PBS):               {wt_aso / wt:.2f}")
    print(f"  log2FC:                               {math.log2(wt_aso / wt):+.2f}")
    print("  Two-sided p-values:")
    print_tests(wt_aso_vals, wt_vals)
    print()


def http_get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def geo_soft_text(accession: str) -> str:
    # view=full retains supplementary_file_N and extract_protocol lines needed here.
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={accession}&targ=self&form=text&view=full"
    )
    return http_get(url).decode("utf-8", errors="replace")


def soft_field_values(text: str, field: str) -> list[str]:
    pattern = re.compile(rf"^{re.escape(field)}\s*=\s*(.*)$", re.MULTILINE)
    return [m.group(1).strip() for m in pattern.finditer(text)]


def soft_field_values_prefix(text: str, field_prefix: str) -> list[str]:
    """Match exact field or numbered GEO variants (e.g. !Sample_supplementary_file_1)."""
    pattern = re.compile(
        rf"^{re.escape(field_prefix)}(?:_\d+)?\s*=\s*(.*)$", re.MULTILINE
    )
    return [m.group(1).strip() for m in pattern.finditer(text)]


def soft_first(text: str, field: str) -> str | None:
    values = soft_field_values(text, field)
    return values[0] if values else None


def ftp_to_https(url: str) -> str:
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://") :]
    return url


def open_maybe_gzip(raw: bytes, url: str) -> io.TextIOBase:
    lower = url.lower()
    if lower.endswith(".gz") or raw[:2] == b"\x1f\x8b":
        return io.TextIOWrapper(io.BytesIO(gzip.decompress(raw)), encoding="utf-8", errors="replace")
    return io.StringIO(raw.decode("utf-8", errors="replace"))


def detect_delimiter(header_line: str) -> str:
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return "\t"


def find_tpm_column(header: list[str], override: str | None) -> tuple[int | None, str]:
    if override is not None:
        for i, name in enumerate(header):
            if name == override or name.lower() == override.lower():
                return i, name
        return None, f"override '{override}' not found in header={header[:12]}"

    for i, name in enumerate(header):
        if name.strip().upper() == "TPM":
            return i, name
    if len(header) == 2:
        return 1, header[1]
    return None, f"no TPM column detected; header={header[:20]}"


def norm_id(identifier: str) -> str:
    return identifier.strip().split("|")[0].split(".")[0].lower()


def extract_gene_tpm_from_table(
    text_handle: io.TextIOBase,
    gene_ids: list[str],
    tpm_column: str | None,
) -> tuple[float | None, str, set[str]]:
    lines = text_handle.read().splitlines()
    if not lines:
        return None, "empty file", set()
    # Skip comment lines if present
    start = 0
    while start < len(lines) and lines[start].startswith("#"):
        start += 1
    if start >= len(lines):
        return None, "only comment lines", set()
    header_line = lines[start]
    delim = detect_delimiter(header_line)
    header = [c.strip() for c in header_line.split(delim)]
    # If first row does not look like a header (no TPM / non-numeric gene id),
    # treat as two-column gene/TPM without header.
    looks_like_header = any(h.upper() == "TPM" for h in header) or (
        tpm_column is not None and any(h.lower() == tpm_column.lower() for h in header)
    )
    data_start = start + 1
    if not looks_like_header and len(header) == 2:
        # synthetic header; include first line as data
        header = ["gene_id", "TPM"]
        data_start = start
    col_idx, col_note = find_tpm_column(header, tpm_column)
    if col_idx is None:
        first_ids = []
        for line in lines[data_start : data_start + 5]:
            parts = line.split(delim)
            if parts:
                first_ids.append(parts[0].strip())
        return None, f"{col_note}; first identifiers={first_ids}", set()

    wanted = {norm_id(g) for g in gene_ids}
    total = 0.0
    hits = 0
    found: set[str] = set()
    first_ids: list[str] = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        parts = line.split(delim)
        if not parts:
            continue
        gene = parts[0].strip()
        if len(first_ids) < 8:
            first_ids.append(gene)
        if norm_id(gene) not in wanted:
            continue
        if col_idx >= len(parts):
            return None, (
                f"row for {gene} lacks column index {col_idx}; "
                f"header={header[:12]}; first identifiers={first_ids}"
            ), found
        try:
            total += float(parts[col_idx])
            hits += 1
            found.add(norm_id(gene))
        except ValueError:
            return None, (
                f"non-numeric TPM for {gene}: {parts[col_idx]!r}; "
                f"header={header[:12]}"
            ), found
    if hits == 0:
        return None, (
            f"no row matched any of {len(gene_ids)} listed IDs; header={header[:12]}; "
            f"first identifiers={first_ids}"
        ), found
    return total, f"summed {hits} row(s) from column '{col_note}'", found


def extract_from_series_matrix(
    text_handle: io.TextIOBase,
    sample_id: str,
    sample_title: str | None,
    gene_ids: list[str],
) -> tuple[float | None, str, set[str]]:
    lines = text_handle.read().splitlines()
    # GEO series matrix: !Sample_geo_accession / !Sample_title then table
    accessions: list[str] = []
    titles: list[str] = []
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            accessions = [c.strip().strip('"') for c in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [c.strip().strip('"') for c in line.split("\t")[1:]]

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("ID_REF") or line.startswith('"ID_REF"'):
            header_idx = i
            break
    if header_idx is None:
        preview = [ln[:120] for ln in lines[:15]]
        return None, f"no ID_REF table header in series matrix; preview={preview}", set()

    header = [c.strip().strip('"') for c in lines[header_idx].split("\t")]
    col = None
    if sample_id in header:
        col = header.index(sample_id)
    elif accessions and sample_id in accessions:
        # matrix columns align with sample order after ID_REF
        col = accessions.index(sample_id) + 1
    elif sample_title and titles and sample_title in titles:
        col = titles.index(sample_title) + 1
    elif sample_title and sample_title in header:
        col = header.index(sample_title)

    if col is None:
        return None, (
            f"could not map {sample_id}/{sample_title!r} to a matrix column; "
            f"header[:8]={header[:8]}; accessions[:8]={accessions[:8]}; "
            f"titles[:4]={titles[:4]}"
        ), set()

    wanted = {norm_id(g) for g in gene_ids}
    total = 0.0
    hits = 0
    found: set[str] = set()
    first_ids: list[str] = []
    for line in lines[header_idx + 1 :]:
        if line.startswith("!") or not line.strip():
            continue
        parts = [c.strip().strip('"') for c in line.split("\t")]
        if not parts:
            continue
        gene = parts[0]
        if len(first_ids) < 8:
            first_ids.append(gene)
        if norm_id(gene) not in wanted:
            continue
        if col >= len(parts):
            return None, f"matrix row short for {gene}; header[:8]={header[:8]}", found
        try:
            total += float(parts[col])
            hits += 1
            found.add(norm_id(gene))
        except ValueError:
            return None, f"non-numeric matrix value for {gene}: {parts[col]!r}", found
    if hits == 0:
        return None, (
            f"no row matched any of {len(gene_ids)} listed IDs in series matrix; "
            f"first identifiers={first_ids}"
        ), found
    return total, (
        f"series matrix column {col} ({header[col] if col < len(header) else 'index'})"
    ), found


def parse_sample_soft(text: str) -> dict:
    return {
        "title": soft_first(text, "!Sample_title"),
        "characteristics": soft_field_values(text, "!Sample_characteristics_ch1"),
        "supplementary": soft_field_values_prefix(text, "!Sample_supplementary_file"),
        "platform": soft_first(text, "!Sample_platform_id"),
        "library_strategy": soft_first(text, "!Sample_library_strategy"),
        "library_source": soft_first(text, "!Sample_library_source"),
        "library_selection": soft_first(text, "!Sample_library_selection"),
        "instrument": soft_first(text, "!Sample_instrument_model"),
        "extract_protocol": soft_field_values(text, "!Sample_extract_protocol_ch1"),
        "treatment_protocol": soft_field_values(text, "!Sample_treatment_protocol_ch1"),
        "data_processing": soft_field_values(text, "!Sample_data_processing"),
    }


def fetch_geo(rows: list[dict], gene_ids: list[str], tpm_column: str | None) -> int:
    print("=" * 64)
    print(f"GEO verification against {SERIES_ACC}")
    print("=" * 64)
    print()

    series_text = geo_soft_text(SERIES_ACC)
    series_supp = soft_field_values_prefix(series_text, "!Series_supplementary_file")
    series_summary = soft_field_values(series_text, "!Series_summary")
    series_overall = soft_field_values(series_text, "!Series_overall_design")
    print("Series-level notes (for documentation):")
    for key, values in [
        ("!Series_supplementary_file", series_supp),
        ("!Series_summary", series_summary[:3]),
        ("!Series_overall_design", series_overall[:3]),
    ]:
        if values:
            for v in values:
                print(f"  {key}: {v}")
        else:
            print(f"  {key}: (not stated)")
    print()

    listed = {norm_id(g): g for g in gene_ids}
    print(f"Matching {len(listed)} listed IDs: {', '.join(gene_ids)}")
    print()

    file_cache: dict[str, bytes] = {}

    def get_bytes(url: str) -> bytes:
        url = ftp_to_https(url)
        if url not in file_cache:
            print(f"  downloading: {url}")
            file_cache[url] = http_get(url)
        return file_cache[url]

    # Prefer a series matrix / quantification file if samples lack per-sample files
    matrix_urls = [
        u for u in series_supp
        if "matrix" in u.lower() or u.lower().endswith(".txt.gz") or "tpm" in u.lower()
    ]

    failures = 0
    meta_seen = False

    for row in rows:
        sid = row["sample_id"]
        csv_val = row["tmprss5_tpm"]
        csv_group = row["group"]
        status = "FAILED"
        geo_val: float | None = None
        detail = ""
        title = ""
        chars = ""
        found: set[str] = set()
        try:
            soft = geo_soft_text(sid)
            meta = parse_sample_soft(soft)
            title = meta["title"] or ""
            chars = "; ".join(meta["characteristics"]) if meta["characteristics"] else ""
            if not meta_seen:
                print("Sample-level library / platform / processing (first sample):")
                for key in (
                    "platform",
                    "library_strategy",
                    "library_source",
                    "library_selection",
                    "instrument",
                ):
                    val = meta.get(key)
                    print(f"  {key}: {val if val else '(not stated)'}")
                for ep in meta.get("extract_protocol") or []:
                    print(f"  extract_protocol: {ep}")
                if not meta.get("extract_protocol"):
                    print("  extract_protocol: (not stated)")
                for tp in meta.get("treatment_protocol") or []:
                    print(f"  treatment_protocol: {tp}")
                if meta["data_processing"]:
                    for dp in meta["data_processing"]:
                        print(f"  data_processing: {dp}")
                else:
                    print("  data_processing: (not stated)")
                print()
                meta_seen = True

            geo_val = None
            note = ""
            if meta["supplementary"]:
                # Try each supplementary file until a gene TPM is extracted
                last_err = ""
                for supp in meta["supplementary"]:
                    raw = get_bytes(supp)
                    handle = open_maybe_gzip(raw, supp)
                    geo_val, note, found = extract_gene_tpm_from_table(
                        handle, gene_ids, tpm_column
                    )
                    if geo_val is not None:
                        break
                    last_err = note
                if geo_val is None:
                    detail = last_err or "no usable per-sample supplementary file"
            else:
                # Fall back to series-level supplementary / matrix
                if not matrix_urls and not series_supp:
                    detail = "no per-sample supplementary file and no series supplementary file"
                else:
                    urls = matrix_urls or series_supp
                    last_err = ""
                    for url in urls:
                        raw = get_bytes(url)
                        # Try as series matrix first
                        text = open_maybe_gzip(raw, url).read()
                        geo_val, note, found = extract_from_series_matrix(
                            io.StringIO(text), sid, title, gene_ids
                        )
                        if geo_val is None:
                            geo_val, note, found = extract_gene_tpm_from_table(
                                io.StringIO(text), gene_ids, tpm_column
                            )
                        if geo_val is not None:
                            break
                        last_err = note
                    if geo_val is None:
                        detail = last_err or "series supplementary files unusable"

            if geo_val is None:
                status = "FAILED"
                failures += 1
            elif abs(geo_val - csv_val) <= MATCH_TOL:
                status = "OK"
                detail = note
            else:
                status = "MISMATCH"
                detail = note
                failures += 1
        except Exception as exc:  # noqa: BLE001 - report and continue
            status = "FAILED"
            detail = f"{type(exc).__name__}: {exc}"
            failures += 1

        geo_disp = f"{geo_val:.4f}" if geo_val is not None else "NA"
        absent = [listed[k] for k in listed if k not in found]
        coverage = f"found {len(listed) - len(absent)}/{len(listed)} listed IDs"
        if absent:
            coverage += f" (absent: {', '.join(absent)})"
        print(
            f"{status:9s} {sid}  CSV={csv_val:.2f}  GEO={geo_disp}  "
            f"CSV_group={csv_group}  title={title!r}  chars={chars!r}  {coverage}"
            + (f"  note={detail}" if detail else "")
        )

    print()
    if failures:
        print(f"GEO check failed: {failures} sample(s) not OK.")
        return 1
    print("GEO check passed: all samples OK.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Recompute GSE103799 Tmprss5 stats and optionally verify against GEO."
    )
    p.add_argument(
        "--fetch-geo",
        action="store_true",
        help="Verify each CSV TPM against NCBI GEO supplementary files.",
    )
    p.add_argument(
        "--gene-id",
        action="append",
        dest="gene_ids",
        help=(
            "Identifier to match on the first column, version suffix ignored "
            "(repeatable). Overrides the default TMPRSS5_TRANSCRIPTS list."
        ),
    )
    p.add_argument(
        "--tpm-column",
        default=None,
        help="Override TPM column name in per-sample tables.",
    )
    p.add_argument(
        "--csv",
        default=CSV_PATH,
        help="Path to CSV (default: data/gse103799_tmprss5_tpm.csv).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = load_rows(args.csv)
    if args.fetch_geo:
        gene_ids = args.gene_ids if args.gene_ids else list(TMPRSS5_TRANSCRIPTS)
        return fetch_geo(rows, gene_ids, args.tpm_column)
    print_analysis(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
