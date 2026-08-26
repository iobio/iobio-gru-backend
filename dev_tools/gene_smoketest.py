#!/usr/bin/env python3
"""Replay smoke-testable backend requests from a HAR against a running gru server.

Examples:
    # Quick suite: deterministic endpoints only.
    python3 dev_tools/gene_smoketest.py --endpoint http://localhost:9001 --har gene.iobio.io.har

    # Replay only geneinfo requests.
    python3 dev_tools/gene_smoketest.py --endpoint http://localhost:9001 --prefix /geneinfo

    # Include slow/sampling-heavy endpoints with structural validation.
    python3 dev_tools/gene_smoketest.py --endpoint http://localhost:9001 --suite all

    # HTTPS endpoint with a self-signed certificate.
    python3 dev_tools/gene_smoketest.py --endpoint https://localhost:9001 --insecure
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import ssl
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable
from urllib.error import HTTPError
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen


@dataclass(frozen=True)
class ValidationResult:
    ok: bool
    details: list[str]


@dataclass(frozen=True)
class Validator:
    name: str
    slow: bool
    validate: Callable[[bytes, bytes], ValidationResult]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Replay backend requests from a HAR and verify deterministic endpoints exactly "
            "and sampling-heavy endpoints structurally."
        )
    )
    parser.add_argument(
        "har_path",
        nargs="?",
        default=None,
        help="HAR file to replay. Can also be supplied with --har.",
    )
    parser.add_argument(
        "--har",
        dest="har_option",
        default=None,
        help="HAR file to replay (overrides positional HAR; default: gene.iobio.io.har).",
    )
    parser.add_argument(
        "--endpoint",
        required=False,
        help=(
            "Base URL of the running server to test, e.g. http://localhost:9001. "
            "If it includes a path prefix such as /gru, HAR paths are replayed under it. "
            "Required unless --list is used."
        ),
    )
    parser.add_argument(
        "--suite",
        choices=("quick", "all"),
        default="quick",
        help=(
            "quick replays deterministic endpoints only. all also replays slow or "
            "sampling-heavy endpoints with structural validation (default: quick)."
        ),
    )
    parser.add_argument(
        "--prefix",
        action="append",
        default=None,
        help=(
            "Only replay HAR paths starting with this prefix. Can be supplied multiple "
            "times. By default all supported backend paths are considered."
        ),
    )
    parser.add_argument(
        "--origin",
        action="append",
        default=None,
        help=(
            "Only replay HAR requests from this origin host. Can be supplied multiple "
            "times. Use 'any' to disable origin filtering (default: backend.iobio.io)."
        ),
    )
    parser.add_argument(
        "--request-timeout",
        type=float,
        default=60.0,
        help="Seconds to wait for each replayed request (default: 60).",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification for https endpoints.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List selected/skipped HAR entries and validators without replaying requests.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print PASS lines in addition to failures.",
    )
    return parser.parse_args()


def decode_har_content(content: dict[str, Any]) -> bytes:
    text = content.get("text", "")
    if content.get("encoding") == "base64":
        return base64.b64decode(text)
    return text.encode("utf-8")


def request_body_from_har(request: dict[str, Any]) -> bytes | None:
    post_data = request.get("postData")
    if not post_data or "text" not in post_data:
        return None

    text = post_data["text"]
    if post_data.get("encoding") == "base64":
        return base64.b64decode(text)
    return text.encode("utf-8")


def target_url_from_har_url(har_url: str, endpoint: str) -> str:
    endpoint_parts = urlsplit(endpoint)
    if endpoint_parts.scheme not in {"http", "https"} or not endpoint_parts.netloc:
        raise ValueError(f"--endpoint must be an http/https URL with a host: {endpoint!r}")

    har_parts = urlsplit(har_url)
    base_path = endpoint_parts.path.rstrip("/")
    har_path = har_parts.path or "/"

    # Treat an endpoint path as a deployment prefix, but don't duplicate it if
    # the HAR path already includes that prefix.
    if base_path and har_path != base_path and not har_path.startswith(base_path + "/"):
        target_path = base_path + har_path
    else:
        target_path = har_path

    return urlunsplit((endpoint_parts.scheme, endpoint_parts.netloc, target_path, har_parts.query, ""))


def first_diff(expected: bytes, actual: bytes) -> tuple[int, bytes, bytes]:
    max_common = min(len(expected), len(actual))
    pos = 0
    while pos < max_common and expected[pos] == actual[pos]:
        pos += 1

    start = max(0, pos - 80)
    end = pos + 160
    return pos, expected[start:end], actual[start:end]


def ok(details: list[str] | None = None) -> ValidationResult:
    return ValidationResult(True, details or [])


def fail(details: list[str]) -> ValidationResult:
    return ValidationResult(False, details)


def text_lines(body: bytes) -> list[str]:
    return body.decode("utf-8", errors="replace").splitlines()


def validate_exact(expected: bytes, actual: bytes) -> ValidationResult:
    if actual == expected:
        return ok()

    pos, expected_slice, actual_slice = first_diff(expected, actual)
    return fail([
        f"body length actual/expected: {len(actual)}/{len(expected)}",
        f"first body diff byte: {pos}",
        f"expected slice: {expected_slice!r}",
        f"actual slice:   {actual_slice!r}",
    ])


def validate_gene_coverage(expected: bytes, actual: bytes) -> ValidationResult:
    errors: list[str] = []
    expected_lines = text_lines(expected)
    actual_lines = text_lines(actual)

    if not actual_lines:
        return fail(["empty geneCoverage response"])

    expected_header = expected_lines[0].split("\t") if expected_lines else []
    actual_header = actual_lines[0].split("\t")
    required_header = ["#id", "region", "min", "max", "q1", "median", "q3", "mean", "sd"]
    if actual_header != required_header:
        errors.append(f"unexpected geneCoverage header: {actual_header!r}")
    if expected_header and actual_header != expected_header:
        errors.append("geneCoverage header differs from HAR")

    expected_rows = [line.split("\t") for line in expected_lines[1:] if line]
    actual_rows = [line.split("\t") for line in actual_lines[1:] if line]
    if expected_rows and not actual_rows:
        errors.append(f"geneCoverage response has no rows; HAR had {len(expected_rows)}")

    for idx, row in enumerate(actual_rows):
        if len(row) != len(required_header):
            errors.append(f"geneCoverage row {idx + 1} has {len(row)} columns")
            continue

        if idx < len(expected_rows):
            exp = expected_rows[idx]
            # The requested bins/regions should be stable when they are present
            # in both responses, but the backend can legitimately return rows
            # where the HAR had only the header (and vice versa) depending on
            # sampling/remote data state.
            if len(exp) >= 2 and row[:2] != exp[:2]:
                errors.append(f"geneCoverage row {idx + 1} id/region differs: {row[:2]!r} != {exp[:2]!r}")

        for col_name, value in zip(required_header[2:], row[2:]):
            if value in {"", "NA", "NaN", "nan"}:
                continue
            try:
                number = float(value)
            except ValueError:
                errors.append(f"geneCoverage row {idx + 1} column {col_name} is not numeric: {value!r}")
                continue
            if not math.isfinite(number):
                errors.append(f"geneCoverage row {idx + 1} column {col_name} is not finite: {value!r}")

    return fail(errors) if errors else ok([f"rows={len(actual_rows)}"])


def validate_empty_gene2pheno_associations(expected: bytes, actual: bytes) -> ValidationResult:
    # Older Koa routing let /gene2pheno/associations/ fall through to the v2
    # /:gene route and returned []. The current fetch handler treats it as the
    # v3 associations route with an empty gene and returns an empty association
    # object. Both responses indicate "no gene/no associations".
    exact = validate_exact(expected, actual)
    if exact.ok:
        return exact

    allowed: list[Any] = [[], {"gene": "", "phenotypes": [], "disorders": []}]
    try:
        parsed = json.loads(actual.decode("utf-8"))
    except json.JSONDecodeError:
        return exact

    if parsed in allowed:
        return ok(["accepted empty gene2pheno associations response"])

    return exact


def parse_sections(body: bytes) -> dict[str, list[list[str]]]:
    sections: dict[str, list[list[str]]] = {}
    current: str | None = None
    for line in text_lines(body):
        if not line:
            continue
        if line.startswith("#"):
            current = line[1:]
            sections.setdefault(current, [])
        elif current is not None:
            sections[current].append(line.split("\t"))
    return sections


def validate_alignment_coverage(expected: bytes, actual: bytes) -> ValidationResult:
    errors: list[str] = []
    expected_sections = parse_sections(expected)
    actual_sections = parse_sections(actual)

    if not actual_sections:
        return fail(["empty or unsectioned alignmentCoverage response"])

    for section in ("specific_points", "reduced_points"):
        if section not in actual_sections:
            errors.append(f"missing alignmentCoverage section #{section}")

    for section, rows in actual_sections.items():
        if section not in expected_sections:
            errors.append(f"unexpected alignmentCoverage section #{section}")
        for row_idx, row in enumerate(rows):
            if len(row) != 2:
                errors.append(f"alignmentCoverage #{section} row {row_idx + 1} has {len(row)} columns")
                continue
            for value in row:
                try:
                    int(value)
                except ValueError:
                    errors.append(f"alignmentCoverage #{section} row {row_idx + 1} non-integer value: {value!r}")

    # The exact reduced sampling can vary. Keep this intentionally structural,
    # but require the explicitly requested specific points to be stable.
    exp_specific = expected_sections.get("specific_points", [])
    act_specific = actual_sections.get("specific_points", [])
    if exp_specific and act_specific != exp_specific:
        errors.append("alignmentCoverage #specific_points differs from HAR")

    exp_reduced = len(expected_sections.get("reduced_points", []))
    act_reduced = len(actual_sections.get("reduced_points", []))
    if exp_reduced and act_reduced == 0:
        errors.append("alignmentCoverage #reduced_points is empty")

    return fail(errors) if errors else ok([
        f"specific_points={len(act_specific)}",
        f"reduced_points={act_reduced}",
    ])


def vcf_stats(body: bytes) -> tuple[list[str], str | None, list[str]]:
    lines = text_lines(body)
    meta = [line for line in lines if line.startswith("##")]
    chrom_header = next((line for line in lines if line.startswith("#CHROM")), None)
    records = [line for line in lines if line and not line.startswith("#")]
    return meta, chrom_header, records


def validate_vcf_like(expected: bytes, actual: bytes) -> ValidationResult:
    errors: list[str] = []
    exp_meta, exp_chrom_header, exp_records = vcf_stats(expected)
    act_meta, act_chrom_header, act_records = vcf_stats(actual)

    if not act_meta or not act_meta[0].startswith("##fileformat=VCF"):
        errors.append("VCF response missing ##fileformat header")
    if act_chrom_header is None:
        errors.append("VCF response missing #CHROM header")

    if exp_chrom_header and act_chrom_header:
        exp_samples = exp_chrom_header.split("\t")[9:]
        act_samples = act_chrom_header.split("\t")[9:]
        if exp_samples != act_samples:
            errors.append(f"VCF samples actual/expected: {act_samples!r}/{exp_samples!r}")

    if exp_records and not act_records:
        errors.append("VCF response has no variant records")

    for idx, record in enumerate(act_records[:100]):
        cols = record.split("\t")
        if len(cols) < 8:
            errors.append(f"VCF record {idx + 1} has {len(cols)} columns")
            continue
        try:
            int(cols[1])
        except ValueError:
            errors.append(f"VCF record {idx + 1} POS is not an integer: {cols[1]!r}")

    # Don't require identical counts; annotate/freebayes outputs can vary with
    # sampling, tool versions, and remote reference data. Do require the current
    # result to be in the same broad ballpark when the HAR had records.
    if exp_records:
        min_records = max(1, int(len(exp_records) * 0.05))
        if len(act_records) < min_records:
            errors.append(f"VCF record count too small actual/min/expected: {len(act_records)}/{min_records}/{len(exp_records)}")

    return fail(errors) if errors else ok([
        f"meta={len(act_meta)}",
        f"records={len(act_records)}",
    ])


def validator_for_path(path: str) -> Validator | None:
    if path == "/gene2pheno/associations/":
        return Validator("gene2pheno-empty-compatible", False, validate_empty_gene2pheno_associations)

    # Deterministic/cacheable JSON and header-style endpoints. These are exact
    # byte comparisons against the HAR body.
    if (
        path.startswith("/geneinfo")
        or path.startswith("/gene2pheno")
        or path.startswith("/hpo/hot/lookup")
        or path == "/getChromosomes"
        or path == "/variantHeader"
        or path == "/alignmentHeader"
    ):
        return Validator("exact", False, validate_exact)

    # Sampling/tool-heavy endpoints. Validate status and output shape, not exact
    # bytes.
    if path == "/geneCoverage":
        return Validator("geneCoverage-structural", True, validate_gene_coverage)
    if path == "/alignmentCoverage":
        return Validator("alignmentCoverage-structural", True, validate_alignment_coverage)
    if path in {"/annotateVariantsV3", "/freebayesJointCallV3"}:
        return Validator("vcf-structural", True, validate_vcf_like)

    return None


def replay_entry(
    entry: dict[str, Any],
    endpoint: str,
    timeout: float,
    ssl_context: ssl.SSLContext | None,
) -> tuple[int, bytes]:
    request = entry["request"]
    method = request.get("method", "GET")
    target_url = target_url_from_har_url(request["url"], endpoint)
    body = request_body_from_har(request)

    headers: dict[str, str] = {}
    for header in request.get("headers", []):
        name = header.get("name")
        value = header.get("value")
        if not name or value is None:
            continue

        # Let urllib set local transport details. Also avoid compressed response
        # bodies so byte comparisons are against the real response payload.
        if name.startswith(":") or name.lower() in {
            "host",
            "connection",
            "content-length",
            "accept-encoding",
        }:
            continue

        headers[name] = value

    req = Request(target_url, data=body, headers=headers, method=method)
    try:
        with urlopen(req, timeout=timeout, context=ssl_context) as response:
            return response.status, response.read()
    except HTTPError as err:
        return err.code, err.read()


def entry_label(idx: int, entry: dict[str, Any]) -> str:
    request = entry["request"]
    har_url_parts = urlsplit(request["url"])
    label = f"#{idx} {request.get('method', 'GET')} {har_url_parts.path}"
    if har_url_parts.query:
        label += f"?{har_url_parts.query}"
    return label


def origin_matches(netloc: str, origins: list[str]) -> bool:
    if "any" in origins:
        return True
    return netloc in origins


def prefix_matches(path: str, prefixes: list[str] | None) -> bool:
    if not prefixes:
        return True
    return any(path.startswith(prefix) for prefix in prefixes)


def classify_entries(har: dict[str, Any], args: argparse.Namespace) -> tuple[list[tuple[int, dict[str, Any], Validator]], Counter[str]]:
    origins = args.origin or ["backend.iobio.io"]
    selected: list[tuple[int, dict[str, Any], Validator]] = []
    skipped: Counter[str] = Counter()

    for idx, entry in enumerate(har.get("log", {}).get("entries", [])):
        request_url = entry["request"]["url"]
        url_parts = urlsplit(request_url)
        if not origin_matches(url_parts.netloc, origins):
            skipped["origin"] += 1
            continue
        if not prefix_matches(url_parts.path, args.prefix):
            skipped["prefix"] += 1
            continue

        validator = validator_for_path(url_parts.path)
        if validator is None:
            skipped["unsupported"] += 1
            continue
        if validator.slow and args.suite != "all":
            skipped["slow"] += 1
            continue

        selected.append((idx, entry, validator))

    return selected, skipped


def main() -> int:
    args = parse_args()

    har_arg = args.har_option or args.har_path or "gene.iobio.io.har"
    har_path = Path(har_arg).resolve()
    if not har_path.exists():
        print(f"HAR not found: {har_path}", file=sys.stderr)
        return 2

    if not args.endpoint and not args.list:
        print("--endpoint is required unless --list is used", file=sys.stderr)
        return 2

    endpoint = args.endpoint.rstrip("/") if args.endpoint else "http://example.test"
    ssl_context = ssl._create_unverified_context() if args.insecure else None

    try:
        # Validate endpoint before doing any replay work.
        target_url_from_har_url("https://example.test/geneinfo/ping", endpoint)
    except ValueError as err:
        print(str(err), file=sys.stderr)
        return 2

    har = json.loads(har_path.read_text(encoding="utf-8"))
    selected, skipped = classify_entries(har, args)

    if not selected:
        print("No HAR entries selected for replay", file=sys.stderr)
        if skipped:
            print("Skipped:", ", ".join(f"{k}={v}" for k, v in sorted(skipped.items())), file=sys.stderr)
        return 2

    validator_counts = Counter(validator.name for _, _, validator in selected)
    print(
        "Selected:",
        f"total={len(selected)}",
        "validators=" + ", ".join(f"{k}={v}" for k, v in sorted(validator_counts.items())),
    )
    if skipped:
        print("Skipped:", ", ".join(f"{k}={v}" for k, v in sorted(skipped.items())))

    if args.list:
        by_validator: dict[str, list[str]] = defaultdict(list)
        for idx, entry, validator in selected:
            by_validator[validator.name].append(entry_label(idx, entry))
        for validator_name, labels in sorted(by_validator.items()):
            print(f"\n[{validator_name}] {len(labels)} entries")
            for label in labels:
                print("  ", label)
        return 0

    failures = 0
    for idx, entry, validator in selected:
        response = entry["response"]
        label = entry_label(idx, entry)

        expected_status = int(response.get("status", 0))
        expected_body = decode_har_content(response.get("content", {}))

        try:
            actual_status, actual_body = replay_entry(
                entry,
                endpoint,
                args.request_timeout,
                ssl_context,
            )
        except Exception as err:
            failures += 1
            print(f"FAIL {label} [{validator.name}]")
            print(f"  request error: {err!r}")
            continue

        errors: list[str] = []
        if actual_status != expected_status:
            errors.append(f"status actual/expected: {actual_status}/{expected_status}")
        if b"GRU_ERROR_SENTINEL" in actual_body:
            errors.append("response contains GRU_ERROR_SENTINEL")

        validation = validator.validate(expected_body, actual_body)
        if not validation.ok:
            errors.extend(validation.details)

        if errors:
            failures += 1
            print(f"FAIL {label} [{validator.name}]")
            for error in errors:
                print(f"  {error}")
            continue

        if args.verbose:
            suffix = ""
            if validation.details:
                suffix = " " + " ".join(validation.details)
            print(f"PASS {label} [{validator.name}] status={actual_status} len={len(actual_body)}{suffix}")

    passed = len(selected) - failures
    print(f"Summary: total={len(selected)} passed={passed} failed={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
