#!/usr/bin/env python3
"""Replay /geneinfo requests from a HAR against a running gru server.

Examples:
    python3 dev_tools/gene_smoketest.py --endpoint http://localhost:9001
    python3 dev_tools/gene_smoketest.py --endpoint https://example.org --har gene.iobio.io.har
    python3 dev_tools/gene_smoketest.py gene.iobio.io.har --endpoint http://localhost:9001
"""

from __future__ import annotations

import argparse
import base64
import json
import ssl
import sys
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay /geneinfo requests from a HAR and compare status/body exactly."
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
        help="HAR file to replay (overrides positional HAR).",
    )
    parser.add_argument(
        "--endpoint",
        required=True,
        help=(
            "Base URL of the running server to test, e.g. http://localhost:9001. "
            "If it includes a path prefix such as /gru, the HAR path is replayed under it."
        ),
    )
    parser.add_argument(
        "--prefix",
        default="/geneinfo",
        help="Only replay requests whose HAR path starts with this prefix (default: /geneinfo)",
    )
    parser.add_argument(
        "--request-timeout",
        type=float,
        default=60.0,
        help="Seconds to wait for each replayed request (default: 60)",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification for https endpoints.",
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


def main() -> int:
    args = parse_args()

    har_arg = args.har_option or args.har_path or "gene.iobio.io.har"
    har_path = Path(har_arg).resolve()
    if not har_path.exists():
        print(f"HAR not found: {har_path}", file=sys.stderr)
        return 2

    endpoint = args.endpoint.rstrip("/")
    ssl_context = ssl._create_unverified_context() if args.insecure else None

    try:
        # Validate endpoint before doing any replay work.
        target_url_from_har_url("https://example.test/geneinfo/ping", endpoint)
    except ValueError as err:
        print(str(err), file=sys.stderr)
        return 2

    har = json.loads(har_path.read_text(encoding="utf-8"))
    entries = []
    for idx, entry in enumerate(har.get("log", {}).get("entries", [])):
        path = urlsplit(entry["request"]["url"]).path
        if path.startswith(args.prefix):
            entries.append((idx, entry))

    if not entries:
        print(f"No HAR entries found with path prefix {args.prefix!r}", file=sys.stderr)
        return 2

    failures = 0
    for idx, entry in entries:
        request = entry["request"]
        response = entry["response"]
        har_url_parts = urlsplit(request["url"])
        label = f"#{idx} {request.get('method', 'GET')} {har_url_parts.path}"
        if har_url_parts.query:
            label += f"?{har_url_parts.query}"

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
            print(f"FAIL {label}")
            print(f"  request error: {err!r}")
            continue

        status_ok = actual_status == expected_status
        body_ok = actual_body == expected_body
        if status_ok and body_ok:
            if args.verbose:
                print(f"PASS {label} status={actual_status} len={len(actual_body)}")
            continue

        failures += 1
        print(f"FAIL {label}")
        print(
            f"  status actual/expected: {actual_status}/{expected_status} "
            f"({'OK' if status_ok else 'DIFF'})"
        )
        print(
            f"  body length actual/expected: {len(actual_body)}/{len(expected_body)} "
            f"({'OK' if body_ok else 'DIFF'})"
        )
        if not body_ok:
            pos, expected_slice, actual_slice = first_diff(expected_body, actual_body)
            print(f"  first body diff byte: {pos}")
            print(f"  expected slice: {expected_slice!r}")
            print(f"  actual slice:   {actual_slice!r}")

    passed = len(entries) - failures
    print(f"Summary: total={len(entries)} passed={passed} failed={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
