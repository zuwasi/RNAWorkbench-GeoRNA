#!/usr/bin/env python3
"""Report the features available to a specific Parasoft cpptestcli."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path


FEATURES = (
    "C++test",
    "Static Analysis",
    "Flow Analysis",
    "Unit Test",
    "Coverage",
    "Automation",
    "DTP Publish",
    "Desktop Command Line",
    "Embedded Support",
    "Coding Standards",
    "Requirements Traceability",
    "RuleWizard",
    "Runtime Error Detection",
    "AUTOSAR Rules",
    "CWE Rules",
    "DISA STIG Rules",
    "HIC++ Rules",
    "JSF Rules",
    "MISRA Rules",
    "MISRA C 2012 Rules",
    "MISRA C 2025 Rules",
    "MISRA C++ 2023 Rules",
    "OWASP Rules",
    "Security Rules",
    "SEI CERT C Rules",
    "SEI CERT C++ Rules",
    "Code Dependency Export",
    "LLM Integration",
)

UNKNOWN_RE = re.compile(
    r"Unknown license feature names in the -check-license option:\s*(.+)",
    re.IGNORECASE,
)
MISSING_RE = re.compile(r"License features not set:\s*(.+)", re.IGNORECASE)
ALL_SET_RE = re.compile(r"All license features are set:", re.IGNORECASE)
FATAL_MARKERS = (
    "Unable to connect to License Server",
    "License Server host is unreachable",
    "There is no license for this product",
    "There is no valid license",
    "License not set - aborting",
    "No more tokens available",
    "Request authentication failed",
    "Request denied",
    "Request time out",
    "license has expired",
    "EULA",
)


@dataclass(frozen=True)
class ProbeResult:
    licensed: tuple[str, ...]
    not_licensed: tuple[str, ...]
    not_recognized: tuple[str, ...]
    version: str


class ProbeError(RuntimeError):
    pass


def _split_features(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


def _display_name(feature: str) -> str:
    return "C++Test" if feature.casefold() == "c++test" else f"C++Test {feature}"


def classify_license_mode(features: tuple[str, ...]) -> str:
    has_automation = "Automation" in features
    has_desktop_cli = "Desktop Command Line" in features
    if has_automation and has_desktop_cli:
        return "Both Automation and Desktop CLI"
    if has_automation:
        return "Automation edition"
    if has_desktop_cli:
        return "Desktop CLI only"
    return "Neither Automation nor Desktop CLI"


def _run(
    executable: Path,
    features: tuple[str, ...],
    workspace: Path,
    settings: Path | None,
    timeout: float,
) -> subprocess.CompletedProcess[str]:
    command = [
        str(executable),
        "-workspace",
        str(workspace),
    ]
    if settings is not None:
        command.extend(("-settings", str(settings)))
    command.extend(("-check-license", ",".join(features)))
    try:
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            errors="replace",
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ProbeError(f"License check timed out after {timeout:g} seconds") from exc
    except OSError as exc:
        raise ProbeError(f"Could not run {executable}: {exc}") from exc


def probe(
    executable: Path,
    settings: Path | None = None,
    timeout: float = 60.0,
) -> ProbeResult:
    if not executable.is_file():
        raise ProbeError(f"cpptestcli was not found: {executable}")
    if settings is not None and not settings.is_file():
        raise ProbeError(f"Settings file was not found: {settings}")

    supported = list(FEATURES)
    unrecognized: list[str] = []
    with tempfile.TemporaryDirectory(prefix="cpptest-license-") as workspace_name:
        workspace = Path(workspace_name)
        while supported:
            completed = _run(executable, tuple(supported), workspace, settings, timeout)
            output = f"{completed.stdout}\n{completed.stderr}"

            unknown_match = UNKNOWN_RE.search(output)
            if unknown_match:
                unknown = _split_features(unknown_match.group(1))
                removed = [feature for feature in supported if feature in unknown]
                if not removed:
                    raise ProbeError(unknown_match.group(0).strip())
                unrecognized.extend(removed)
                supported = [feature for feature in supported if feature not in unknown]
                continue

            fatal = next(
                (
                    line.strip()
                    for line in output.splitlines()
                    if any(
                        marker.casefold() in line.casefold()
                        for marker in FATAL_MARKERS
                    )
                ),
                None,
            )
            if fatal:
                raise ProbeError(fatal)

            version = next(
                (line.strip() for line in output.splitlines() if line.startswith("Parasoft C/C++test")),
                "Unknown",
            )
            missing_match = MISSING_RE.search(output)
            if missing_match:
                missing_names = _split_features(missing_match.group(1))
                not_licensed = tuple(feature for feature in supported if feature in missing_names)
                licensed = tuple(feature for feature in supported if feature not in missing_names)
                return ProbeResult(licensed, not_licensed, tuple(unrecognized), version)
            if completed.returncode == 0 and ALL_SET_RE.search(output):
                return ProbeResult(tuple(supported), (), tuple(unrecognized), version)

            detail = next(
                (
                    line.strip()
                    for line in reversed(output.splitlines())
                    if line.strip()
                ),
                "no diagnostic output",
            )
            raise ProbeError(
                f"cpptestcli returned exit code {completed.returncode}: {detail}"
            )

    raise ProbeError("This cpptestcli does not recognize any known C/C++test features")


def _print_wrapped(label: str, features: tuple[str, ...]) -> None:
    value = ", ".join(_display_name(feature) for feature in features) or "None"
    prefix = f"{label:<20}: "
    print(
        textwrap.fill(
            value,
            width=100,
            initial_indent=prefix,
            subsequent_indent=" " * len(prefix),
        )
    )


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query a specific cpptestcli license without analyzing a project."
    )
    parser.add_argument("cpptestcli", type=Path, help="Path to cpptestcli or cpptestcli.exe")
    parser.add_argument("--settings", type=Path, help="Optional Parasoft .properties file")
    parser.add_argument("--timeout", type=float, default=60.0, help="Timeout in seconds (default: 60)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        result = probe(args.cpptestcli.resolve(), args.settings, args.timeout)
    except ProbeError as exc:
        print(f"License query failed: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(
            json.dumps(
                {
                    "cpptestcli": str(args.cpptestcli.resolve()),
                    "version": result.version,
                    "license_mode": classify_license_mode(result.licensed),
                    "features": [_display_name(item) for item in result.licensed],
                    "not_licensed": [_display_name(item) for item in result.not_licensed],
                    "not_recognized": [
                        _display_name(item) for item in result.not_recognized
                    ],
                },
                indent=2,
            )
        )
    else:
        print(f"C/C++test           : {result.version}")
        print(f"License mode        : {classify_license_mode(result.licensed)}")
        _print_wrapped("Features", result.licensed)
        _print_wrapped("Not licensed", result.not_licensed)
        _print_wrapped("Not recognized", result.not_recognized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
