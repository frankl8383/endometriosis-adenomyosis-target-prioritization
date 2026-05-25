#!/usr/bin/env python3
"""Report status for detached long download jobs."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

from start_long_downloads import JOBS, LOG_DIR, PID_FILE, PROJECT_ROOT


def is_running(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def describe_path(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"exists": False}
    stat = path.stat()
    allocated_bytes = getattr(stat, "st_blocks", 0) * 512
    return {
        "exists": True,
        "logical_bytes": stat.st_size,
        "allocated_bytes": allocated_bytes,
        "modified_at": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
    }


def tail_text(path: Path, lines: int = 12) -> list[str]:
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    return text.splitlines()[-lines:]


def load_records() -> list[dict[str, object]]:
    if not PID_FILE.exists():
        return []
    try:
        return json.loads(PID_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []


def main() -> int:
    records = load_records()
    by_name = {str(record.get("name")): record for record in records}
    print(f"PID_FILE\t{PID_FILE}\t{'present' if PID_FILE.exists() else 'missing'}")

    for job in JOBS:
        name = str(job["name"])
        target = PROJECT_ROOT / str(job["target"])
        partial = target.with_suffix(target.suffix + ".partial")
        aria2_state = Path(str(partial) + ".aria2")
        record = by_name.get(name, {})
        pid = record.get("pid")
        running = is_running(int(pid)) if pid else False

        print(f"\nJOB\t{name}")
        print(f"PID\t{pid if pid else 'NA'}\t{'running' if running else 'not_running'}")
        print(f"TARGET\t{json.dumps(describe_path(target), ensure_ascii=False)}")
        print(f"PARTIAL\t{json.dumps(describe_path(partial), ensure_ascii=False)}")
        print(f"ARIA2_STATE\t{json.dumps(describe_path(aria2_state), ensure_ascii=False)}")

        stdout_log = Path(str(record.get("stdout_log", ""))) if record.get("stdout_log") else None
        stderr_log = Path(str(record.get("stderr_log", ""))) if record.get("stderr_log") else None
        if stdout_log:
            print(f"STDOUT_LOG\t{stdout_log}")
            for line in tail_text(stdout_log):
                print(f"  stdout> {line}")
        if stderr_log:
            print(f"STDERR_LOG\t{stderr_log}")
            for line in tail_text(stderr_log):
                print(f"  stderr> {line}")

    if not LOG_DIR.exists():
        print(f"\nLOG_DIR\t{LOG_DIR}\tmissing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
