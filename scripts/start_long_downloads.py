#!/usr/bin/env python3
"""Start resumable long downloads as detached background jobs."""

from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = PROJECT_ROOT / "logs" / "long_downloads"
PID_FILE = LOG_DIR / "active_downloads.json"

JOBS = [
    {
        "name": "zenodo_17078290_h5ad",
        "cmd": ["python3", "scripts/download_phase1_files.py", "phase1_optional", "aria2_relaxed"],
        "target": "data/raw_downloads/Zenodo_17078290__celllable.diff_PRO.h5ad",
    },
    {
        "name": "gse179640_raw_tar",
        "cmd": ["python3", "scripts/download_phase1_files.py", "inspect_first", "aria2"],
        "target": "data/raw_downloads/GSE179640__GSE179640_RAW.tar",
    },
]


def is_running(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def load_existing() -> list[dict[str, object]]:
    if not PID_FILE.exists():
        return []
    try:
        return json.loads(PID_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []


def main() -> int:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    existing = load_existing()
    active_by_name = {
        str(item.get("name")): item
        for item in existing
        if item.get("pid") and is_running(int(item["pid"]))
    }

    records: list[dict[str, object]] = []
    for job in JOBS:
        target = PROJECT_ROOT / str(job["target"])
        partial = target.with_suffix(target.suffix + ".partial")
        aria2_state = Path(str(partial) + ".aria2")
        if target.exists() and target.stat().st_size > 0 and not aria2_state.exists():
            record = {
                "name": job["name"],
                "pid": None,
                "cmd": job["cmd"],
                "target": str(target),
                "stdout_log": "",
                "stderr_log": "",
                "started_at": "",
                "status": "target_exists_no_active_download",
            }
            records.append(record)
            print(f"SKIP_COMPLETED {job['name']} target={target}")
            continue

        if job["name"] in active_by_name:
            records.append(active_by_name[job["name"]])
            print(f"ALREADY_RUNNING {job['name']} pid={active_by_name[job['name']]['pid']}")
            continue

        timestamp = time.strftime("%Y%m%d_%H%M%S")
        stdout_path = LOG_DIR / f"{job['name']}.{timestamp}.log"
        stderr_path = LOG_DIR / f"{job['name']}.{timestamp}.err"
        stdout_handle = stdout_path.open("ab")
        stderr_handle = stderr_path.open("ab")
        proc = subprocess.Popen(
            job["cmd"],
            cwd=PROJECT_ROOT,
            stdout=stdout_handle,
            stderr=stderr_handle,
            stdin=subprocess.DEVNULL,
            start_new_session=True,
        )
        stdout_handle.close()
        stderr_handle.close()
        record = {
            "name": job["name"],
            "pid": proc.pid,
            "cmd": job["cmd"],
            "target": str(PROJECT_ROOT / job["target"]),
            "stdout_log": str(stdout_path),
            "stderr_log": str(stderr_path),
            "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        records.append(record)
        print(f"STARTED {job['name']} pid={proc.pid}")

    PID_FILE.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"PID_FILE {PID_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
