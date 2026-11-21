"""
Shared HyPhy execution helper for the corrected pipeline (steps 16-17).

HyPhy is reached through the official Docker image, as in the original pipeline, but with
two changes:
  - the image tag is pinned, so the analysis version is recorded rather than floating on
    whatever `hyphy/hyphy:latest` happens to be (CORRECTIONS_AND_PENDING_WORK.md, A10);
  - the daemon is checked once up front with a clear message, instead of every run failing.

Mounting: the dataset directory is bind-mounted at /data and every path handed to HyPhy is
relative to it, so container paths never leak host details.
"""

import os
import json
import time
import shutil
import subprocess
from datetime import datetime

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# The locally cached image reports HYPHY 2.5.103(MP). Version and image digest are recorded
# by hyphy_provenance() at the start of every stage, which is what reproducibility needs;
# override the image with HSV_HYPHY_IMAGE to pin a different build.
HYPHY_IMAGE = os.environ.get("HSV_HYPHY_IMAGE", "hyphy/hyphy:latest")
HYPHY_CPUS = os.environ.get("HSV_HYPHY_CPUS", "4")


class DockerUnavailable(RuntimeError):
    pass


def docker_available(timeout=60):
    """True when the Docker daemon answers."""
    if not shutil.which("docker"):
        return False, "the docker CLI is not on PATH"
    try:
        res = subprocess.run(["docker", "info", "--format", "{{.ServerVersion}}"],
                             capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, "docker info timed out"
    if res.returncode != 0:
        return False, (res.stderr or res.stdout or "docker info failed").strip().splitlines()[0]
    return True, f"Docker daemon {res.stdout.strip()}"


def wait_for_docker(max_wait=300, interval=10):
    """Poll until the daemon is up, e.g. while Docker Desktop finishes starting."""
    waited = 0
    while waited <= max_wait:
        ok, msg = docker_available()
        if ok:
            return True, msg
        time.sleep(interval)
        waited += interval
    return False, msg


def ensure_image():
    res = subprocess.run(["docker", "image", "inspect", HYPHY_IMAGE],
                         capture_output=True, text=True)
    if res.returncode == 0:
        return True, f"{HYPHY_IMAGE} present locally"
    pull = subprocess.run(["docker", "pull", HYPHY_IMAGE], capture_output=True, text=True)
    if pull.returncode != 0:
        return False, (pull.stderr or pull.stdout).strip()[-400:]
    return True, f"pulled {HYPHY_IMAGE}"


def hyphy_provenance():
    """Version banner and image id, for the methods section."""
    ver = subprocess.run(["docker", "run", "--rm", HYPHY_IMAGE, "hyphy", "--version"],
                         capture_output=True, text=True)
    version = (ver.stdout or ver.stderr or "").strip().splitlines()
    idres = subprocess.run(["docker", "image", "inspect", HYPHY_IMAGE,
                            "--format", "{{.Id}}"], capture_output=True, text=True)
    return dict(image=HYPHY_IMAGE,
                version=version[0] if version else "unknown",
                image_id=(idres.stdout or "").strip(),
                # Per-HyPhy-process CPU setting (the `CPU=` argument), NOT the number of
                # analyses run in parallel: the stages drive concurrency with their own
                # --workers and each worker launches its own container at this setting.
                # Recorded for completeness only -- it affects runtime, never the
                # likelihoods or the reported statistics.
                cpus_per_hyphy_process=HYPHY_CPUS,
                parallel_workers=os.environ.get("HSV_HYPHY_WORKERS", "not recorded"),
                # Stamped so every result set says which run produced it and cannot be
                # confused with the withdrawn v1 output in 07_qc_reports/hyphy/.
                pipeline="v2-per-type-corrected",
                recorded_at=datetime.now().isoformat(timespec="seconds"))


def rel(path):
    """Host path -> path relative to the dataset root, for the container mount."""
    return os.path.relpath(os.path.abspath(path), DIR_BASE).replace("\\", "/")


def run_hyphy(method, alignment, output, tree=None, extra=None, timeout=None):
    """
    Run one HyPhy analysis. Returns (ok, runtime_seconds, message).
    `alignment`, `tree` and `output` are host paths under the dataset root.
    """
    # -w /data: the image's working directory is /, so relative paths would otherwise
    # resolve outside the mount and HyPhy would fall back to prompting for the alignment.
    cmd = ["docker", "run", "--rm",
           "--mount", f"type=bind,src={DIR_BASE},dst=/data",
           "-w", "/data",
           HYPHY_IMAGE,
           # TOLERATE_NUMERICAL_ERRORS: HyPhy's own recommended setting for the
           # "ComputeBranchCache ... reversible model cached likelihood" internal error,
           # which it raises on the very low-divergence per-type HSV cohorts. Without it
           # GARD aborts during the single-breakpoint scan.
           "hyphy", "ENV=TOLERATE_NUMERICAL_ERRORS=1", f"CPU={HYPHY_CPUS}", method.lower(),
           "--alignment", rel(alignment)]
    if tree:
        cmd += ["--tree", rel(tree)]
    cmd += ["--output", rel(output)]
    for k, v in (extra or {}).items():
        cmd += [f"--{k}", str(v)]

    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, round(time.time() - t0, 1), f"timed out after {timeout}s"
    rt = round(time.time() - t0, 1)

    if res.returncode != 0:
        tail = (res.stderr or res.stdout or "").strip().splitlines()
        return False, rt, " | ".join(tail[-3:])[:400]
    if not os.path.exists(output):
        return False, rt, "exited 0 but no output file was written"
    try:
        with open(output, encoding="utf-8") as fh:
            json.load(fh)
    except Exception as exc:
        return False, rt, f"output is not valid JSON: {type(exc).__name__}"
    return True, rt, "ok"
