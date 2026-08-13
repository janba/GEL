#!/usr/bin/env python3
"""Build and run every GEL demo under src/demo.

For each immediate subdirectory:

* If it contains CMakeLists.txt, create a ``build`` directory, configure
  and compile with CMake, then execute the resulting binary (or binaries).
* If it contains Python scripts, run each of them with the ``python3``
  (or ``python``) found on PATH, using the pip-installed ``pygel3d``.

Interactive OpenGL / GLUT demos are stopped after a timeout so the harness
can continue.  A timeout is reported separately from a crash.

Usage:
    python run_demos.py
    python run_demos.py --only MeshDistance kDTree
    python run_demos.py --timeout 60 --jobs 8
    python run_demos.py --list
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

DEMO_ROOT = Path(__file__).resolve().parent
REPO_ROOT = DEMO_ROOT.parents[1]
GEL_BUILD_DIR = REPO_ROOT / "build"
GEL_INSTALL_DIR = Path.home() / ".local"
IN_TREE_PYGEL = REPO_ROOT / "src" / "PyGEL"

ADD_EXECUTABLE_RE = re.compile(r"add_executable\s*\(\s*([A-Za-z0-9_.-]+)")

SKIP_PYTHON_NAMES = {"run_demos.py"}
SKIP_DIR_NAMES = {"build", "__pycache__", ".ipynb_checkpoints"}


# ---------------------------------------------------------------------------
# Result bookkeeping
# ---------------------------------------------------------------------------

@dataclass
class StepResult:
    demo: str
    kind: str
    name: str
    status: str
    detail: str = ""
    seconds: float = 0.0


@dataclass
class Harness:
    results: list[StepResult] = field(default_factory=list)

    def add(self, result: StepResult) -> None:
        self.results.append(result)
        mark = {
            "ok": "PASS",
            "timeout": "TIME",
            "fail": "FAIL",
            "skip": "SKIP",
        }.get(result.status, result.status.upper())
        timing = f" ({result.seconds:.1f}s)" if result.seconds else ""
        extra = f" — {result.detail}" if result.detail else ""
        print(f"  [{mark}] {result.kind}: {result.name}{timing}{extra}")

    def summary(self) -> int:
        counts = {key: 0 for key in ("ok", "timeout", "fail", "skip")}
        for item in self.results:
            counts[item.status] = counts.get(item.status, 0) + 1

        width = 70
        print()
        print("=" * width)
        print("DEMO HARNESS SUMMARY")
        print("=" * width)
        print(f"  passed:   {counts['ok']}")
        print(f"  timeout:  {counts['timeout']}  (interactive demos stopped by --timeout)")
        print(f"  failed:   {counts['fail']}")
        print(f"  skipped:  {counts['skip']}")

        failed = [item for item in self.results if item.status == "fail"]
        if failed:
            print()
            print("Failures:")
            for item in failed:
                print(f"  - {item.demo} / {item.kind} / {item.name}: {item.detail}")

        print("=" * width)
        return 1 if failed else 0


# ---------------------------------------------------------------------------
# Process helpers
# ---------------------------------------------------------------------------

def kill_process_tree(proc: subprocess.Popen) -> None:
    """Stop a demo and any children it spawned (GLUT windows, etc.)."""
    if proc.poll() is not None:
        return
    try:
        if os.name == "nt":
            proc.kill()
        else:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        try:
            proc.kill()
        except OSError:
            pass


def run_cmd(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
    timeout: float | None = None,
    stream: bool = False,
) -> tuple[int | None, str]:
    """Run ``cmd``.

    Returns ``(returncode, output)``.  ``returncode`` is ``None`` on timeout.
    """
    popen_kwargs: dict = {
        "cwd": str(cwd),
        "env": env,
        "stdout": None if stream else subprocess.PIPE,
        "stderr": None if stream else subprocess.STDOUT,
        "text": True,
    }
    if os.name == "nt":
        popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        popen_kwargs["start_new_session"] = True

    proc = subprocess.Popen(cmd, **popen_kwargs)
    try:
        stdout, _ = proc.communicate(timeout=timeout)
        return proc.returncode, stdout or ""
    except subprocess.TimeoutExpired:
        kill_process_tree(proc)
        stdout, _ = proc.communicate()
        return None, stdout or ""


def indent_output(text: str, *, limit: int = 40) -> None:
    lines = text.splitlines()
    if not lines:
        return
    if len(lines) > limit:
        lines = ["…"] + lines[-limit:]
    for line in lines:
        print(f"        {line}")


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def demo_dirs(only: list[str] | None = None) -> list[Path]:
    names = set(only) if only else None
    dirs = []
    for entry in sorted(DEMO_ROOT.iterdir(), key=lambda p: p.name.lower()):
        if not entry.is_dir() or entry.name.startswith("."):
            continue
        if entry.name in SKIP_DIR_NAMES:
            continue
        if names is not None and entry.name not in names:
            continue
        dirs.append(entry)
    return dirs


def cmake_executables(cmake_lists: Path) -> list[str]:
    text = cmake_lists.read_text(encoding="utf-8", errors="replace")
    return ADD_EXECUTABLE_RE.findall(text)


def python_scripts(demo_dir: Path) -> list[Path]:
    scripts = []
    for path in sorted(demo_dir.glob("*.py")):
        if path.name in SKIP_PYTHON_NAMES:
            continue
        scripts.append(path)
    return scripts


def find_built_binary(build_dir: Path, name: str) -> Path | None:
    candidates = [
        build_dir / name,
        build_dir / f"{name}.exe",
        build_dir / "Release" / name,
        build_dir / "Release" / f"{name}.exe",
        build_dir / "Debug" / name,
        build_dir / "Debug" / f"{name}.exe",
        build_dir / "RelWithDebInfo" / name,
        build_dir / "RelWithDebInfo" / f"{name}.exe",
    ]
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return path

    # Fallback: search the build tree, ignoring CMake compiler-id binaries.
    matches = []
    for path in build_dir.rglob(name):
        if not path.is_file() or not os.access(path, os.X_OK):
            continue
        if "CMakeFiles" in path.parts:
            continue
        matches.append(path)
    for path in build_dir.rglob(f"{name}.exe"):
        if path.is_file():
            matches.append(path)
    return matches[0] if matches else None


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

def gel_package_config() -> Path | None:
    """Return GELConfig.cmake from the local install prefix, if present."""
    for libdir in ("lib", "lib64"):
        candidate = GEL_INSTALL_DIR / libdir / "cmake" / "GEL" / "GELConfig.cmake"
        if candidate.is_file():
            return candidate
    return None


def path_python() -> str | None:
    """Return the first ``python3`` or ``python`` on PATH."""
    for name in ("python3", "python"):
        found = shutil.which(name)
        if found:
            return found
    return None


def demo_env() -> dict[str, str]:
    """Environment for running demos.

    Uses the process environment, but drops an in-tree ``src/PyGEL`` entry
    from ``PYTHONPATH`` so the pip-installed ``pygel3d`` is used.
    """
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    if not existing:
        return env

    in_tree = IN_TREE_PYGEL.resolve()
    kept = []
    for part in existing.split(os.pathsep):
        if not part:
            continue
        try:
            if Path(part).resolve() == in_tree:
                continue
        except OSError:
            pass
        kept.append(part)

    if kept:
        env["PYTHONPATH"] = os.pathsep.join(kept)
    else:
        env.pop("PYTHONPATH", None)
    return env


def ensure_pygel(python: str) -> str | None:
    """Return the imported ``pygel3d`` path if PATH's Python has it installed."""
    code, out = run_cmd(
        [python, "-c", "import pygel3d; print(pygel3d.__file__)"],
        cwd=DEMO_ROOT,
        env=demo_env(),
    )
    if code != 0:
        return None
    location = out.strip().splitlines()[-1] if out.strip() else ""
    return location or python


def ensure_gel_package(jobs: int) -> bool:
    """Build GEL and install it as a CMake package under ~/.local."""
    if gel_package_config() is not None:
        return True

    print("Installed GEL package not found; building and installing locally …")
    GEL_BUILD_DIR.mkdir(parents=True, exist_ok=True)
    code, out = run_cmd(
        [
            "cmake",
            "-S",
            str(REPO_ROOT),
            "-B",
            str(GEL_BUILD_DIR),
            f"-DCMAKE_INSTALL_PREFIX={GEL_INSTALL_DIR}",
        ],
        cwd=REPO_ROOT,
        stream=True,
    )
    if code != 0:
        print("Failed to configure the main GEL project.")
        if out:
            indent_output(out)
        return False

    code, out = run_cmd(
        ["cmake", "--build", str(GEL_BUILD_DIR), "--parallel", str(jobs)],
        cwd=REPO_ROOT,
        stream=True,
    )
    if code != 0:
        print("Failed to build the main GEL project.")
        if out:
            indent_output(out)
        return False

    code, out = run_cmd(
        ["cmake", "--install", str(GEL_BUILD_DIR)],
        cwd=REPO_ROOT,
        stream=True,
    )
    if code != 0:
        print("Failed to install the GEL CMake package.")
        if out:
            indent_output(out)
        return False
    return gel_package_config() is not None


# ---------------------------------------------------------------------------
# CMake + run
# ---------------------------------------------------------------------------

def build_cmake_demo(demo_dir: Path, jobs: int) -> tuple[bool, str]:
    build_dir = demo_dir / "build"
    build_dir.mkdir(parents=True, exist_ok=True)

    cmake_prefix = str(GEL_INSTALL_DIR)
    existing_prefix = os.environ.get("CMAKE_PREFIX_PATH")
    if existing_prefix:
        cmake_prefix = os.pathsep.join((cmake_prefix, existing_prefix))

    code, out = run_cmd(
        [
            "cmake",
            "-S",
            str(demo_dir),
            "-B",
            str(build_dir),
            f"-DCMAKE_PREFIX_PATH={cmake_prefix}",
        ],
        cwd=demo_dir,
        stream=True,
    )
    if code != 0:
        return False, "cmake configure failed" + (f": {out.strip().splitlines()[-1]}" if out.strip() else "")

    code, out = run_cmd(
        ["cmake", "--build", str(build_dir), "--config", "Release", "--parallel", str(jobs)],
        cwd=demo_dir,
        stream=True,
    )
    if code != 0:
        return False, "cmake build failed" + (f": {out.strip().splitlines()[-1]}" if out.strip() else "")
    return True, ""


def run_binary(binary: Path, demo_dir: Path, timeout: float | None) -> StepResult:
    started = time.perf_counter()
    code, out = run_cmd([str(binary)], cwd=demo_dir, env=demo_env(), timeout=timeout)
    elapsed = time.perf_counter() - started

    if code is None:
        return StepResult(demo_dir.name, "binary", binary.name, "timeout",
                          f"stopped after {timeout:g}s", elapsed)
    if code != 0:
        tail = out.strip().splitlines()[-1] if out.strip() else f"exit {code}"
        indent_output(out)
        return StepResult(demo_dir.name, "binary", binary.name, "fail",
                          tail, elapsed)
    if out.strip():
        indent_output(out)
    return StepResult(demo_dir.name, "binary", binary.name, "ok", seconds=elapsed)


def run_python(script: Path, demo_dir: Path, python: str, timeout: float | None) -> StepResult:
    started = time.perf_counter()
    code, out = run_cmd(
        [python, str(script)],
        cwd=demo_dir,
        env=demo_env(),
        timeout=timeout,
    )
    elapsed = time.perf_counter() - started

    if code is None:
        return StepResult(demo_dir.name, "python", script.name, "timeout",
                          f"stopped after {timeout:g}s", elapsed)
    if code != 0:
        tail = out.strip().splitlines()[-1] if out.strip() else f"exit {code}"
        indent_output(out)
        return StepResult(demo_dir.name, "python", script.name, "fail",
                          tail, elapsed)
    if out.strip():
        indent_output(out)
    return StepResult(demo_dir.name, "python", script.name, "ok", seconds=elapsed)


# ---------------------------------------------------------------------------
# Per-directory driver
# ---------------------------------------------------------------------------

def process_demo(
    demo_dir: Path,
    harness: Harness,
    *,
    jobs: int,
    timeout: float | None,
    skip_cmake: bool,
    skip_python: bool,
    python: str,
) -> None:
    cmake_lists = demo_dir / "CMakeLists.txt"
    scripts = [] if skip_python else python_scripts(demo_dir)
    has_cmake = cmake_lists.is_file() and not skip_cmake

    print()
    print(f"── {demo_dir.name} ──")

    if not has_cmake and not scripts:
        harness.add(StepResult(demo_dir.name, "demo", demo_dir.name, "skip",
                               "no CMakeLists.txt and no Python scripts"))
        return

    if has_cmake:
        print("  configuring and building …")
        started = time.perf_counter()
        ok, detail = build_cmake_demo(demo_dir, jobs)
        elapsed = time.perf_counter() - started
        if not ok:
            harness.add(StepResult(demo_dir.name, "cmake", "build", "fail", detail, elapsed))
        else:
            harness.add(StepResult(demo_dir.name, "cmake", "build", "ok", seconds=elapsed))
            names = cmake_executables(cmake_lists)
            if not names:
                harness.add(StepResult(demo_dir.name, "binary", "(none)", "fail",
                                       "no add_executable() in CMakeLists.txt"))
            for name in names:
                binary = find_built_binary(demo_dir / "build", name)
                if binary is None:
                    harness.add(StepResult(demo_dir.name, "binary", name, "fail",
                                           "executable not found after build"))
                    continue
                print(f"  running {binary.name} …")
                harness.add(run_binary(binary, demo_dir, timeout))

    for script in scripts:
        print(f"  running {script.name} …")
        harness.add(run_python(script, demo_dir, python, timeout))


def list_demos(only: list[str] | None) -> None:
    print(f"Demos in {DEMO_ROOT}:")
    for demo_dir in demo_dirs(only):
        kinds = []
        if (demo_dir / "CMakeLists.txt").is_file():
            exes = cmake_executables(demo_dir / "CMakeLists.txt")
            kinds.append("cmake [" + ", ".join(exes) + "]" if exes else "cmake")
        scripts = python_scripts(demo_dir)
        if scripts:
            kinds.append("python [" + ", ".join(p.name for p in scripts) + "]")
        if not kinds:
            kinds.append("nothing to run")
        print(f"  {demo_dir.name:32} {', '.join(kinds)}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str]) -> argparse.Namespace:
    cpu = os.cpu_count() or 4
    parser = argparse.ArgumentParser(
        description="Build and run every demo under src/demo.",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="DIR",
        help="restrict the run to these demo directory names",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=120.0,
        help="seconds to wait for each binary/script before stopping it "
             "(default: 120; use 0 for no timeout)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=cpu,
        help=f"parallel build jobs (default: {cpu})",
    )
    parser.add_argument(
        "--skip-cmake",
        action="store_true",
        help="do not configure, build, or run C++ demos",
    )
    parser.add_argument(
        "--skip-python",
        action="store_true",
        help="do not run Python demos",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list discovered demos and exit",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    timeout = None if args.timeout <= 0 else args.timeout

    if args.list:
        list_demos(args.only)
        return 0

    if args.only:
        known = {p.name for p in demo_dirs()}
        missing = [name for name in args.only if name not in known]
        if missing:
            print(f"Unknown demo director{'y' if len(missing) == 1 else 'ies'}: "
                  + ", ".join(missing), file=sys.stderr)
            return 2

    python = path_python()
    if python is None:
        print("No python3 or python interpreter found on PATH.", file=sys.stderr)
        return 2

    print(f"GEL demo harness")
    print(f"  repo:    {REPO_ROOT}")
    print(f"  demos:   {DEMO_ROOT}")
    print(f"  python:  {python}")
    print(f"  timeout: {'none' if timeout is None else f'{timeout:g}s'}")
    print(f"  jobs:    {args.jobs}")

    if not args.skip_cmake:
        if not ensure_gel_package(args.jobs):
            print("Cannot build demos without an installed GEL package. "
                  "Build and install the main project first, for example:\n"
                  "  cmake -S . -B build\n"
                  "  cmake --build build && cmake --install build\n"
                  "GEL installs to ~/.local by default.",
                  file=sys.stderr)
            return 2

    pygel_loc = None
    if not args.skip_python:
        pygel_loc = ensure_pygel(python)
        if pygel_loc is None:
            print("PATH Python cannot import the pip-installed pygel3d package.\n"
                  "Install it first, for example:\n"
                  "  sh build_pygel.sh\n"
                  "  pip install PyGEL3D",
                  file=sys.stderr)
            return 2
        print(f"  pygel3d: {pygel_loc}")

    harness = Harness()
    for demo_dir in demo_dirs(args.only):
        process_demo(
            demo_dir,
            harness,
            jobs=args.jobs,
            timeout=timeout,
            skip_cmake=args.skip_cmake,
            skip_python=args.skip_python,
            python=python,
        )
    return harness.summary()


if __name__ == "__main__":
    sys.exit(main())
