"""Build PyGEL3D as a ctypes wheel tagged py3-none-<platform>.

The native library is produced with CMake and loaded at import time; it does
not link against libpython, so one wheel per OS/arch covers every Python 3.11+.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution

try:
    from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:  # pragma: no cover
    from wheel.bdist_wheel import bdist_wheel

ROOT = Path(__file__).resolve().parent
PKG_DIR = ROOT / "src" / "PyGEL" / "pygel3d"
# Keep CMake output out of setuptools' build/lib and build/bdist.* trees.
NATIVE_BUILD_DIR = ROOT / "build" / "native"
_NATIVE_NAMES = ("libPyGEL.so", "libPyGEL.dylib", "PyGEL.dll")


def native_lib_name() -> str:
    if sys.platform == "darwin":
        return "libPyGEL.dylib"
    if sys.platform == "win32":
        return "PyGEL.dll"
    return "libPyGEL.so"


def _is_pygel_lib(path: Path) -> bool:
    name = path.name
    if name.startswith("libglfw") or name.startswith("glfw"):
        return False
    return name == native_lib_name()


def _find_built_lib(search_root: Path) -> Path | None:
    name = native_lib_name()
    candidates = [
        search_root / name,
        search_root / "Release" / name,
        search_root / "RelWithDebInfo" / name,
        search_root / "Debug" / name,
    ]
    for path in candidates:
        if path.is_file() and _is_pygel_lib(path):
            return path
    return None


def _remove_other_native_libs(directory: Path) -> None:
    keep = native_lib_name()
    if not directory.is_dir():
        return
    for name in _NATIVE_NAMES:
        path = directory / name
        if name != keep and path.is_file():
            path.unlink()


def _cmake() -> str:
    exe = shutil.which("cmake")
    if not exe:
        raise RuntimeError(
            "cmake is required to build PyGEL3D from source. "
            "Install CMake >= 3.25, or pip-install the 'cmake' package."
        )
    return exe


def _configure_and_build() -> Path:
    build_dir = NATIVE_BUILD_DIR
    build_dir.mkdir(parents=True, exist_ok=True)
    cmake = _cmake()
    configure = [
        cmake,
        "-S",
        str(ROOT),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_TESTING=OFF",
        "-DUse_GLGraphics=ON",
    ]
    if sys.platform.startswith("linux"):
        configure.append("-DOpenGL_GL_PREFERENCE=LEGACY")
    # Ninja on Windows often picks MinGW from PATH; CPython wheels must use
    # MSVC, so leave the Visual Studio generator as CMake's default there.
    if (
        sys.platform != "win32"
        and not (build_dir / "CMakeCache.txt").exists()
        and shutil.which("ninja")
    ):
        configure.extend(["-G", "Ninja"])
    print("PyGEL3D: configuring native library:", " ".join(configure))
    subprocess.run(configure, check=True)

    build_cmd = [
        cmake,
        "--build",
        str(build_dir),
        "--target",
        "PyGEL",
        "--config",
        "Release",
    ]
    print("PyGEL3D: building native library:", " ".join(build_cmd))
    subprocess.run(build_cmd, check=True)

    built = _find_built_lib(build_dir)
    if built is None:
        raise RuntimeError(
            f"CMake succeeded but {native_lib_name()} was not found under {build_dir}"
        )
    return built


def ensure_native_lib() -> Path:
    """Return the current-platform PyGEL library, building it if needed."""
    preset = os.environ.get("PYGEL_NATIVE_LIB")
    if preset:
        path = Path(preset)
        if not path.is_file():
            raise FileNotFoundError(f"PYGEL_NATIVE_LIB={preset} does not exist")
        return path

    skip = os.environ.get("PYGEL_SKIP_CMAKE", "") == "1"
    if skip:
        # build_install.sh compiles into build/ and then packages with
        # PYGEL_SKIP_CMAKE=1. That library is the one to ship. build/native
        # is a separate tree used only when setuptools runs CMake itself,
        # and a later pip build can make it newer than build/libPyGEL.dylib.
        built = _find_built_lib(ROOT / "build")
        if built is not None:
            return built
        built = _find_built_lib(NATIVE_BUILD_DIR)
        if built is not None:
            return built
    if skip:
        raise RuntimeError(
            "PYGEL_SKIP_CMAKE=1 but "
            f"{native_lib_name()} was not found in {NATIVE_BUILD_DIR} or {ROOT / 'build'}"
        )
    return _configure_and_build()


class BinaryDistribution(Distribution):
    """Treat the ctypes library as platform-specific so the wheel is platlib."""

    def has_ext_modules(self):
        return True


class BuildPyGEL(build_py):
    def run(self):
        lib = ensure_native_lib()
        dest = PKG_DIR / native_lib_name()
        dest.parent.mkdir(parents=True, exist_ok=True)
        _remove_other_native_libs(PKG_DIR)
        if lib.resolve() != dest.resolve():
            print(f"PyGEL3D: copying {lib} -> {dest}")
            shutil.copy2(lib, dest)
        super().run()
        _remove_other_native_libs(Path(self.build_lib) / "pygel3d")


class PlatformWheel(bdist_wheel):
    """Tag as py3-none-<platform> (ctypes; no libpython ABI)."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        python, abi, plat = super().get_tag()
        python, abi = "py3", "none"
        if sys.platform == "darwin":
            # CMakeLists.txt always builds a universal2 dylib (arm64;x86_64)
            # with MACOSX_DEPLOYMENT_TARGET 11.0.
            plat = re.sub(
                r"macosx_\d+_\d+_(arm64|x86_64|intel|universal2)",
                "macosx_11_0_universal2",
                plat,
            )
            if "macosx_" not in plat:
                plat = "macosx_11_0_universal2"
        return python, abi, plat


if sys.platform == "darwin":
    os.environ.setdefault("MACOSX_DEPLOYMENT_TARGET", "11.0")

setup(
    distclass=BinaryDistribution,
    cmdclass={
        "build_py": BuildPyGEL,
        "bdist_wheel": PlatformWheel,
    },
    package_data={"pygel3d": [native_lib_name()]},
)
