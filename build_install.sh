#!/bin/sh
# Build GEL/PyGEL, install the C++ package to ~/.local, then create and
# pip-install the PyGEL3D wheel.
#
# The native build is Release unless --debug or -D is given. --release and
# -R select the default explicitly. An existing CMake cache is used only
# when its build type already matches; otherwise CMake is configured again.
# A matching cache is not cleaned, so an up-to-date tree is not rebuilt.
set -eu
cd "$(dirname "$0")"

requested=
while [ $# -gt 0 ]; do
    case "$1" in
        -D|--debug)
            if [ "$requested" = Release ]; then
                echo "build_install.sh: --debug and --release cannot be combined" >&2
                exit 1
            fi
            requested=Debug
            ;;
        -R|--release)
            if [ "$requested" = Debug ]; then
                echo "build_install.sh: --debug and --release cannot be combined" >&2
                exit 1
            fi
            requested=Release
            ;;
        *)
            echo "build_install.sh: unknown argument: $1" >&2
            echo "usage: build_install.sh [-R|--release] [-D|--debug]" >&2
            exit 1
            ;;
    esac
    shift
done
if [ -z "$requested" ]; then
    requested=Release
fi

cache=build/CMakeCache.txt
cache_type=
cache_prefix=
generator=
if [ -f "$cache" ]; then
    cache_type=$(sed -n 's/^CMAKE_BUILD_TYPE:[^=]*=//p' "$cache")
    cache_prefix=$(sed -n 's/^CMAKE_INSTALL_PREFIX:[^=]*=//p' "$cache")
    generator=$(sed -n 's/^CMAKE_GENERATOR:INTERNAL=//p' "$cache")
fi

multi=0
case "$generator" in
    "Visual Studio"*|Xcode)
        multi=1
        ;;
esac

# Single-config generators bake the build type into compile flags. Reconfigure
# when that cached type is missing or different. Leave a matching cache alone
# so cmake --build stays incremental. Multi-config generators select the type
# with --config at build time.
if [ "$multi" -eq 0 ] && { [ ! -f "$cache" ] || [ "$cache_type" != "$requested" ] || [ "$cache_prefix" != "${HOME}/.local" ]; }; then
    if [ -f "$cache" ] && [ "$cache_type" != "$requested" ]; then
        echo "CMake cache build type is ${cache_type:-unset}; configuring ${requested}"
    else
        echo "Configuring ${requested}"
    fi
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE="${requested}" \
        -DCMAKE_INSTALL_PREFIX="${HOME}/.local"
    generator=$(sed -n 's/^CMAKE_GENERATOR:INTERNAL=//p' "$cache")
    multi=0
    case "$generator" in
        "Visual Studio"*|Xcode)
            multi=1
            ;;
    esac
else
    echo "CMake build type is already ${requested}"
fi

if [ "$multi" -eq 1 ]; then
    cmake --build build -j 12 --config "${requested}"
    cmake --install build --config "${requested}"
else
    cmake --build build -j 12
    cmake --install build
fi

rm -fr dist
# Reuse the CMake tree just built; the Python build copies only this
# platform's library into a py3-none-<platform> wheel.
python -m pip install -q build setuptools cmake
PYGEL_SKIP_CMAKE=1 python -m build --wheel --no-isolation
pip uninstall --yes PyGEL3D
pip install dist/*whl
