from setuptools import setup, Extension
import os
import shutil
import sys
import subprocess

if os.name != "nt":
    # Do not force GCC on macOS; prefer Apple clang to avoid SDK/header issues.
    if sys.platform == "darwin":
        os.environ.setdefault("CC", "clang")
        os.environ.setdefault("CXX", "clang++")
    else:
        os.environ.setdefault("CC", "gcc")

VERSION = "1.0.0"

# Define the sources and the corresponding C flags
sources = [
    "src/booster/src/hashtables_bfields.c",
    "src/booster/src/transfer_distance.c",
    "src/booster/src/M1M2_hashmap.c",
    "src/booster/src/tree.c",
    "src/booster/src/stats.c",
    "src/booster/src/prng.c",
    "src/booster/src/hashmap.c",
    "src/booster/src/version.c",
    "src/booster/src/sort.c",
    "src/booster/src/io.c",
    "src/booster/src/tree_utils.c",
    "src/booster/src/bitset_index.c",
    "src/booster/src/rapid_transfer.c",
    "src/booster/src/heavy_paths.c",
    "src/booster/src/debug.c",
    "src/booster/src/kludge.c",
    "src/booster/src/node_stack.c",
    "src/booster/src/split.c",
    "src/booster/src/booster.c",  # main file
]

def _env_flag(name: str, default: bool) -> bool:
    """Parse a boolean environment variable."""
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


# OpenMP frequently fails with Apple clang on macOS (`clang: unsupported option '-fopenmp'`).
# Default to no OpenMP on macOS, but allow users to force-enable it if they have a compatible toolchain.
use_openmp = _env_flag("CONSENSUS_USE_OPENMP", default=(sys.platform != "darwin"))

base_compile_flags = ["-Wall", "-O3", "-g", "-fPIC", "-DNDEBUG", f'-DVERSION="{VERSION}"']
compile_flags = list(base_compile_flags)
link_flags = ["-lm"]

if use_openmp:
    if sys.platform == "darwin":
        libomp_prefix = None
        if shutil.which("brew"):
            try:
                libomp_prefix = subprocess.check_output(
                    ["brew", "--prefix", "libomp"],
                    text=True,
                ).strip()
            except subprocess.SubprocessError:
                libomp_prefix = None
        if not libomp_prefix:
            raise RuntimeError(
                "CONSENSUS_USE_OPENMP=1 on macOS requires Homebrew libomp. "
                "Install with: brew install libomp"
            )
        compile_flags.extend(["-Xpreprocessor", "-fopenmp", f"-I{libomp_prefix}/include"])
        link_flags.extend([f"-L{libomp_prefix}/lib", "-lomp"])
    else:
        compile_flags.append("-fopenmp")
        link_flags.append("-fopenmp")

if sys.platform == "darwin" and not use_openmp:
    print("PhyloCRISP: building booster extension without OpenMP on macOS (set CONSENSUS_USE_OPENMP=1 to force-enable).")

build_booster = _env_flag("CONSENSUS_BUILD_BOOSTER", default=(os.name != "nt"))

ext_modules = []
if build_booster:
    booster_extension = Extension(
        "booster",  # Name of the extension
        sources=sources,
        extra_compile_args=compile_flags,
        extra_link_args=link_flags,
    )
    ext_modules = [booster_extension]
elif sys.platform.startswith("win"):
    print(
        "PhyloCRISP: skipping booster extension build on Windows by default. "
        "Fast prune-only transfer methods will be unavailable unless booster is provided separately."
    )
else:
    print("PhyloCRISP: skipping booster extension build (CONSENSUS_BUILD_BOOSTER=0).")

setup(
    name="PhyloCRISP",
    version="0.0.1",
    ext_modules=ext_modules,  # Booster is optional on some platforms (e.g., Windows stage-0 install)
)
