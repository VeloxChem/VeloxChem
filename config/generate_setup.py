# -*- coding: utf-8 -*-

#!/usr/bin/env python3

import os
import platform
import re
import subprocess
import sys
from distutils.sysconfig import get_config_var, get_python_inc
from pathlib import Path

# we import these modules to get their include directories
import mpi4py
import numpy
import pybind11


class SearchReplace(dict):
    """All-in-one multiple-string-substitution class."""

    def _make_regex(self):
        """Build re object based on the keys of the current dictionary."""
        return re.compile("|".join(map(re.escape, self.keys())))

    def __call__(self, match):
        """Handler invoked for each regex match."""
        return self[match.group(0)]

    def replace(self, text):
        """Translate text, returns the modified text."""
        return self._make_regex().sub(self, text)


def find_exe(executables):
    for exe in executables:
        for path in os.environ["PATH"].split(os.pathsep):
            fname = os.path.join(path, exe)
            if os.path.isfile(fname) and os.access(fname, os.X_OK):
                return exe, path
    return None, None


def get_command_output(command):
    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        print()
        print("*** Error: Unable to execute '{}'".format(" ".join(command)))
        sys.exit(1)
    return output.decode("utf-8")


def find_avx_linux():
    cpuinfo = os.path.join(os.sep, "proc", "cpuinfo")
    if os.path.isfile(cpuinfo):
        for avx in ["avx512", "avx2", "avx"]:
            with open(cpuinfo, "r") as fh:
                for line in fh:
                    if line[:5] == "flags" and avx in line:
                        return avx
    return None


def find_avx_macos():
    output = get_command_output(["sysctl", "-a"])
    lines = output.split(os.linesep)
    for line in lines:
        if "machdep.cpu.leaf7_features" in line and "AVX2" in line:
            return "avx2"
    for line in lines:
        if "machdep.cpu.features" in line and "AVX1" in line:
            return "avx"
    return None


def find_mkl_avx(is_linux, is_macos):

    print("*** Checking avx... ", end="")
    if is_linux:
        avx = find_avx_linux()
        if check_pdc_beskow():
            avx = "avx2"
    elif is_macos:
        avx = find_avx_macos()
    print(avx)

    if avx is not None:
        mkl_avx = "-lmkl_{}".format(avx)
    else:
        mkl_avx = "-lmkl_def"
    return mkl_avx


def check_ubuntu():
    for name in ["lsb-release", "os-release"]:
        fname = os.path.join(os.sep, "etc", name)
        if os.path.isfile(fname):
            with open(fname, "r") as fh:
                for line in fh:
                    if "ubuntu" in line.lower():
                        return True
    return False


def check_pdc_beskow():
    is_pdc = "SNIC_SITE" in os.environ and os.environ["SNIC_SITE"] == "pdc"
    is_beskow = (
        "SNIC_RESOURCE" in os.environ and os.environ["SNIC_RESOURCE"] == "beskow"
    )
    return is_pdc and is_beskow


def check_dir(dirname, label):
    if not os.path.isdir(dirname):
        print("*** Error: {} dir {} does not exist!".format(label, dirname))
        sys.exit(1)


def check_file(filename, label):
    if not os.path.isfile(filename):
        print("*** Error: {} file {} does not exist!".format(label, filename))
        sys.exit(1)


def generate_setup(
    template_file, setup_file, user_flag=None, build_lib=Path("build/lib")
):

    if isinstance(template_file, str):
        template_file = Path(template_file)

    if isinstance(setup_file, str):
        setup_file = Path(setup_file)

    ext_suffix = get_config_var("EXT_SUFFIX") or get_config_var("SO")

    # OS information

    print("*** Checking operating system... ", end="")

    is_linux = "Linux" == platform.system()
    is_macos = "Darwin" == platform.system()
    if not (is_linux or is_macos):
        print()
        print("*** Error: Unsupported OS!")
        print("***        Only Linux and MacOS are supported.")
        sys.exit(1)

    is_ubuntu = is_linux and check_ubuntu()

    if is_linux:
        print("Linux")
    elif is_macos:
        print("MacOS")

    # compiler information

    print("*** Checking c++ compiler... ", end="")

    if "CRAYPE_VERSION" in os.environ and "CXX" in os.environ:
        cxx, cxx_path = find_exe([os.environ["CXX"]])
    else:
        if isinstance(user_flag, str) and user_flag.lower() == "gnu":
            cxx, cxx_path = find_exe(["mpicxx", "mpiicpc", "mpiCXX"])
        else:
            cxx, cxx_path = find_exe(["mpiicpc", "mpicxx", "mpiCXX"])

    print(cxx)

    if cxx is None:
        print("*** Error: Unable to find c++ compiler!")
        if "CRAYPE_VERSION" in os.environ and "CXX" in os.environ:
            print("***        Please make sure that CXX is correctly set.")
        else:
            print("***        Please make sure that mpiicpc, mpicxx, or")
            print("***        mpiCXX is in your PATH.")
        sys.exit(1)

    if cxx in ["icpc", "g++", "clang++"]:
        print("*** Error: {} is not a MPI compiler!".format(cxx))
        sys.exit(1)

    if cxx in ["mpiicpc", "mpicxx", "mpiCXX"]:
        cxxname = get_command_output([cxx, "-show"])
    else:
        cxxname = get_command_output([cxx, "--version"])
    if "Cray clang" in cxxname:
        cxxname = cxxname.replace("Cray clang", "Crayclang")
    cxxname = cxxname.split()[0]

    if cxxname in ["icc", "gcc", "clang"]:
        print("*** Error: {} is not a c++ compiler!".format(cxx))
        sys.exit(1)

    use_intel = cxxname == "icpc"
    use_gnu = re.match(r"(.*(c|g|gnu-c)\+\+)", cxxname)
    use_clang = cxxname in ["clang++", "Crayclang"] or re.match(
        r".*-clang\+\+", cxxname
    )

    if not (use_intel or use_gnu or use_clang):
        print("*** Error: Unrecognized c++ compiler!")
        print("***        Only Intel, GNU, and Clang compilers are supported.")
        sys.exit(1)

    # cxx and omp flags

    if use_intel:
        cxx_flags = "-xHost -qopenmp"
        if check_pdc_beskow():
            cxx_flags = "-qopenmp"
        omp_flag = "-liomp5"
    elif use_gnu:
        cxx_flags = "-fopenmp"
        omp_flag = "-lgomp"
    elif use_clang:
        cxx_flags = "-Xpreprocessor -fopenmp"
        omp_flag = "-lomp"

    # math library

    print("*** Checking math library... ", end="")

    # check conda environment
    is_conda = os.path.exists(os.path.join(sys.prefix, "conda-meta"))

    # check whether MKL is in conda environment
    if "MKLROOT" not in os.environ:
        has_lib = (
            Path(sys.prefix, "lib/libmkl_core.so").is_file()
            or Path(sys.prefix, "lib/libmkl_core.dylib").is_file()
        )
        has_header = Path(sys.prefix, "include/mkl.h").is_file()
        if is_conda and has_lib and has_header:
            os.environ["MKLROOT"] = sys.prefix

    # check whether OpenBLAS is in conda environment
    if "OPENBLASROOT" not in os.environ:
        has_lib = (
            Path(sys.prefix, "lib/libopenblas.so").is_file()
            or Path(sys.prefix, "lib/libopenblas.dylib").is_file()
        )
        has_header = (
            Path(sys.prefix, "include/lapacke.h").is_file()
            and Path(sys.prefix, "include/cblas.h").is_file()
        )
        if is_conda and has_lib and has_header:
            os.environ["OPENBLASROOT"] = sys.prefix

    use_mkl = "MKLROOT" in os.environ
    use_openblas = "OPENBLASROOT" in os.environ
    use_craylibsci = "CRAY_LIBSCI_VERSION" in os.environ

    if not (use_mkl or use_openblas or use_craylibsci):
        print()
        print("*** Error: Unable to find math library!")
        print("***        Please make sure that you have set MKLROOT or")
        print("***        OPENBLASROOT. OpenBLAS can be downloaded from")
        print("***        https://github.com/xianyi/OpenBLAS")
        sys.exit(1)

    # mkl flags

    if use_mkl:
        print("MKL")

        mkl_inc = os.path.join(os.environ["MKLROOT"], "include")
        check_dir(mkl_inc, "mkl include")

        mkl_dir = os.path.join(os.environ["MKLROOT"], "lib", "intel64")
        if not os.path.isdir(mkl_dir):
            mkl_dir = os.path.join(os.environ["MKLROOT"], "lib")
        check_dir(mkl_dir, "mkl lib")

        math_lib = f"MATH_INC := -isystem {mkl_inc}"
        math_lib += os.linesep + f"MATH_LIB := -L{mkl_dir}"
        if is_macos:
            math_lib += (
                os.linesep + f"MATH_LIB += -Wl,-rpath,{mkl_dir} -lmkl_rt -lpthread -lm -ldl"
            )
        else:
            math_lib += (
                os.linesep + f"MATH_LIB += -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl"
            )

        if is_ubuntu and not use_intel:
            conf = {
                "@_mkl_interface_layer_@": "MKL_INTERFACE_LP64+MKL_INTERFACE_GNU",
                "@_mkl_threading_layer_@": "MKL_THREADING_GNU",
            }
        elif use_intel:
            conf = {
                "@_mkl_interface_layer_@": "MKL_INTERFACE_LP64",
                "@_mkl_threading_layer_@": "MKL_THREADING_INTEL",
            }
        elif use_clang:
            conf = {
                "@_mkl_interface_layer_@": "MKL_INTERFACE_LP64",
                "@_mkl_threading_layer_@": "MKL_THREADING_INTEL",
            }
        replacer = SearchReplace(conf)

        # read in src/general/ConfigMKL.hpp.in
        conf_mkl_in = Path("src/general/ConfigMKL.hpp.in")
        conf_mkl = Path("src/general/ConfigMKL.hpp")

        with conf_mkl_in.open("r") as f:
            contents = "".join(f.readlines())

        with conf_mkl.open("w") as f:
            f.write(replacer.replace(contents))

    # openblas flags

    elif use_openblas:
        print("OpenBLAS")

        openblas_inc = os.path.join(os.environ["OPENBLASROOT"], "include")
        check_dir(openblas_inc, "openblas include")

        openblas_dir = os.path.join(os.environ["OPENBLASROOT"], "lib")
        check_dir(openblas_dir, "openblas lib")

        math_lib = "MATH_INC := -isystem {}".format(openblas_inc)
        math_lib += os.linesep + "MATH_LIB := -L{}".format(openblas_dir)
        math_lib += os.linesep + "MATH_LIB += -Wl,-rpath,{}".format(openblas_dir)
        openblas_flag = "-lopenblas"
        if use_intel:
            openblas_flag += " -lifcore"
        math_lib += os.linesep + "MATH_LIB += {} {} {}".format(
            openblas_flag, omp_flag, "-lpthread -lm -ldl"
        )

    # cray-libsci flags

    elif use_craylibsci:
        print("Cray LibSci")

        math_lib = "MATH_INC := "
        math_lib += os.linesep + "MATH_LIB := {} {}".format(
            omp_flag, "-lpthread -lm -ldl"
        )

    # extra flags for mac

    maclibs = ""
    if is_macos:
        maclibs = "-undefined dynamic_lookup"

    # lto flag

    lto_flag = ""
    if use_gnu:
        lto_flag = "-fno-lto"

    # xtb package

    use_xtb = False
    xtb_root = os.getenv("XTBHOME", sys.prefix)
    # include
    xtb_inc = Path(xtb_root, "include/xtb")
    has_xtb_header = xtb_inc.is_dir() and (xtb_inc / "xtb.h").is_file()
    # library
    has_xtb_lib = False
    _lib64 = Path(xtb_root, "lib64")
    if (_lib64 / "libxtb.so").is_file() or (_lib64 / "libxtb.dylib").is_file():
        xtb_dir = _lib64
        has_xtb_lib = True
    else:
        xtb_dir = Path(xtb_root, "lib")
        has_xtb_lib = (xtb_dir / "libxtb.so").is_file() or (
            xtb_dir / "libxtb.dylib"
        ).is_file()
    # support files
    xtb_path = Path(xtb_root, "share/xtb")
    xtb_params = ["param_gfn0-xtb.txt", "param_gfn1-xtb", "param_gfn2-xtb.txt"]
    has_xtb_share = xtb_path.is_dir() and all(
        [(xtb_path / f"param_gfn{x}-xtb.txt").is_file() for x in range(3)]
    )

    if has_xtb_header and has_xtb_lib and has_xtb_share:
        use_xtb = True

        xtb_lib = f"XTB_INC := -I{str(xtb_inc)}\n"
        xtb_lib += f"XTB_LIB := -L{str(xtb_dir)} -Wl,-rpath,{str(xtb_dir)} -lxtb\n"
        xtb_lib += f"XTB_PATH := {str(xtb_path)}\n"

    # google test lib

    gtest_root = os.getenv("GTESTROOT", sys.prefix)
    # include
    gtest_inc = Path(gtest_root, "include")
    has_gtest_header = gtest_inc.is_dir() and (gtest_inc / "gtest/gtest.h").is_file()
    # library
    _lib = Path(f"{gtest_root}/lib")
    has_gtest_lib = False
    gtest_lib = None
    if _lib.is_dir():
        if (_lib / "libgtest.a").is_file():
            has_gtest_lib = True
            gtest_lib = _lib / "libgtest.a"
        elif (_lib / "libgtest.so").is_file():
            has_gtest_lib = True
            gtest_lib = _lib / "libgtest.so"
        elif (_lib / "libgtest.dylib").is_file():
            has_gtest_lib = True
            gtest_lib = _lib / "libgtest.dylib"

    if not (has_gtest_header and has_gtest_lib):
        gtest_root = None
        gtest_inc = None
        gtest_lib = None

    # print Makefile.setup

    with template_file.open("r", encoding='utf-8') as f_temp:
        lines = f_temp.readlines()

    with setup_file.open("w", encoding='utf-8') as f_mkfile:
        for line in lines:
            if "====placeholder====" in line:
                print("# Automatically generated settings", file=f_mkfile)
                print("", file=f_mkfile)

                print("# Python-related headers files", file=f_mkfile)
                print(f"PYTHON_INC = -isystem {get_python_inc()}", file=f_mkfile)
                print(f"PYTHON_INC += -isystem {mpi4py.get_include()}", file=f_mkfile)
                print(f"PYTHON_INC += -isystem {numpy.get_include()}", file=f_mkfile)
                print(f"PYTHON_INC += -isystem {pybind11.get_include()}", file=f_mkfile)

                build_lib_str = (build_lib / "veloxchem").resolve()
                vlx_target_str = "$(BUILD_LIB)/veloxchemlib" + ext_suffix
                print("BUILD_LIB := {}".format(build_lib_str), file=f_mkfile)
                print("VLX_TARGET := {}".format(vlx_target_str), file=f_mkfile)
                print("", file=f_mkfile)

                print("USE_MPI := true", file=f_mkfile)
                print(
                    "USE_MKL := {}".format("true" if use_mkl else "false"),
                    file=f_mkfile,
                )
                print(
                    "USE_XTB := {}".format("true" if use_xtb else "false"),
                    file=f_mkfile,
                )
                print("", file=f_mkfile)

                print(math_lib, file=f_mkfile)
                print("", file=f_mkfile)

                print(
                    "PYTHON :=",
                    "python{}.{}".format(sys.version_info[0], sys.version_info[1]),
                    file=f_mkfile,
                )
                print("", file=f_mkfile)

                print("CXX :=", cxx, file=f_mkfile)
                print("", file=f_mkfile)

                print("CXX_REL_FLG :=", cxx_flags, file=f_mkfile)
                print("CXX_DEB_FLG :=", cxx_flags, file=f_mkfile)
                print("", file=f_mkfile)

                print("MACLIBS :=", maclibs, file=f_mkfile)
                print("", file=f_mkfile)

                print("LTOFLAG :=", lto_flag, file=f_mkfile)
                print("", file=f_mkfile)

                if use_xtb:
                    print(xtb_lib, file=f_mkfile)
                    print("", file=f_mkfile)

                if gtest_root is not None and gtest_lib is not None:
                    print("GST_ROOT :=", gtest_root, file=f_mkfile)
                    print("GST_LIB :=", gtest_lib, file=f_mkfile)
                    print("", file=f_mkfile)

                print("# Generic settings", file=f_mkfile)
            else:
                print(line, end="", file=f_mkfile)

    print("*** Successfully generated {}".format(setup_file))
    sys.stdout.flush()


if __name__ == "__main__":

    template_file = "Setup.template"
    setup_file = os.path.join(os.pardir, "src", "Makefile.setup")

    if not os.path.isfile(template_file):
        template_file = os.path.join("config", "Setup.template")
        setup_file = os.path.join("src", "Makefile.setup")

    if not os.path.isfile(template_file):
        print("*** Error: Cannot find template file {}".format(template_file))
        sys.exit(1)

    user_flag = None
    if len(sys.argv) > 1:
        user_flag = sys.argv[1]

    generate_setup(template_file, setup_file, user_flag)
