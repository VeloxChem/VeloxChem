#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

# -*- coding: utf-8 -*-

import os
import platform
import re
import subprocess
import sys
import sysconfig
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


def is_executable(exe):
    return Path(exe).is_file() and os.access(exe, os.X_OK)


def find_exe(executables):
    for exe in executables:
        if Path(exe).is_absolute() and is_executable(exe):
            return exe
        for p in os.environ["PATH"].split(os.pathsep):
            fname = str(Path(p) / exe)
            if is_executable(fname):
                return exe
    return None


def get_command_output(command):
    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        cmd_line = " ".join(command)
        print(f"\n*** Error: Failed to execute '{cmd_line}'")
        sys.exit(1)
    return output.decode("utf-8")


def check_cray():
    if "CRAYPE_VERSION" in os.environ and "CRAYPE_DIR" in os.environ:
        return True
    return False


def check_dir(dir_path, label):
    if not dir_path.is_dir():
        print(f"*** Error: {label} dir {dir_path} does not exist!")
        sys.exit(1)


def generate_setup(template_file, setup_file, build_lib=Path("build", "lib")):

    if isinstance(template_file, str):
        template_file = Path(template_file)

    if isinstance(setup_file, str):
        setup_file = Path(setup_file)

    ext_suffix = (sysconfig.get_config_var("EXT_SUFFIX") or
                  sysconfig.get_config_var("SO"))

    # ==> OS information <==

    print("*** Checking operating system... ", end="")

    is_linux = (platform.system() == "Linux")
    is_macos = (platform.system() == "Darwin")

    if is_linux:
        print("Linux")
    elif is_macos:
        print("MacOS")
    else:
        print()
        print("*** Error: Unsupported OS!")
        print("***        Only Linux and MacOS are supported.")
        sys.exit(1)

    # ==> compiler information <==

    print("*** Checking C++ compiler... ", end="")

    if check_cray():
        if "CXX" in os.environ and "MPICXX" not in os.environ:
            os.environ["MPICXX"] = os.environ["CXX"]

    if "MPICXX" in os.environ:
        cxx = find_exe([os.environ["MPICXX"]])
    else:
        cxx = find_exe(["mpiicpc", "mpicxx", "mpiCXX"])

    print(cxx)

    if cxx is None:
        print("*** Error: Unable to find C++ compiler!")
        print("***        Please make sure that MPICXX is correctly set.")
        sys.exit(1)

    if Path(cxx).name in ["icpc", "g++", "clang++"]:
        print(f"*** Error: {cxx} is not a MPI compiler!")
        sys.exit(1)

    if check_cray():
        cxxname = get_command_output([cxx, "--version"])
        if (cxxname.startswith("Cray clang") or
                cxxname.startswith("AMD clang")):
            cxxname = "clang++" if Path(cxx).name == "CC" else (
                "clang" if Path(cxx).name == "cc" else "Unknown")
    else:
        cxxname = get_command_output([cxx, "-show"])
    cxxname = cxxname.split()[0]

    if cxxname in ["icc", "gcc", "clang"]:
        print(f"*** Error: {cxx} is not a C++ compiler!")
        sys.exit(1)

    use_intel = (cxxname == "icpc")
    use_gnu = (cxxname == "g++" or
               re.match(r".*-(g|gnu-c)\+\+", cxxname) is not None)
    use_clang = (cxxname == "clang++" or
                 re.match(r".*-clang\+\+", cxxname) is not None)

    if not (use_intel or use_gnu or use_clang):
        print("*** Error: Unrecognized C++ compiler!")
        print("***        Only Intel, GNU, and Clang compilers are supported.")
        sys.exit(1)

    elif [use_intel, use_gnu, use_clang].count(True) != 1:
        print(f"*** Error: Unexpected C++ compiler: {cxxname}")
        print(f"***        use_intel = {use_intel}")
        print(f"***        use_gnu   = {use_gnu}")
        print(f"***        use_clang = {use_clang}")
        sys.exit(1)

    # ==> openmp flags <==

    if use_intel:
        if check_cray():
            cxx_flags = "-qopenmp"
        else:
            cxx_flags = "-xHost -qopenmp"
        omp_flag = "-liomp5"
    elif use_gnu:
        cxx_flags = "-fopenmp"
        omp_flag = "-lgomp"
    elif use_clang:
        cxx_flags = "-Xclang -fopenmp"
        omp_flag = "-lomp"

    # ==> math library <==

    print("*** Checking math library... ", end="")

    # check conda environment
    is_conda = Path(sys.prefix, "conda-meta").is_dir()

    # check whether MKL is in conda environment
    if is_conda and ("MKLROOT" not in os.environ):
        has_lib = (Path(sys.prefix, "lib", "libmkl_core.so").is_file() or
                   Path(sys.prefix, "lib", "libmkl_core.dylib").is_file())
        has_header = Path(sys.prefix, "include", "mkl.h").is_file()
        if has_lib and has_header:
            os.environ["MKLROOT"] = sys.prefix

    # check whether OpenBLAS is in conda environment
    if is_conda and ("OPENBLASROOT" not in os.environ):
        has_lib = (Path(sys.prefix, "lib", "libopenblas.so").is_file() or
                   Path(sys.prefix, "lib", "libopenblas.dylib").is_file())
        has_header = (Path(sys.prefix, "include", "lapacke.h").is_file() and
                      Path(sys.prefix, "include", "cblas.h").is_file())
        if has_lib and has_header:
            os.environ["OPENBLASROOT"] = sys.prefix

    use_mkl = ("MKLROOT" in os.environ)
    use_openblas = (("OPENBLASROOT" in os.environ) or
                    ("OPENBLAS_INCLUDE_DIR" in os.environ and
                     "OPENBLAS_LIBRARY" in os.environ))
    use_craylibsci = check_cray() and ("CRAY_LIBSCI_VERSION" in os.environ and
                                       "CRAY_LIBSCI_DIR" in os.environ)

    if not (use_mkl or use_openblas or use_craylibsci):
        print()
        print("*** Error: Unable to find math library!")
        print("***        Please make sure that you have set MKLROOT or")
        print("***        OPENBLASROOT (or OPENBLAS_INCLUDE_DIR and")
        print("***        OPENBLAS_LIBRARY). OpenBLAS can be downloaded")
        print("***        from https://github.com/xianyi/OpenBLAS")
        sys.exit(1)

    # ==> mkl flags <==

    if use_mkl:
        print("MKL")

        mkl_inc = Path(os.environ["MKLROOT"], "include")
        check_dir(mkl_inc, "mkl include")

        mkl_dir = Path(os.environ["MKLROOT"], "lib", "intel64")
        if not mkl_dir.is_dir():
            mkl_dir = Path(os.environ["MKLROOT"], "lib")
        check_dir(mkl_dir, "mkl lib")

        math_lib = f"MATH_INC := -I{mkl_inc}"
        math_lib += f"\nMATH_LIB := -L{mkl_dir}"
        if is_macos:
            math_lib += f"\nMATH_LIB += -Wl,-rpath,{mkl_dir} -lmkl_rt"
        else:
            math_lib += "\nMATH_LIB += -lmkl_rt -Wl,--no-as-needed"
        math_lib += " -lpthread -lm -ldl"

        if is_linux and not use_intel:
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
        conf_mkl_in = Path("src", "general", "ConfigMKL.hpp.in")
        conf_mkl = Path("src", "general", "ConfigMKL.hpp")

        with conf_mkl_in.open("r") as f:
            contents = "".join(f.readlines())

        with conf_mkl.open("w") as f:
            f.write(replacer.replace(contents))

    # ==> openblas flags <==

    elif use_openblas:
        print("OpenBLAS")

        if "OPENBLASROOT" in os.environ:
            openblas_inc = Path(os.environ["OPENBLASROOT"], "include")
        else:
            openblas_inc = Path(os.environ["OPENBLAS_INCLUDE_DIR"])
        check_dir(openblas_inc, "openblas include")

        if "OPENBLASROOT" in os.environ:
            openblas_dir = Path(os.environ["OPENBLASROOT"], "lib")
        else:
            openblas_dir = Path(os.environ["OPENBLAS_LIBRARY"])
        check_dir(openblas_dir, "openblas lib")

        math_lib = f"MATH_INC := -I{openblas_inc}"
        math_lib += f"\nMATH_LIB := -L{openblas_dir}"
        math_lib += f"\nMATH_LIB += -Wl,-rpath,{openblas_dir}"
        openblas_flag = "-lopenblas"
        if use_intel:
            openblas_flag += " -lifcore"
        math_lib += f"\nMATH_LIB += {openblas_flag} {omp_flag} -lpthread -lm -ldl"

    # ==> cray-libsci flags <==

    elif use_craylibsci:
        print("Cray LibSci")

        math_lib = "MATH_INC := "
        math_lib += f"\nMATH_LIB := {omp_flag} -lpthread -lm -ldl"

    # ==> extra flags for mac <==

    if is_macos:
        maclibs = "-undefined dynamic_lookup"
    else:
        maclibs = ""

    # ==> lto flag <==

    if use_gnu:
        lto_flag = "-fno-lto"
    else:
        lto_flag = ""

    # ==> xtb package <==

    xtb_root = os.getenv("XTBHOME", sys.prefix)

    # xtb include
    xtb_inc = Path(xtb_root, "include")
    if not xtb_inc.is_dir():
        xtb_inc = Path(xtb_root, "include", "xtb")
    has_xtb_header = (xtb_inc / "xtb.h").is_file()

    # xtb library
    xtb_dir = Path(xtb_root, "lib64")
    if not xtb_dir.is_dir():
        xtb_dir = Path(xtb_root, "lib")
    has_xtb_lib = ((xtb_dir / "libxtb.so").is_file() or
                   (xtb_dir / "libxtb.dylib").is_file())

    # xtb parameter files
    xtb_path = Path(xtb_root, "share", "xtb")
    xtb_params = [
        "param_gfn0-xtb.txt", "param_gfn1-xtb.txt", "param_gfn2-xtb.txt"
    ]
    has_xtb_share = all([(xtb_path / x).is_file() for x in xtb_params])

    use_xtb = (has_xtb_header and has_xtb_lib and has_xtb_share)
    if use_xtb:
        xtb_lib = f"XTB_INC := -I{xtb_inc}"
        xtb_lib += f"\nXTB_LIB := -L{xtb_dir}"
        xtb_lib += f"\nXTB_LIB += -Wl,-rpath,{xtb_dir} -lxtb"
        xtb_lib += f"\nXTB_PATH := {xtb_path}"
        print(f"*** Checking XTB... {xtb_root}")

    # ==> google test <==

    # use GTESTROOT for local from-source installations
    # this env-var takes priority!
    gtest_root = os.getenv("GTESTROOT", None)
    # use GTESTLIB for system-package (APT) installations
    gtest_lib = os.getenv("GTESTLIB", None)

    gtest_incdir = ""
    gtest_libdir = ""
    if gtest_root is not None:
        if Path(gtest_root, "include").is_dir():
            gtest_incdir = Path(gtest_root, "include")
        if Path(gtest_root, "lib").is_dir():
            gtest_libdir = Path(gtest_root, "lib")
    else:
        if gtest_lib is not None:
            gtest_libdir = gtest_lib

    print(f"*** Checking GoogleTest... {gtest_incdir} {gtest_libdir}")


    # ==> write Makefile.setup <==

    with template_file.open("r", encoding='utf-8') as f_temp:
        lines = f_temp.readlines()

    with setup_file.open("w", encoding='utf-8') as f_mkfile:
        for line in lines:
            if "====placeholder====" in line:
                print("# Automatically generated settings", file=f_mkfile)
                print("", file=f_mkfile)

                # MPI, MKL and XTB
                print("USE_MPI := true", file=f_mkfile)
                print("USE_MKL := {}".format("true" if use_mkl else "false"),
                      file=f_mkfile)
                print("USE_XTB := {}".format("true" if use_xtb else "false"),
                      file=f_mkfile)
                print("", file=f_mkfile)

                # build path
                build_lib_str = str((build_lib / "veloxchem").resolve())
                vlx_target_str = f"$(BUILD_LIB){os.sep}veloxchemlib{ext_suffix}"
                print(f"BUILD_LIB := {build_lib_str}", file=f_mkfile)
                print(f"VLX_TARGET := {vlx_target_str}", file=f_mkfile)
                print("", file=f_mkfile)

                # C++
                print("CXX :=", cxx, file=f_mkfile)
                print("", file=f_mkfile)
                print("CXX_REL_FLG :=", cxx_flags, file=f_mkfile)
                print("CXX_DEB_FLG :=", cxx_flags, file=f_mkfile)
                print("", file=f_mkfile)

                # math library
                print(math_lib, file=f_mkfile)
                print("", file=f_mkfile)

                # python
                python_version = f"{sys.version_info[0]}.{sys.version_info[1]}"
                print(f"PYTHON := python{python_version}", file=f_mkfile)
                print("", file=f_mkfile)
                python_include_path = sysconfig.get_path("include")
                print(f"PYTHON_INC := -I{python_include_path}", file=f_mkfile)
                print(f"PYTHON_INC += -I{mpi4py.get_include()}", file=f_mkfile)
                print(f"PYTHON_INC += -I{numpy.get_include()}", file=f_mkfile)
                print(f"PYTHON_INC += -I{pybind11.get_include()}",
                      file=f_mkfile)
                print("", file=f_mkfile)

                # additional settings
                print("MACLIBS :=", maclibs, file=f_mkfile)
                print("", file=f_mkfile)
                print("LTOFLAG :=", lto_flag, file=f_mkfile)
                print("", file=f_mkfile)
                if use_xtb:
                    print(xtb_lib, file=f_mkfile)
                    print("", file=f_mkfile)
                if gtest_incdir:
                    print(f"GST_INC := {gtest_incdir}", file=f_mkfile)
                if gtest_libdir:
                    print(f"GST_LIB := -L{gtest_libdir} -lgtest", file=f_mkfile)
                    print("", file=f_mkfile)

                print("# Generic settings", file=f_mkfile)
            else:
                print(line, end="", file=f_mkfile)

    print(f"*** Successfully generated {setup_file}")
    sys.stdout.flush()


if __name__ == "__main__":

    config_dir = Path(__file__).parent
    template_file = config_dir / "Setup.template"
    setup_file = config_dir.parent / "src" / "Makefile.setup"

    if not template_file.is_file():
        print(f"*** Error: Cannot find template file {template_file}")
        sys.exit(1)

    generate_setup(template_file, setup_file)
