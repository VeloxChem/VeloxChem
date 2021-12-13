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

import sys
from multiprocessing import cpu_count
from os import environ
from pathlib import Path
from subprocess import STDOUT, check_call

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Newer packaging standards may recommend removing the current dir from the
# path, add it back if needed.
if "" not in sys.path:
    sys.path.insert(0, "")
from config.generate_setup import generate_setup


def get_version(rel_path):
    abs_path = Path(__file__).parent / rel_path
    with abs_path.open("r") as fh:
        lines = fh.readlines()
    for line in lines:
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


class MakefileExtension(Extension):

    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = Path(sourcedir).resolve()


class MakefileBuild(build_ext):

    def run(self):
        build_lib = Path(self.build_lib)
        if not build_lib.exists():
            build_lib.mkdir(parents=True)

        setup_file = Path("src", "Makefile.setup")
        template_file = Path("config", "Setup.template")
        generate_setup(template_file=template_file,
                       setup_file=setup_file,
                       build_lib=build_lib)

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        if hasattr(self, "parallel") and self.parallel:
            njobs = self.parallel
        else:
            njobs = environ.get("VLX_NUM_BUILD_JOBS", default=cpu_count())

        build_cmd = [
            "make",
            "-C",
            str(ext.sourcedir / "src"),
            "release",
            "-s",
            f"-j{int(njobs):d}",
        ]

        check_call(
            build_cmd,
            cwd=self.build_lib,
            stdout=sys.stdout,
            stderr=STDOUT,
        )


setup(
    version=get_version(Path("src", "pymodule", "__init__.py")),
    package_data={
        # glob and sort
        "veloxchem": [
            str(p) for p in sorted(Path("basis").resolve().glob("**"))
        ],
    },
    package_dir={
        "veloxchem": str(Path("src", "pymodule").resolve()),
    },
    ext_modules=[
        MakefileExtension("src"),
    ],
    cmdclass={
        "build_ext": MakefileBuild,
    },
)
