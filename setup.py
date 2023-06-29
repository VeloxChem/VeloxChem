#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
from os import environ, cpu_count

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

# number of cores to use for compiling
# if unspecified, use all available cores to compile
ncores = cpu_count()
if ncores is None:
    ncores = 1
compile_jobs = int(environ.get("VLX_NUM_BUILD_JOBS", default=ncores))

# we keep most of the configuration in setup.cfg
# usage of setup.py is deprecated
# https://setuptools.pypa.io/en/latest/userguide/quickstart.html#setup-py
setup(
    cmake_args=[
        "-DCMAKE_JOB_POOL_COMPILE:STRING=compile",
        "-DCMAKE_JOB_POOL_LINK:STRING=link",
        # always use 2 cores for linking
        f"-DCMAKE_JOB_POOLS:STRING=compile={compile_jobs:d};link=2",
    ],
)
