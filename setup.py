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


class MakefileExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = Path(sourcedir).resolve()


class MakefileBuild(build_ext):
    def run(self):
        build_lib = Path(self.build_lib)
        if not build_lib.exists():
            build_lib.mkdir(parents=True)

        f = Path("src", "Makefile.setup")
        template = Path("config", "Setup.template")
        generate_setup(template_file=template, setup_file=f, build_lib=build_lib)

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
            f"{ext.sourcedir}/src",
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
    version="1.0rc1.post1",
    package_data={
        "veloxchem":
        # glob and sort
        [str(_) for _ in sorted(Path("basis").resolve().glob("**"))],
    },
    package_dir={
        "veloxchem": str(Path("src/pymodule").resolve()),
    },
    ext_modules=[
        MakefileExtension("src"),
    ],
    cmdclass={
        "build_ext": MakefileBuild,
    },
)
