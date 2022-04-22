#!/usr/bin/env python3

"""Crawl through `src` and update `CMakeLists.txt` to include all sources."""

import re
from pathlib import Path
from textwrap import indent

UPDATE = [
    {
        # location of sources
        "folder": "src",
        # regexes to replace contents
        "regexes": [
            # first pattern:
            # target_sources(vlxobjs
            #   PRIVATE
            #   <contents-to-replace>
            #   )
            re.compile(
                r"(?P<HEAD>target_sources\(vlxobjs\s+PRIVATE\s)(.*?)(\s+\))",
                re.DOTALL | re.MULTILINE,
            ),
            # second pattern:
            # pybind11_add_module(veloxchemlib
            #   NO_EXTRAS
            #   <contents-to-replace>
            #   )
            re.compile(
                r"(?P<HEAD>pybind11_add_module\(veloxchemlib\s+NO_EXTRAS\s)(.*?)(\s+\))",
                re.DOTALL | re.MULTILINE,
            ),
        ],
        # folders to exclude
        "exclude": ["pymodule", "device"],
    },
    {
        # location of sources
        "folder": "unit_tests",
        # regexes to replace contents
        "regexes": [
            # target_sources(utests
            #   PRIVATE
            #   <contents-to-replace>
            #   )
            re.compile(
                r"(?P<HEAD>target_sources\(utests\s+PRIVATE\s)(.*?)(\s+\))",
                re.DOTALL | re.MULTILINE,
            ),
        ],
        # folders to exclude
        "exclude": [],
    },
]


def update_cmakelists():

    for up in UPDATE:
        where = (Path(__file__).parents[1] / up["folder"]).resolve()
        folders = sorted(
            (_ for _ in where.rglob("*") if _.is_dir() and _.name not in up["exclude"])
        )
        print(folders)

        for f in folders:
            cmakelists = f / "CMakeLists.txt"

            # get old contents
            with cmakelists.open("r") as fh:
                old = fh.read()

            # get sorted list of all source files
            cpps = sorted(
                (_.name for _ in f.iterdir() if _.is_file() and _.suffix == ".cpp")
            )
            # format list of sources as a string
            repl = indent("\n".join(cpps), prefix="    ")

            # generate new contents of CMakeLists.txt
            new = ""
            for p in up["regexes"]:
                if p.search(old):
                    new = p.sub(fr"\g<{'HEAD'}>{repl}\n  )", old)
                    # use only the regex that matches
                    break

            # overwrite CMakeLists.txt
            with cmakelists.open("w") as fh:
                fh.write(new)


if __name__ == "__main__":
    update_cmakelists()
