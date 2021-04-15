#!/usr/bin/env python3

"""Crawl through `src` and update `CMakeLists.txt` to include all sources."""

import re
from pathlib import Path
from textwrap import indent


def update_cmakelists():
    # location of sources
    src = (Path(__file__).parents[1] / "src").resolve()

    # regex to replace contents
    # first pattern:
    # target_sources(vlxobjs
    #   PRIVATE
    #   <contents-to-replace>
    #   )
    p1 = re.compile(
        r"(?P<srcs>target_sources\(vlxobjs\s+PRIVATE\s)(.*?)(\s+\))",
        re.DOTALL | re.MULTILINE,
    )
    # second pattern:
    # pybind11_add_module(veloxchemlib
    #   NO_EXTRAS
    #   <contents-to-replace>
    #   )
    p2 = re.compile(
        r"(?P<pb11>pybind11_add_module\(veloxchemlib\s+NO_EXTRAS\s)(.*?)(\s+\))",
        re.DOTALL | re.MULTILINE,
    )

    # all folders inside `src`, except:
    # - `pymodule`
    exclude = ["pymodule"]
    folders = sorted((_ for _ in src.iterdir() if _.name not in exclude and _.is_dir()))

    for f in folders:
        cmakelists = f / "CMakeLists.txt"
        print(f)
        # get old contents
        with cmakelists.open("r") as fh:
            old = fh.read()

        cpps = sorted(
            (_.name for _ in f.iterdir() if _.is_file() and _.suffix == ".cpp")
        )
        repl = indent("\n".join(cpps), prefix="    ")
        if p1.match(old):
            new = p1.sub(fr"\g<{'srcs'}>{repl}\n  )", old)
        else:
            new = p2.sub(fr"\g<{'pb11'}>{repl}\n  )", old)

        # overwrite CMakeLists.txt
        with cmakelists.open("w") as fh:
            fh.write(new)


if __name__ == "__main__":
    update_cmakelists()
