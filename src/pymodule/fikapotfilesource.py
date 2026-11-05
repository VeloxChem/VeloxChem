#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pathlib import Path

from .inputparser import InputParser


def normalize_potfile_source(potfile):
    """Return a validated, explicit potfile descriptor.

    A string or path is the legacy spelling for a file descriptor. Explicit
    descriptors use ``{"type": "file", "path": ...}`` or
    ``{"type": "string", "text": ...}``.
    """

    if isinstance(potfile, (str, Path)):
        path = str(potfile)
        if not path:
            raise ValueError("potfile path must not be empty")
        return {"type": "file", "path": path}

    if not isinstance(potfile, dict):
        raise TypeError("potfile must be a path or a descriptor dictionary")

    source_type = potfile.get("type")
    if source_type == "file":
        unknown = set(potfile) - {"type", "path"}
        if unknown:
            raise ValueError(
                f"file potfile descriptor has unknown entries {sorted(unknown)}"
            )
        path = potfile.get("path")
        if not isinstance(path, (str, Path)) or not str(path):
            raise ValueError("file potfile descriptor needs a non-empty path")
        return {"type": "file", "path": str(path)}

    if source_type == "string":
        unknown = set(potfile) - {"type", "text"}
        if unknown:
            raise ValueError(
                "string potfile descriptor has unknown entries "
                f"{sorted(unknown)}")
        text = potfile.get("text")
        if not isinstance(text, str) or not text.strip():
            raise ValueError("string potfile descriptor needs non-empty text")
        return {"type": "string", "text": text}

    raise ValueError("potfile type must be 'file' or 'string'")


def potfile_source_path(potfile):
    """Return the file path of a potfile source, or ``None`` if inline."""

    source = normalize_potfile_source(potfile)
    return source["path"] if source["type"] == "file" else None


def potfile_source_description(potfile):
    """Return a concise description suitable for output messages."""

    source = normalize_potfile_source(potfile)
    if source["type"] == "file":
        return source["path"]
    return "inline string"


def parse_potfile_source(potfile):
    """Parse a file or inline potfile source and return data and its label."""

    source = normalize_potfile_source(potfile)
    if source["type"] == "file":
        label = source["path"]
        data = InputParser(label).input_dict
    else:
        label = "inline potfile"
        data = InputParser.from_string(source["text"], source=label).input_dict
    return data, label
