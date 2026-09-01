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
from typing import List, Sequence, Any


def fetch_pdbfile(
    dir_name: str,
    keywords: Sequence[str],
    nokeywords: Sequence[str],
    ostream: Any,
) -> List[str]:
    """Finds PDB files in a directory matching keyword filters.

    Recursively searches for `.pdb` files under the specified directory whose filenames
    contain all substrings in `keywords` and do not contain any substring from `nokeywords`.

    Args:
        dir_name (str): The directory path to search within.
        keywords (Sequence[str]): List of substrings that must be present in the filename.
        nokeywords (Sequence[str]): List of substrings that must NOT be present in the filename.
        ostream (Any): Output stream object with `print_info` method for logging information.

    Returns:
        List[str]: List containing the name(s) of matching PDB file(s).

    Raises:
        ValueError: If no matching PDB file is found.

    Example:
        >>> fetch_pdbfile('mydir', ['ABC'], ['bad'], ostream)
        ["ABC1.pdb"]

    Note:
        If multiple files match, all are returned. If one match, a single-element list is returned.
    """
    candidates: List[str] = []
    for pdb in Path(dir_name).rglob("*.pdb"):
        name = pdb.name
        if all(i in name for i in keywords) and all(j not in name for j in nokeywords):
            candidates.append(pdb.name)

    if len(candidates) == 0:
        raise ValueError(f"Cannot find a file including '{keywords}'")
    elif len(candidates) == 1:
        ostream.print_info(f"Found the file including {keywords}: {candidates[0]}")
        return candidates
    else:  # len(candidates) > 1
        ostream.print_info(f"Found many files including {keywords}: {candidates}")
        return candidates
