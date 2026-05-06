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

import numpy as np
import re
from typing import List, Sequence, Any, Tuple, Union


def nn(s: str) -> str:
    """Remove all digits from a string.

    Args:
        s (str): Input string.

    Returns:
        str: String with all digits removed.

    Example:
        >>> nn("Fe123")
        'Fe'
    """
    return re.sub(r"\d+", "", s)


def nl(s: str) -> str:
    """Remove all non-digit characters from a string.

    Args:
        s (str): Input string.

    Returns:
        str: String containing only digits.

    Example:
        >>> nl("ABC123")
        '123'
    """
    return re.sub(r"\D+", "", s)


def pname(s: str) -> str:
    """Extract the prefix before the first underscore.

    Args:
        s (str): Input string of the form 'prefix_something'.

    Returns:
        str: Substring before the first underscore.

    Example:
        >>> pname("Cu_1 (0.0 2.2 3.3)")
        'Cu'
    """
    return s.split("_")[0]


def lname(s: str) -> np.ndarray:
    """Extract coordinates embedded in a string after the first underscore.

    Args:
        s (str): String with format 'atom_(x y z)'.

    Returns:
        np.ndarray: Array of floats with coordinates. If not present, returns array [0.0, 0.0, 0.0].

    Example:
        >>> lname("Cu_(1.0 2.0 3.0)")
        array([1., 2., 3.])
        >>> lname("Cu")
        array([0., 0., 0.])
    """
    if len(s.split("_")) < 2:
        lis = np.array([0.0, 0.0, 0.0])
    else:
        lis = np.asanyarray(s.split("_")[1][1:-1].split(), dtype=float)
    return lis


def arr_dimension(arr: np.ndarray) -> int:
    """Determine if an array is 1D or 2D.

    Args:
        arr (np.ndarray): Input numpy array.

    Returns:
        int: 2 if arr.ndim > 1, else 1.

    Example:
        >>> arr_dimension(np.ones((4, 3)))
        2
        >>> arr_dimension(np.ones(3))
        1
    """
    if arr.ndim > 1:
        return 2
    else:
        return 1


def is_list_A_in_B(A: Sequence[Any], B: Sequence[Any]) -> bool:
    """Check if all elements of A are equal (allclose) to elements of B (pairwise).

    Args:
        A (Sequence[Any]): First list/array.
        B (Sequence[Any]): Second list/array.

    Returns:
        bool: True if all pairs are close within tol, else False.

    Example:
        >>> is_list_A_in_B([np.array([1.0])], [np.array([1.0])])
        True
    """
    return all([np.allclose(a, b, atol=1e-9) for a, b in zip(A, B)])


def remove_blank_space(value: str) -> str:
    """Remove all whitespace from a string.

    Args:
        value (str): Input string.

    Returns:
        str: String with all whitespace removed.

    Example:
        >>> remove_blank_space("a b c")
        'abc'
    """
    return re.sub(r"\s", "", value)


def remove_empty_lines(lines: List[str]) -> List[str]:
    """Remove empty or whitespace-only lines from a list of strings.

    Args:
        lines (List[str]): Input list of lines.

    Returns:
        List[str]: New list without empty lines.

    Example:
        >>> remove_empty_lines(["a", "   ", "b"])
        ['a', 'b']
    """
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip() != "":
            newlines.append(lines[i])
    return newlines


def remove_bracket(value: str) -> float:
    """Remove anything inside parentheses from a string and convert to float.

    Args:
        value (str): Input string with numeric value and optional bracket.

    Returns:
        float: Value with parenthetical part removed.

    Example:
        >>> remove_bracket("42(3)")
        42.0
    """
    value_float = float(re.sub(r"\(.*?\)", "", value))
    return value_float


def remove_tail_number(value: str) -> str:
    """Remove all digits from a string.

    Args:
        value (str): String input.

    Returns:
        str: String with digits removed.

    Example:
        >>> remove_tail_number("AB12")
        'AB'
    """
    return re.sub(r"\d", "", value)


def add_quotes(value: str) -> str:
    """Add single quotes around a string.

    Args:
        value (str): Input string.

    Returns:
        str: String wrapped in single quotes.

    Example:
        >>> add_quotes("Fe")
        "'Fe'"
    """
    return "'" + value + "'"


def remove_note_lines(lines: List[str]) -> List[str]:
    """Remove lines containing an underscore '_' from a list of strings.

    Args:
        lines (List[str]): List of input lines.

    Returns:
        List[str]: List of lines with lines containing '_' removed.

    Example:
        >>> remove_note_lines(["foo", "bar_baz"])
        ['foo']
    """
    newlines = []
    for i in range(len(lines)):
        m = re.search(r"_", lines[i])
        if m is None:
            newlines.append(lines[i])
    return newlines


def extract_quote_lines(lines: List[str]) -> List[str]:
    """Extract lines that start with a single quote.

    Args:
        lines (List[str]): List of input lines.

    Returns:
        List[str]: List of lines whose first non-whitespace character is a single quote.

    Example:
        >>> extract_quote_lines(["'a'", "b"])
        ["'a'"]
    """
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip()[0] == "'":
            newlines.append(lines[i])
    return newlines


def extract_xyz_lines(lines: List[str]) -> List[str]:
    """Process lines: remove whitespace, add quotes, and select lines not starting with '_'.

    Args:
        lines (List[str]): List of input lines.

    Returns:
        List[str]: Processed and filtered lines.

    Note:
        Lines starting with '_' are ignored. Whitespace removed and line quoted.

    Example:
        >>> extract_xyz_lines(["abc", "_not", "def"])
        ["'abc'", "'def'"]
    """
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip()[0] != "_":
            quote_value = add_quotes(remove_blank_space(lines[i]).strip())
            newlines.append(quote_value)
    newlines = remove_empty_lines(newlines)
    return newlines


def remove_quotes(value: str) -> str:
    """Extract value within single or double quotes from a string.

    Args:
        value (str): Input string containing a quoted value.

    Returns:
        str: Extracted string inside quotes.

    Example:
        >>> remove_quotes("'Fe'")
        'Fe'
    """
    pattern = r"[\"']([^\"']+)[\"']"
    extracted_values = re.findall(pattern, value)
    return extracted_values[0]


def convert_fraction_to_decimal(expression: str) -> str:
    """Convert all fractions in a string to their decimal form.

    Args:
        expression (str): String containing fractions, e.g., '3/4'.

    Returns:
        str: Expression with fractions converted to decimals.

    Example:
        >>> convert_fraction_to_decimal("a 3/4 b")
        'a 0.75 b'
    """

    def replace_fraction(match):
        numerator, denominator = map(int, match.groups())
        return str(numerator / denominator)

    fraction_pattern = r"(-?\d+)/(\d+)"
    converted_expression = re.sub(fraction_pattern, replace_fraction,
                                  expression)
    return converted_expression


def find_keyword(keyword: str, s: str) -> bool:
    """Check whether a keyword is present in a string.

    Args:
        keyword (str): Keyword or regex pattern to search for.
        s (str): String to search in.

    Returns:
        bool: True if keyword found, False otherwise.

    Example:
        >>> find_keyword("foo", "barfoo")
        True
    """
    m = re.search(keyword, s)
    if m:
        return True
    else:
        return False


def locate_min_idx(matrix: np.ndarray) -> Tuple[int, int]:
    """Find the indices of the minimum value of a 2D numpy array.

    Args:
        matrix (np.ndarray): Input 2D array.

    Returns:
        Tuple[int, int]: Row and column index of the minimum value.

    Example:
        >>> mat = np.array([[1,2],[3,0]])
        >>> locate_min_idx(mat)
        (1, 1)
    """
    min_idx = np.unravel_index(matrix.argmin(), matrix.shape)
    return min_idx[0], min_idx[1]
