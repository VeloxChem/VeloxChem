import unittest
import sys
from pathlib import Path

from veloxchem.outputstream import OutputStream


class TestStreams(unittest.TestCase):

    def test_file(self):

        here = Path(__file__).parent
        outfile = str(here / 'inputs/dummy.out')

        ostream = OutputStream(outfile)
        ostream.print_line("")
        ostream.print_title("")
        ostream.print_separator()
        ostream.print_blank()
        self.assertTrue(ostream.get_state())

    def test_stdout(self):

        ostream = OutputStream(sys.stdout)
        ostream.print_blank()
        self.assertTrue(ostream.get_state())

    def test_none(self):

        ostream = OutputStream()
        ostream.print_line("")
        ostream.print_title("")
        ostream.print_separator()
        ostream.print_blank()
        self.assertFalse(ostream.get_state())


if __name__ == "__main__":
    unittest.main()
