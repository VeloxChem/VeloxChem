from pathlib import Path
import sys

from veloxchem.outputstream import OutputStream


class TestStreams:

    def run_ostream_test(self, ostream, state):

        ostream.print_line('')
        ostream.print_title('')
        ostream.print_separator()
        ostream.print_blank()
        assert ostream.get_state() is state

    def test_file(self):

        here = Path(__file__).parent
        outfile = here / 'inputs' / 'dummy.out'

        ostream = OutputStream(str(outfile))
        self.run_ostream_test(ostream, True)

        ostream.close()
        assert not ostream.get_state()

    def test_path(self):

        here = Path(__file__).parent
        outfile = here / 'inputs' / 'dummy.out'

        ostream = OutputStream(outfile)
        self.run_ostream_test(ostream, True)

        ostream.close()
        assert not ostream.get_state()

    def test_stdout(self):

        ostream = OutputStream(sys.stdout)
        self.run_ostream_test(ostream, True)

    def test_default(self):

        ostream = OutputStream()
        self.run_ostream_test(ostream, True)

    def test_none(self):

        ostream = OutputStream(None)
        self.run_ostream_test(ostream, False)
