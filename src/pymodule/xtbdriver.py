import os

from .veloxchemlib import XTBDriver


def _XTBDriver_compute(self, molecule, ostream):
    """
    Computes DFT-B energy using XTB driver.

    :param molecule:
        The molecule.
    :param ostream:
        The output stream.
    """

    ostream.print_blank()
    ostream.print_header('XTB Driver')
    ostream.print_header(12 * '=')
    ostream.flush()

    self._compute(molecule)

    if self.is_master_node():
        for line in self.get_output():
            ostream.print_line(line)
        ostream.flush()
        if os.path.isfile(self.get_output_filename()):
            os.remove(self.get_output_filename())


XTBDriver.compute = _XTBDriver_compute
