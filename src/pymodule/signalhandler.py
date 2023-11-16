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

import signal


class SignalHandler:
    """
    Impements signal handler.

    Instance variable
        - func: The exit function.
        - args: The arguments for the exit function.
        - original_sigterm_handler: The original SIGTERM handler.
    """

    def __init__(self):
        """
        Initializes signal handler.
        """

        self.func = None
        self.args = None

        self.original_sigterm_handler = signal.getsignal(signal.SIGTERM)

    def sigterm_handler(self, signum, frame):
        """
        Handles the SIGTERM signal.

        :param signum:
            The signum parameter as in the handler function.
        :param frame:
            The frame parameter as in the handler function.
        """

        if self.func is not None:
            self.func(*self.args)

    def add_sigterm_function(self, func, *args):
        """
        Sets the handler for SIGTERM.

        :param func:
            The exit function.
        :param args:
            The arguments for the exit function.
        """

        self.func = func
        self.args = args

        signal.signal(signal.SIGTERM, self.sigterm_handler)

    def remove_sigterm_function(self):
        """
        Resets the handler for SIGTERM.
        """

        self.func = None
        self.args = None

        signal.signal(signal.SIGTERM, self.original_sigterm_handler)
