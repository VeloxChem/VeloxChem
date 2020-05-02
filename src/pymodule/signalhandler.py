import signal
import sys
import os


class SignalHandler:

    def __init__(self):

        self.func = None
        self.args = None

        self.terminated = False
        self.original_sigterm_handler = signal.getsignal(signal.SIGTERM)

    def sigterm_handler(self, signum, frame):

        if self.terminated:
            return

        if not self.terminated:
            self.terminated = True

        if self.func is not None:
            if self.args is not None:
                sys.exit(self.func(*self.args))
            else:
                sys.exit(self.func())
        else:
            sys.exit(1)

    def raise_signal(self, sigstr):

        if sigstr == 'SIGTERM':
            # available in Python 3.8
            # signal.raise_signal(signal.SIGTERM)
            os.kill(os.getpid(), signal.SIGTERM)

    def add_sigterm_function(self, func, *args):

        self.func = func
        self.args = args

        signal.signal(signal.SIGTERM, self.sigterm_handler)

    def remove_sigterm_function(self):

        self.func = None
        self.args = None

        signal.signal(signal.SIGTERM, self.original_sigterm_handler)
