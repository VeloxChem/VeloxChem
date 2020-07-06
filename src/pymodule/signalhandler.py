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
