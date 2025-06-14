# -*- coding: utf-8 -*-


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


program_header = "= Building the [" + bcolors.OKBLUE + "biokg" + bcolors.ENDC + "] knowledge graph"
dwn_sym = "(" + bcolors.FAIL + "⤓" + bcolors.ENDC + ") "
done_sym = " (" + bcolors.OKGREEN + "✓" + bcolors.ENDC + ")"
fail_sym = " (" + bcolors.FAIL + "✘" + bcolors.ENDC + ")"
prc_sym = "(" + bcolors.OKBLUE + "⏣" + bcolors.ENDC + ") "
hsh_sym = " (" + bcolors.OKBLUE + "#" + bcolors.ENDC + ") "
inf_sym = "(" + bcolors.WARNING + bcolors.BOLD + "‣" + bcolors.ENDC + ") "


def print_line():
    """ print stdout line
    """
    print("------------------------------------------------")


def print_bold_line():
    """ print stdout line
    """
    print("================================================")


def print_section_header(header_txt):
    """

    Parameters
    ----------
    header_txt: str
        header text to print
    """
    print(">>>  %s ... " % header_txt)
    print_line()

class SetWriter:
    """
    Utility class for writing DrugBank statements
    Enforces uniqueness of statements between written between flushes
    Set clear_on_flush to false to enforce uniquness for on all writes
    (should not be set for very large datasets)
    """
    def __init__(self, path):
        """
        Initialize a new SetWriter

        Parameters
        ----------
        """
        self._lines = []
        self._lineset = set()
        self._fd = open(path, 'w', encoding='utf-8')
        self._clear_on_flush = True
        self._closed = False

    @property
    def clear_on_flush(self):
        return self._clear_on_flush

    @clear_on_flush.setter
    def clear_on_flush(self, value):
        self._clear_on_flush = value 

    def write(self, line):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        
        if line in self._lineset:
            return
        self._lineset.add(line)
        self._lines.append(line)

    def flush(self):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        self._fd.writelines(self._lines)
        self._lines = []
        if self._clear_on_flush:
            self._lineset = set()

    def close(self):
        if len(self._lines) > 0:
            self.flush()
        self._lineset = set()
        self._fd.close()
