import collections
from collections import abc

_Pos = collections.namedtuple("Position", ("start", "stop"))
_Index = collections.namedtuple("Index", ("start", "stop"))

__slots__ = ("ProteinLocation", )

class ProteinLocation():
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    # protein:      ELVISLIVES
    #  - positions: 1234567890
    # peptide:          -----
    #  - positions      5   9

    seq = "ELVISLIVES"
    peptide = ProteinLocation(5, 9, seq)
    peptide.pos  # returns (5, 9)
    peptide.slice  # slice(4, 9, None)
    seq[peptide.slice]  # "LIVE"
    seq[peptide.index.start]  # "L"
    seq[peptide.index.stop]  # "E"
    peptide2 = ProteinLocaion(1, 5)
    peptide2 < peptide 1  # True, because (1, 5) < (5, 9)
    """

    def __init__(self, start, stop=None, seq=None, protein_sequence=None):
        """peptide cordinate, counting from 1, with inclusive end"""
        if stop is None:
            stop = start

        self.pos = _Pos(int(start), int(stop)) 
        if self.pos.start < 1:
            raise ValueError("start < 1")
        if self.pos.stop < self.pos.start:
            raise ValueError("stop < start")

        self.index = _Index(self.pos.start - 1, self.pos.stop - 1)
        self.slice = slice(self.pos.start - 1, self.pos.stop)

        if seq:
            self.seq = seq
        elif protein_sequence:
            self.seq = protein_sequence[self.slice]

        self.length = self.pos.stop - self.pos.start + 1

    # alternative constructors
    @classmethod
    def from_index_and_length(cls, start, length):
        start += 1
        end = length + start - 1
        return cls(start, end)

    @classmethod
    def from_indexes(cls, start, end):
        return cls(start + 1, end + 1)

    @classmethod
    def from_slice(cls, start, stop=None):
        if isinstance(start, slice):
            return cls(start.start + 1, start.stop)
        return cls(start + 1, stop)

    @classmethod
    def from_sequence(cls, protein_sequence, peptide_sequence):
        start_index = protein_sequence.find(peptide_sequence)
        if start_index == -1:
            raise IndexError("{} not in {}".format(peptide_sequence, protein_sequence))

        start = start_index + 1
        end = start_index + len(peptide_sequence)
        return cls(start, end, seq=peptide_sequence, protein_sequence=protein_sequence)

    # dunders
    def __len__(self):
        return self.length

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, *self.pos)

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.pos < other.pos
        elif isinstance(other, abc.Sized) and isinstance(other, abc.Iterable):
            if len(other) == 2:
                return self.pos < other[:2]
        return False

    def __gt__(self, other):
        return not self < other

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pos == other.pos
        elif isinstance(other, abc.Sized) and isinstance(other, abc.Iterable):
            if len(other) == 2:
                return self.pos == other[:2]
        return False

    def __hash__(self):
        return hash(repr(self))
