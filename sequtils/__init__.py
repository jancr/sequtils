
# core imports
import collections
from collections import abc
from abc import ABCMeta, abstractmethod

__slots__ = ("SequencePoint", "SequenceRange")


# named tuples used by SequenceRange
_Pos = collections.namedtuple("Position", ("start", "stop"))
_Index = collections.namedtuple("Index", ("start", "stop"))


class BaseSequenceLocation:
    __metaclass__ = ABCMeta

    @abstractmethod
    def __lt__(self, other):
        raise NotImplementedError("Please Implement this method")

    def __gt__(self, other):
        return not self < other

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError("Please Implement this method")

class SequencePoint(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    # protein:      ELVISLIVES
    #  - positions: 1234567890
    # peptide:          -----
    #  - positions      5   9

    seq = "ELVISLIVES"
    mutation = ProteinLocation(5)
    mutation.pos  # returns 5
    seq[mutation.index]  # "L"
    """

    def __init__(self, position, *, validate=True):
        if isinstance(position, SequenceRange):
            if len(position) != 1:
                raise("can only Convert {} to {} if len({}) = 1".format(
                    type(position), type(self), type(position)))
            position = position.pos.start

        self.pos = int(position)
        if self.pos < 1:
            raise ValueError("position < 1")
        self.index = self.pos - 1

    # alternative constructors
    @classmethod
    def from_index(cls, index):
        return cls(index + 1)

    # dunders
    def __str__(self):
        return str(self.pos)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.pos)

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.pos < other.pos
        try:
            other = self.__class__(other)
            return self.pos < other.pos
        except:
            raise TypeError("cannot compare type {} and {}".format(type(self), type(other)))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pos == other.pos
        try:
            other = self.__class__(other, validate=False)
            return self.pos == other.pos
        except:
            return False

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            other = self.__class__(other, validate=False)
        return self.__class__(self.pos + other.pos)

    def __sub__(self, other):
        if not isinstance(other, self.__class__):
            other = self.__class__(other, validate=False)
        return self.__class__(self.pos - other.pos)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return self.__class__(other, validate=False) - self


class SequenceRange(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    # protein:      ELVISLIVES
    #  - positions: 1234567890
    # peptide:          -----
    #  - positions      5   9

    seq = "ELVISLIVES"
    peptide = SequenceRange(5, 9, seq)
    peptide.pos  # returns (5, 9)
    peptide.slice  # slice(4, 9, None)
    seq[peptide.slice]  # "LIVE"
    seq[peptide.index.start]  # "L"
    seq[peptide.index.stop]  # "E"
    peptide2 = SequenceRange(1, 5)
    peptide2 < peptide 1  # True, because (1, 5) < (5, 9)
    """

    _str_separator = ':'

    def __init__(self, start, stop=None, seq=None, protein_sequence=None, *, validate=True):
        """sequence cordinate, counting from 1, with inclusive stop"""

        if isinstance(start, SequencePoint):
            start = start.pos

        if stop is None:
            stop = start
        elif isinstance(stop, SequencePoint):
            stop = stop.pos

        self.pos = _Pos(int(start), int(stop)) 
        if validate:
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
        stop = length + start - 1
        return cls(start, stop)

    @classmethod
    def from_indexes(cls, start, stop):
        return cls(start + 1, stop + 1)

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
        stop = start_index + len(peptide_sequence)
        return cls(start, stop, seq=peptide_sequence, protein_sequence=protein_sequence)

    # dunders
    def __len__(self):
        return self.length

    def __str__(self):
        return "{}{}{}".format(self.pos.start, self._str_separator, self.pos.stop)

    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, *self.pos)

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.pos < other.pos
        elif isinstance(other, abc.Sized) and isinstance(other, abc.Iterable):
            if len(other) == 2:
                return self.pos < other[:2]
        try:
            return self < self.__class__(other)
        except:
            raise TypeError("cannot compare types {} and {}".format(type(self), type(other)))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pos == other.pos
        elif isinstance(other, abc.Sized) and isinstance(other, abc.Iterable):
            if len(other) == 2:
                return self.pos == other[:2]
        try:
            return self == self.__class__(other)
        except:
            return False

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            other = self.__class__(other, validate=False)
        return self.__class__(self.pos.start + other.pos.start, 
                                self.pos.stop + other.pos.stop)

    def __sub__(self, other):
        if not isinstance(other, self.__class__):
            other = self.__class__(other, validate=False)
        return self.__class__(self.pos.start - other.pos.start, 
                                self.pos.stop - other.pos.stop)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return self.__class__(other, validate=False) - self

    def __hash__(self):
        return hash(repr(self))
