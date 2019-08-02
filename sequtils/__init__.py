# core imports
import collections
from collections.abc import Sequence
from abc import ABCMeta, abstractmethod
import operator
import functools

__slots__ = ("SequencePoint", "SequenceRange")


# named tuples used by SequenceRange
_Pos = collections.namedtuple("Position", ("start", "stop"))
_Index = collections.namedtuple("Index", ("start", "stop"))


@functools.total_ordering
class BaseSequenceLocation:
    __metaclass__ = ABCMeta

    def __new__(cls, arg, *args, **kwargs):
        if isinstance(arg, cls):
            # if arg is already of the correct type, then keep it because subclasses are mutable
            # just like other immutables: x = 213124512421312; x is int(x)
            return arg
        return super().__new__(cls)

    # read only attributes
    @property
    def pos(self):
        return self._pos

    @property
    def index(self):
        return self._index

    @property
    def slice(self):
        return self._slice

    def is_valid(self):
        try:
            self.validate()
            return True
        except ValueError:
            return False

    def __lt__(self, other):
        try:
            return self.pos < self.__class__(other, validate=False).pos
        except ValueError:
            raise TypeError("cannot compare type {} and {}".format(type(self), type(other)))

    def __eq__(self, other):
        try:
            return self.pos == self.__class__(other, validate=False).pos
        except ValueError:
            return NotImplemented

    #  # using only __lt__ and __eq__
    #  def __le__(self, other):
        #  return self < other or self == other
#
    #  # using only __lt__ and __eq__
    #  def __ge__(self, other):
        #  return other < self or self == other
#
    #  def __gt__(self, other):
        #  return not self <= other

    def __hash__(self):
        return hash(self.pos)

    @abstractmethod
    def _join(self, other, operator):
        "helper methood, needed to make __add__ and __sub__ work"
        raise NotImplementedError("Please Implement this method")

    def _arithmetic(self, other, operator):
        """
        helper method that calls _join, ensures that NotImplemented is returned
        when a cast fails, this is important to ensure that stuff like the
        following woks:
            SequencePoint(2) + SequenceRecord(3, 5)
        this works because:
            SequencePoint.__add__ calls
                SequencePoint.__init__ who raises an ValueError, thus
            SequencePoint_arthmetic returns NotImplemented,
                This signals to the Python Interperter to use
            SequenceRange.__radd__ who calls
                SequenceRange.__init__, which does not fail :)
                    because you can convert Points to Ranges but not vice versa!
        """

        if not isinstance(other, self.__class__):
            try:
                other = self.__class__.from_index(other, validate=False)
            except TypeError:
                return NotImplemented
        return self._join(other, operator)

    def __add__(self, other):
        return self._arithmetic(other, operator.add)

    def __sub__(self, other):
        return self._arithmetic(other, operator.sub)

    def __radd__(self, other):
        return self._arithmetic(other, operator.add)
        #  return self.__class__.from_index(other, validate=False) + self

    def __rsub__(self, other):
        # other - self
        return self.__class__.from_index(other, validate=False) - self


class SequencePoint(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    # protein:      ELVISLIVES
    #  - positions: 1234567890
    # peptide:           LIVE
    #  - positions       6  9

    seq = "ELVISLIVES"
    mutation = ProteinLocation(6)
    mutation.pos  # returns 6
    seq[mutation.index]  # "L"
    """

    def __init__(self, position, *, validate=True):
        if isinstance(position, self.__class__.mro()[1]):  # isinstance of parent
            if isinstance(position, self.__class__):
                return
            if isinstance(position, SequenceRange):
                if len(position) != 1:
                    raise TypeError("can only Convert {} to {} if len({}) = 1".format(
                        type(position), type(self), type(position)))
                position = position.pos.start

        self._pos = int(position)
        if validate:
            self.validate()
        self._index = self._pos - 1
        self._slice = slice(self.index, self.index + 1)

    def validate(self):
        if self.pos < 1:
            raise ValueError("position < 1")

    # alternative constructors
    @classmethod
    def from_index(cls, index, *, validate=True):
        return cls(index + 1, validate=validate)

    # implementation of abstract methods
    def _join(self, other, operator):
        return self.__class__.from_index(operator(self.index, other.index), validate=False)

    # other dunders
    def __str__(self):
        return str(self.pos)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.pos)

    #  def iter_pos(self):
        #  yield self.pos


class SequenceRange(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    # protein:      ELVISLIVES
    #  - positions: 1234567890
    # peptide:           -----
    #  - positions       6  9

    seq = "ELVISLIVES"
    peptide = SequenceRange(6, 9, seq)
    peptide.pos  # returns (6, 9)
    peptide.slice  # slice(5, 9, None)
    seq[peptide.slice]  # "LIVE"
    seq[peptide.index.start]  # "L"
    seq[peptide.index.stop]  # "E"
    peptide2 = SequenceRange(1, 5)
    peptide2 < peptide 1  # True, because (1, 5) < (5, 9)
    """

    _str_separator = ':'

    def __init__(self, start, stop=None, seq=None, protein_sequence=None, *, validate=True):
        """
        sequence cordinate, counting from 1, with inclusive stop
        viable calls to the constructor includes:
        SequenceRange(1, 2)
        SequenceRange(SequencePoint(1), SequencePoint(2))
        SequenceRange((1, 2))
        SequenceRange(SequenceRange(1, 2))
        """

        if isinstance(start, self.__class__.mro()[1]):  # isinstance of parent
            if isinstance(start, self.__class__):
                return
            elif isinstance(start, SequencePoint):
                start = start.pos

        if isinstance(stop, SequencePoint):
            stop = stop.pos
        elif stop is None:
            if isinstance(start, Sequence) and len(start) == 2:
                start, stop = start[:2]
            elif seq is not None:
                stop = self._stop_from_start_and_length(start, len(seq))
            else:
                stop = start

        self._pos = _Pos(int(start), int(stop))
        if validate:
            self.validate()

        self._index = _Index(self.pos.start - 1, self.pos.stop - 1)
        self._slice = slice(self.pos.start - 1, self.pos.stop)

        self._seq = self._get_seq(seq, protein_sequence)
        self.length = self.pos.stop - self.pos.start + 1

    # __init__ helper methods
    def _get_seq(self, seq, protein_sequence):
        if seq:
            return str(seq)
        elif protein_sequence:
            return protein_sequence[self.slice]
        return None

    def validate(self):
        if self.pos.start < 1:
            raise ValueError("start < 1")
        if self.pos.stop < self.pos.start:
            raise ValueError("stop({}) < start({})".format(self.pos.stop, self.pos.start))

    @classmethod
    def from_index(cls, start_index, stop_index=None, *, validate=True):
        if isinstance(start_index, cls.mro()[1]):  # instance of parent
            return cls(start_index)
        if isinstance(start_index, Sequence) and len(start_index) == 2:
            start_index, stop_index = start_index
        if stop_index is None:
            stop_index = start_index
        return cls(start_index + 1, stop_index + 1, validate=validate)

    @classmethod
    def _stop_from_start_and_length(cls, start, length):
        return start + length - 1

    # alternative constructors
    @classmethod
    def from_index_and_length(cls, start_index, length):
        stop = cls._stop_from_start_and_length(start_index + 1, length)
        return cls(start_index + 1, stop)

    @classmethod
    def from_slice(cls, start_slice, stop_slice=None):
        if isinstance(start_slice, slice):
            return cls(start_slice.start + 1, start_slice.stop)
        return cls(start_slice + 1, stop_slice)

    @classmethod
    def from_sequence(cls, protein_sequence, peptide_sequence):
        start_index = protein_sequence.find(peptide_sequence)
        if start_index == -1:
            raise IndexError("{} not in {}".format(peptide_sequence, protein_sequence))

        start = start_index + 1
        stop = start_index + len(peptide_sequence)
        return cls(start, stop, seq=peptide_sequence, protein_sequence=protein_sequence)

    # implementation of abstract methods
    def _join(self, other, operator):
        return self.__class__.from_index(operator(self.index.start, other.index.start),
                                         operator(self.index.stop, other.index.stop),
                                         validate=False)

    # other dunders
    def __len__(self):
        return self.length

    def __str__(self):
        if self.pos.start == self.pos.stop:
            return str(self.pos.start)
        return "{}{}{}".format(self.pos.start, self._str_separator, self.pos.stop)

    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, *self.pos)

    #  def __iter__(self):
        #  yield from self.iter_pos()

    #  def iter_pos(self):
        #  yield from self.pos

    def __iter__(self):
        for pos in range(self.pos.start, self.pos.stop + 1):
            yield SequencePoint(pos, validate=False)

    def __getitem__(self, item):
        return self.pos[item]

    def _contains(self, sequence_point):
        "Helper method that checks if a SequencePoint is in self"

        return self.pos.start <= sequence_point.pos <= self.pos.stop

    def __contains__(self, item):
        " returns True if all of item is inside self"
        return self.contains(item, part=all)

    def contains(self, item, part=all):
        """
        Check wheter item is inside 'self'.
        by default (part=all) all amino acids has to be inside self:
        by changing to (part=any), then only one of the amino acids have to be inside:
        self = -----ELVISLIVES
        item = ----------L----        <--- part=all -> True,  part=any -> True
        item = ----------LIVE-        <--- part=all -> True,  part=any -> True
        item = ----------LIVESANDDIES <--- part=all -> False, part=any -> True
        item = ELVENELVISLIVESANDDIES <--- part=all -> False, part=any -> True
        """

        if isinstance(item, self.__class__.mro()[0]):
            return part(map(self._contains, item))
        try:
            return self._contains(SequencePoint(item))
        except (ValueError, TypeError):
            try:
                return part(map(self._contains, SequenceRange(item)))
            except ValueError:
                return False
        return False

    # properties, to make it read-only
    @property
    def seq(self):
        return self._seq
