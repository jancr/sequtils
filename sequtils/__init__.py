# core imports
import collections as _collections
from collections.abc import Sequence as _Sequence
from abc import ABCMeta as _ABCMeta
from abc import abstractmethod as _abstractmethod
import operator as _operator
from functools import total_ordering as _total_ordering
import math as _math
from warnings import warn as _warn


__slots__ = ("SequencePoint", "SequenceRange")
__all__ = ("SequencePoint", "SequenceRange")


# named tuples used by SequenceRange
#  _Pos = _collections.namedtuple("Position", ("start", "stop"))
#  _Index = _collections.namedtuple("Index", ("start", "stop"))


class _Positions(_collections.namedtuple("Pos", ("_1", "_2"), rename=True)):
    _warning = ("SequenceRange.{name}.{position} is deprecated, "
                "use SequenceRange.{position}.{name} instead")

    @property
    def start(self):
        return self._warn_get('start', self._0)

    @property
    def stop(self):
        return self._warn_get('stop', self._1)

    @classmethod
    def _warn_get(cls, position, field):
        msg = cls._warning.format(name=cls.name, position=position)
        _warn(msg, DeprecationWarning)
        return field


class _Pos(_Positions):
    name = "pos"


class _Index(_Positions):
    name = "index"


@_total_ordering
class BaseSequenceLocation:
    __metaclass__ = _ABCMeta

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

    def _comparison_cast(self, other):
        """
        method used by __eq__ and __lt__ that returns True if you want it to perform implicit
        casting...
        SequencePoint(2) == SequencePoint("2")
        but
        SequencePoint(2) == "2"
        only if SequencePoint._comparison_cast(str) returns True
        """
        return isinstance(other, (self.__class__, int))

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

    # abstract methods
    @_abstractmethod
    def _join(self, other, operator):
        "helper methood, needed to make __add__ and __sub__ work"
        raise NotImplementedError("Please Implement this method")

    @_abstractmethod
    def validate(self):
        raise NotImplementedError("Please Implement this method")

    #  math dunders
    def __add__(self, other):
        return self._arithmetic(other, _operator.add)

    def __sub__(self, other):
        return self._arithmetic(other, _operator.sub)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        # other - self
        return self.__class__.from_index(other, validate=False) - self

    # comparison dunders
    def __lt__(self, other):
        if self._comparison_cast(other):
            try:
                return self.pos < self.__class__(other, validate=False).pos
            except ValueError:
                raise TypeError("cannot compare type {} and {}".format(type(self), type(other)))
        return NotImplemented

    def __eq__(self, other):
        if self._comparison_cast(other):
            try:
                return self.pos == self.__class__(other, validate=False).pos
            except (ValueError, TypeError):
                pass
        return NotImplemented

    # other dunders
    def __hash__(self):
        return hash(self.pos)


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
                position = position.start

        self._pos = int(position)
        if validate:
            self.validate()
        self._index = self._pos - 1
        self._slice = slice(self.index, self.index + 1)

    # alternative constructors
    @classmethod
    def from_index(cls, index, *, validate=True):
        return cls(index + 1, validate=validate)

    # implementation of abstract methods
    def validate(self):
        if self.pos < 1:
            raise ValueError("position({}) < 1".format(self.pos))

    def _join(self, other, operator):
        return self.__class__.from_index(operator(self.index, other.index), validate=False)

    # other dunders
    def __str__(self):
        return str(self.pos)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.pos)


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

    def __init__(self, start, stop=None, seq=None, full_sequence=None, *, validate=True,
                 length=None):
        """
        sequence cordinate, counting from 1, with inclusive stop
        viable calls to the constructor includes:
        SequenceRange(1, 2)
        SequenceRange(SequencePoint(1), SequencePoint(2))
        SequenceRange((1, 2))
        SequenceRange(SequenceRange(1, 2))

        A "start" has to be always provided, but stop can be infered from the other arguments
        thus the following creates the same object:
            SequenceRange(1, 3)
            SequenceRange((1, 3))
            SequenceRange(1, seq="ABC")
            SequenceRange(1, length=3)
        if stop is None
        """

        if isinstance(start, self.__class__.mro()[1]):  # isinstance of parent
            if isinstance(start, self.__class__):
                return
            elif isinstance(start, SequencePoint):
                start = start

        if isinstance(stop, SequencePoint):
            stop = stop
        elif stop is None:
            if isinstance(start, _Sequence) and len(start) == 2:
                start, stop = start[:2]
            elif length is not None:
                stop = self._stop_from_start_and_length(start, length)
            elif seq is not None:
                stop = self._stop_from_start_and_length(start, len(seq))
            else:
                stop = start

        self._start = SequencePoint(start, validate=validate)
        self._stop = SequencePoint(stop, validate=validate)

        if validate:
            self.validate()

        self._slice = slice(self.start.pos - 1, self.stop.pos)
        self._seq = self._get_seq(seq, full_sequence)
        self.length = self.stop.pos - self.start.pos + 1

    # alternate constructors
    @classmethod
    def from_index(cls, start_index, stop_index=None, **kwargs):
        if isinstance(start_index, cls.mro()[1]):  # instance of parent
            return cls(start_index)
        if isinstance(start_index, _Sequence) and len(start_index) == 2:
            start_index, stop_index = start_index

        if stop_index is None:
            return cls(start_index + 1, **kwargs)
        return cls(start_index + 1, stop_index + 1, **kwargs)

    @classmethod
    def from_center_and_window(cls, center, window_size, max_length=_math.inf, **kwargs):
        start = max(1, center - window_size)
        stop = min(max_length, center + window_size)
        return cls(start, stop, **kwargs)

    @classmethod
    def _stop_from_start_and_length(cls, start, length):
        return SequencePoint(start + length - 1)

    @classmethod
    def from_slice(cls, start_slice, stop_slice=None):
        if isinstance(start_slice, slice):
            return cls(start_slice.start + 1, start_slice.stop)
        return cls(start_slice + 1, stop_slice)

    @classmethod
    def from_sequence(cls, full_sequence, peptide_sequence):
        start_index = full_sequence.find(peptide_sequence)
        if start_index == -1:
            raise IndexError("{} not in {}".format(peptide_sequence, full_sequence))

        start = start_index + 1
        stop = start_index + len(peptide_sequence)
        return cls(start, stop, seq=peptide_sequence, full_sequence=full_sequence)

    def _get_seq(self, seq, full_sequence):
        if seq:
            return str(seq)
        elif full_sequence:
            return full_sequence[self.slice]
        return None

    # implementation of abstract methods
    def validate(self):
        if self.start.pos < 1:
            raise ValueError("start < 1")
        if self.stop.pos < self.start.pos:
            raise ValueError("stop({}) < start({})".format(self.stop.pos, self.start.pos))

    def _comparison_cast(self, other):
        if super()._comparison_cast(other):
            return True
        elif isinstance(other, _Sequence):
            return len(other) == 2 and isinstance(other[0], int) and isinstance(other[1], int)
        return False

    def _join(self, other, operator):
        start = operator(self.start, other.start)
        stop = operator(self.stop, other.stop)
        return self.__class__(start, stop, validate=False)

    # dunders
    def __len__(self):
        return self.length

    def __str__(self):
        if self.start == self.stop:  # call SequencePoints.__str__
            return str(self.start)
        return "{}{}{}".format(self.start.pos, self._str_separator, self.stop.pos)

    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, self.start.pos, self.stop.pos)

    def __iter__(self):
        for pos in range(self.start.pos, self.stop.pos + 1):
            yield SequencePoint(pos, validate=False)

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

    def _contains(self, sequence_point):
        "Helper method that checks if a SequencePoint is in self"
        return self.start.pos <= sequence_point.pos <= self.stop.pos

    # properties, to make it read-only
    @property
    def seq(self):
        return self._seq

    @property
    def start(self):
        return SequencePoint(self._start)

    @property
    def stop(self):
        return SequencePoint(self._stop)

    @property
    def index(self):
        return _Index(self.start.index, self.stop.index)

    @property
    def pos(self):
        return _Pos(self.start.pos, self.stop.pos)
