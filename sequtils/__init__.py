"""
#. :code:`SequencePoint`, useful for emulating Mutations, SNPs, PTM's etc., it's two
   most important attributes are:

   - :code:`SequencePoint.pos`, the human readable number, counting from 1
   - :code:`SequencePoint.index`, the python readable number counting from 0
	
#. :code:`SequenedRange`, useful for emulating Proteins, domains, secondary structure etc.

   - Its 3 most important attributes are:

       - :code:`SequenceRange.start` is a :code:`SequencePoint` pointing to the first amino acid
       - :code:`SequenceRange.stop` is a :code:`SequencePoint` pointing to the last amino acid
       - :code:`SequenceRange.slice[start, stop]`: The python slice object, to index strings

   - It also has the following two properties for easy conversion to tuple

       - :code:`SequencePoint.pos.[start, stop]`: tuple with
         (:code:`self.start.pos`, :code:`self.stop.pos`)
       - :code:`SequencePoint.index[start, stop]`: tuple with
         (:code:`self.start.index`, :code:`self.stop.index`)

   - SequencePoint.slice[start, stop]: The python slice object, to index strings

For Developers:

:code:`SequencePoint` and :code:`SequenceRange` both subclass The base class
:code:`BaseSequenceLocation`, which defines most the dunders and :code:`_arithmetic` and
:code:`_comparison_cast` which takes care of most of the math and comparason for the subclasses
"""


# __future__ imports
from __future__ import annotations

# core imports
import collections as _collections
from collections.abc import Sequence as _Sequence
#  from abc import ABCMeta as _ABCMeta
from abc import abstractmethod as _abstractmethod
import operator as _operator
from functools import total_ordering as _total_ordering
import math as _math
from warnings import warn as _warn
import pathlib as _pathlib
from typing import Union as _Union

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
    #  __metaclass__ = _ABCMeta

    # read only attributes
    @property
    def pos(self):
        raise NotImplementedError("Please Implement this method")
    
    @property
    def index(self):
        raise NotImplementedError("Please Implement this method")

    @property
    def slice(self):
        raise NotImplementedError("Please Implement this method")

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
                SequencePoint.__init__ who raises an TypeError, thus
            SequencePoint_arthmetic returns NotImplemented,
                This signals to the Python Interperter to use
            SequenceRange.__radd__ 
                which works, because you can convert Points to Ranges but not vice versa!
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


# Type definitons for typehints, have to be defined after BaseSequenceLocation
# SequencePoint has to be a string because otherwise: SequencePoint was refrence before defined
_range_types = _Union[str, int, float, _Sequence, BaseSequenceLocation]
_point_types = _Union[str, int, float, 'SequencePoint']


class SequencePoint(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent

    :param position: Human Readable Sequence Position counting from 1

    .. code-block:: python

        # protein:      ELVISLIVES
        #  - positions: 1234567890
        # peptide:           LIVE
        #  - positions       6  9

        >>> seq = "ELVISLIVES"
        >>> mutation = SequencePoint(6)
        >>> mutation.pos 
        6

        >>> seq[mutation.index]
        'L'
    """

    def __new__(cls, arg, *args, **kwargs):
        if isinstance(arg, cls):
            # if arg is already of the correct type, then keep it because subclasses are mutable
            # just like other immutables: x = 213124512421312; x is int(x)
            return arg
        return super().__new__(cls)

    def __init__(self, position: _point_types, *, validate=True):
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
        """
        Alternative Constructure, using python indexes

        :param index: python index of position
        """
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

    # needed for pickle protocol 3+
    def __getnewargs__(self):
        return (self.pos,)

    @property
    def pos(self):
        return self._pos

    @property
    def index(self):
        return self._index

    @property
    def slice(self):
        return self._slice

class SequenceRange(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent
    
    :param start:
        | Human readable beginning of sequence
        | Can be of type:
        | * any type that can be convereted to int, 
        | * any :code:`collections.abc.Sequence` of length 2 
        | * any :code:`SequenceRange` or :code:`SequencePoint`
    :param stop:
        | Human readable end of sequence, can be infered from :code:`seq` and :code:`length`
        | Can be of type:
        * any type that can be convereted to int, or :code:`SequencePoint`
        * None it can be infered from :code:`start`, :code:`seq` or :code:`length`.
    :param seq:
        The biological sequence usally peptide or motif (eg :code:`CAT` or :code:`ELVISLIVES`)
    :param full_sequence:
        The biological sequence index by start and stop, :code:`seq` will be sliced out of
        :code:`full_sequence`
    :param length: length of the sequence
    :param validate:
        | raise exception if stop < start, is often set to false internally when performing math
            because:
        .. code-block:: python

            SequenceRange(5, 10) - (3, 2)
        | temporarely creates the :code:`SequenceRange(3, 2)` object
    :type start:
    :type stop: str, int or SequencePoint
    :type seq: str
    :type full_sequence: str
    :type validate: bool
    :type length: int

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

    Example

    .. code-block:: python

        # protein:     ELVISLIVES
        # - positions: 1234567890
        # peptide:          -----
        # - positions       6  9

        seq = "ELVISLIVES"
        peptide = SequenceRange(6, 9, seq)
        peptide.pos  # returns (6, 9)
        peptide.slice  # slice(5, 9, None)
        seq[peptide.slice]  # "LIVE"
        seq[peptide.index.start]  # "L"
        seq[peptide.index.stop]  # "E"
        peptide2 = SequenceRange(1, 5)
        peptide2 < peptide1  # True, because (1, 5) < (5, 9)

    """

    _str_separator = ':'

    def __init__(self, start: _range_types, stop: _point_types=None, seq: _Union[None, str]=None,
                 full_sequence: _Union[None, str]=None, *, validate: bool=True,
                 length: _Union[int, bool]=None):
        if isinstance(start, self.__class__.mro()[1]):  # isinstance of parent
            if isinstance(start, self.__class__):
                if seq is None and full_sequence is None:
                    seq = start.seq
                start, stop = start.pos
            elif isinstance(start, SequencePoint):
                start = start
        if isinstance(start, _Sequence) and len(start) == 2 and \
                not isinstance(start, (bytes, str)):
            if stop is not None:
                raise ValueError(
                    "if start is a Sequence (e.g (2, 5), then stop has to be None\n"
                    " - Thise is ok : SequenceRange(2, 5)\n"
                    " -             : SequenceRange((2, 5))\n"
                    " - but not this: SequenceRange((2, 5), 5)\n")
            start, stop = start[:2]
        elif isinstance(start, str) and self._str_separator in start and stop is None:
            start, stop = map(int, start.split(':'))

        stop = self._resolve_stop(start, stop, length, seq)
        self._seq = seq

        self._start = SequencePoint(start, validate=validate)
        self._stop = SequencePoint(stop, validate=validate)

        if validate:
            self.validate()

        self._slice = slice(self.start.pos - 1, self.stop.pos)
        self._seq = self._get_seq(seq, full_sequence)

    def _resolve_stop(self, start, stop, length, seq):
        if isinstance(stop, SequencePoint):
            stop = stop
        elif stop is None:
            if length is not None:
                stop = self._stop_from_start_and_length(start, length)
            elif seq is not None:
                stop = self._stop_from_start_and_length(start, len(seq))
            else:
                stop = start
        return stop

    # alternate constructors
    @classmethod
    def from_index(cls, start_index: _range_types, stop_index: _point_types=None, **kwargs):
        """
        Alternative Constructure, using python indexes

        :param start_index: python index of start position
        :param stop_index: python index of stop position
        """
        if isinstance(start_index, cls.mro()[1]):  # instance of parent
            return cls(start_index)
        if isinstance(start_index, _Sequence) and len(start_index) == 2:
            start_index, stop_index = start_index

        if stop_index is None:
            return cls(start_index + 1, **kwargs)
        return cls(start_index + 1, stop_index + 1, **kwargs)

    #  @classmethod
    #  def _from_string(cls, string):
    #      if cls._str_separator in string:
    #          start, stop = map(SequencePoint, string.split(cls._str_separator))
    #      else:
    #          start = stop = SequencePoint(string, validate=False)
    #      return start, stop

    @classmethod
    def from_center_and_window(cls, center: _Union[int, SequencePoint], window: int,
                               max_length: _point_types=_math.inf, **kwargs):
        """
        Alternative Constructur, :math:`center\pm{}window`

        :param center: central position
        :param window: extension to the left and right of :code:`center`
        :param max_length: `stop` cannot be extended pass this number

        example:

        .. code-block:: python

            >>> SequenceRange.from_center_and_window(10, 5)
            SequenceRange(5, 15, seq=None)

            # cannot extend past pos=1
            >>> SequenceRange.from_center_and_window(10, 15)
            SequenceRange(1, 25, seq=None)

            # cannot extend past max_length
            >>> SequenceRange.from_center_and_window(100, 10, max_length=105)
            SequenceRange(90, 105, seq=None)
        """

        start = max(1, center - window)
        stop = min(max_length, center + window)
        return cls(start, stop, **kwargs)

    @classmethod
    def _stop_from_start_and_length(cls, start, length):
        return SequencePoint(start + length - 1)

    @classmethod
    def from_slice(cls, start_slice: _Union[int, _Sequence, slice],
                   stop_slice: _Union[int, None]=None):
        """
        Alternative Constructor, from python slice, or slice coordinates

        :param start_sclice

        """

        if isinstance(start_slice, slice):
            return cls(start_slice.start + 1, start_slice.stop)
        elif isinstance(start_slice, _Sequence) and len(stop_slice) == 2:
            start_slice, stop_slice = start_slice
        return cls(start_slice + 1, stop_slice)

    @classmethod
    def from_sequence(cls, full_sequence: str, sequence: str):
        """
        Alternative Constructure, cut out the sequence from a full_sequence (protein, gene, ect)

        :param full_sequence: a biological sequence
        :param sequence: a biological sequence contained within :code:`full_sequence`

        Example:

        .. code-block:: python

        >>> SequenceRange.from_sequence('EVILELVISLIVES', 'ELVIS')
        SequenceRange(5, 9, seq="ELVIS")

        **Warning:** if code:`sequence` is found multiple times in :code:`full_sequence`, then the
        first occurance will be returned
        """

        start_index = full_sequence.find(sequence)
        if start_index == -1:
            raise IndexError("{} not in {}".format(sequence, full_sequence))

        start = start_index + 1
        stop = start_index + len(sequence)
        return cls(start, stop, seq=sequence, full_sequence=full_sequence)

    def _get_seq(self, seq, full_sequence):
        if seq:
            seq = str(seq)
            if len(seq) != len(self):
                msg = "The sequence {} length does not match the one implied by {}"
                raise ValueError(msg.format(full_sequence, self))
        elif full_sequence:
            seq = full_sequence[self.slice]
            if len(seq) != len(self):
                msg = "The sequence {} is to short to contain {}"
                raise ValueError(msg.format(full_sequence, self))
        return seq
        #  return None

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
        if other.start == other.stop:
            # if length is unchanged after math, then keep seq
            return self.__class__(start, stop, seq=self.seq, validate=False)
        return self.__class__(start, stop, validate=False)

    # dunders
    def __len__(self):
        return self.length

    def __str__(self):
        if self.start == self.stop:  # call SequencePoints.__str__
            return str(self.start)
        return "{}{}{}".format(self.start.pos, self._str_separator, self.stop.pos)

    def __repr__(self):
        base = "{}({{}})".format(type(self).__name__)
        if self.seq is None:
            return base.format("{}, {}, seq=None".format(self.start, self.stop))

        seq = self.seq
        if 40 < len(seq):
            seq = self.seq[:5] + '..' + self.seq[-5:]
        return base.format('{}, {}, seq="{}"'.format(self.start, self.stop, seq))

        #  if self.seq is None:
        #      return '{}({}, {}, seq=None)'.format(type(self).__name__, self.start, self.stop)
        #  return '{}({}, {}, seq="{}")'.format(type(self).__name__, self.start, self.stop, self.seq)

    def __iter__(self):
        for pos in range(self.start.pos, self.stop.pos + 1):
            yield SequencePoint(pos, validate=False)

    def __contains__(self, item):
        " returns True if all of item is inside self"
        return self.contains(item, part=all)

    def contains(self, item, part=all):
        """
        Check wheter item is inside 'self'.

        | by default (part=all) all amino acids has to be inside self:
        | by changing to (part=any), then only one of the amino acids have to be inside:

        .. code-block::

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
    def seq(self) -> _Union[str, None]:
        return self._seq

    @property
    def start(self) -> SequencePoint:
        return self._start

    @property
    def stop(self) -> SequencePoint:
        return self._stop

    @property
    def slice(self) -> _Index:
        return self._slice

    @property
    def index(self) -> _Index:
        return _Index(self.start.index, self.stop.index)

    @property
    def pos(self) -> _Pos:
        return _Pos(self.start.pos, self.stop.pos)

    @property
    def length(self) -> int:
        return self.stop.pos - self.start.pos + 1

    def __eq__(self, other):
        return self._eq_helper(other)

    def _eq_helper(self, other, *, compare_seq=True):
        if self._comparison_cast(other):
            try:
                other = self.__class__(other, validate=False)
                #  return self.pos == other.pos and self.seq == other.seq
                return self.pos == other.pos and (not compare_seq or self.seq == other.seq)
            except (ValueError, TypeError):
                pass
        return NotImplemented

    @classmethod
    def _hash(cls, pos, seq):
        return hash((pos, seq))

    def __hash__(self):
        return self._hash(self.pos, self.seq)

    def equals(self, other, compare_seq=True, cast=True):
        """
        With default arguments it works much like :code:`==`, but the behavior can be altered

        .. code-block::

            >>> sr_seq = SequenceRange(10, seq="A"*11)
            >>> sr_non = SequenceRange(10, 20)

            >>> # default == behaviour
            >>> sr_seq == (10, 20)  # returns False - because seq is different
            False
            >>> sr_non == (10, 20)  # returns True - because both seq are None
            True

            >>> sr_non.equals((10, 20), compare_seq=True, cast=False)  # no cast
            False
            >>> sr_non.equals(sr_seq, compare_seq=True, cast=False)  # different seq
            False
            >>> sr_non.equals(sr_seq, compare_seq=False, cast=False)  # same coordinates
            True
        """

        if cast == True:
            return self._eq_helper(other, compare_seq=compare_seq)
        elif isinstance(other, self.__class__.mro()[1]):  # isinstance of parent
            return self._eq_helper(other, compare_seq=compare_seq)
        return False





