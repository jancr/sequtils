
# core imports
import collections
from collections.abc import Sequence
import math
import warnings
from typing import Union

# local imports
from ._base import BaseSequenceLocation
from ._point import SequencePoint, point_types


range_types = Union[str, int, float, Sequence, BaseSequenceLocation]


class _Positions(collections.namedtuple("Pos", ("_1", "_2"), rename=True)):
    """
    Wrapped namedtuple, to issue deprecation warnings when the interface is used the old way

    The old SequenceRange object had two named tuples one named index and one named pos, so it was
    accessed like this

    .. code:: python

        >>> sr = SequenceRange(5, 10)
        >>> sr.index  # <-- this is fine
        (4, 9)
        >>> sr.index.start  # <-- depricated
        4

    This behaviour is still valid, but throws a :code:`DeprecationWarning` via this class, the
    recomended syntax is to first do :code:`start` or :code:`stop` which returns a
    :code:`SequencePoint`

    .. code:: python

        >>> sr = SequenceRange(5, 10)
        >>> sr.start
        SequencePoint(5)
        >>> sr.start.index  # <-- new syntax
        4
    """

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
        warnings.warn(msg, DeprecationWarning)
        return field

    def __repr__(self):
        return repr(tuple(self))


class _Pos(_Positions):
    name = "pos"


class _Index(_Positions):
    name = "index"


# Type definitons for typehints, have to be defined after BaseSequenceLocation
# SequencePoint has to be a string because otherwise: SequencePoint was refrence before defined
class SequenceRange(BaseSequenceLocation):
    """
    helper class that converts between "normal" sequence numbers and pythons equivalent

    sequence cordinate, counting from 1, with inclusive stop
    viable calls to the constructor includes:
        | * :code:`SequenceRange(1, 2)`
        | * :code:`SequenceRange(SequencePoint(1), SequencePoint(2))`
        | * :code:`SequenceRange((1, 2))`
        | * :code:`SequenceRange(SequenceRange(1, 2))`

    A "start" has to be always provided, but stop can be infered from the other arguments
    thus the following creates the same object:
        | * :code:`SequenceRange(1, 3)`
        | * :code:`SequenceRange((1, 3))`
        | * :code:`SequenceRange(1, seq="ABC")`
        | * :code:`SequenceRange(1, length=3)`

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
            because the two lines are equivalent:
        .. code-block:: python

            SequenceRange(5, 10) - (3, 2)
            SequenceRange(5, 10) - SequenceRange.from_index(3, 2, validate=False)
    :type start:
    :type stop: str, int or SequencePoint
    :type seq: str
    :type full_sequence: str
    :type validate: bool
    :type length: int

    Example

    .. code-block:: python

        # protein:     ELVISLIVES
        # - positions: 1234567890
        # peptide:          -----
        # - positions       6  9

        >>> seq = "ELVISLIVES"
        >>> peptide = SequenceRange(6, 9, full_sequence=seq)
        >>> peptide.pos  
        (6, 9)
        >>> peptide.slice
        slice(5, 9, None)
        >>> seq[peptide.slice]
        'LIVE'

        >>> seq[peptide.start.index]
        'L'
        >>> seq[peptide.stop.index]  
        'E'

        >>> peptide2 = SequenceRange(1, 5)
        >>> peptide2 < peptide  # True, because (1, 5) < (5, 9)
        True

    """

    _str_separator = ':'
    #  _bytes_seperator = b':'  # this should be a class decorator created from _str_seperator

    def __init__(self, start: range_types, stop: point_types=None, seq: Union[None, str]=None,
                 full_sequence: Union[None, str]=None, *, validate: bool=True,
                 length: Union[int, bool]=None, _special=None):
        if isinstance(start, self.__class__.mro()[1]):  # isinstance of parent
            if isinstance(start, self.__class__):
                if stop is not None:
                    raise ValueError("either:\n"
                                     "    1. start and stop are scalars or\n"
                                     "    2. start is asequence of length 2 and stop is None")
                if seq is None and full_sequence is None:
                    seq = start.seq
                start, stop = start.pos
            #  elif isinstance(start, SequencePoint):
            #      start = start.pos
        elif self._valid_range(start):
            start, stop = self._parse_range(start, stop)

        if isinstance(start, (str, bytes)):
            start = int(start)
        if isinstance(stop, (str, bytes)):
            stop = int(stop)
        elif stop is None:
            stop = self._resolve_none_stop(start, stop, length, seq)

        self._seq = seq

        if _special == 'index':
            start_offset, stop_offset = (1, 1)
        elif _special == 'slice':
            start_offset, stop_offset = (1, 0)
        else:
            start_offset, stop_offset = (0, 0)

        self._start = SequencePoint(start + start_offset, validate=validate)
        self._stop = SequencePoint(stop + stop_offset, validate=validate)

        if validate:
            self.validate()

        self._slice = slice(self.start.pos - 1, self.stop.pos)
        self._seq = self._get_seq(seq, full_sequence)

    def _resolve_none_stop(self, start, stop, length, seq):
        #  if isinstance(stop, (str, bytes)):
        #      return int(stop)
        if length is not None:
            return self._stop_from_start_and_length(start, length)
        elif seq is not None:
            return self._stop_from_start_and_length(start, len(seq))
        return start

    @classmethod
    def _valid_range(cls, positions):
        if not isinstance(positions, Sequence):
            return False
        if isinstance(positions, str):
            return cls._str_separator in positions
        elif isinstance(positions, bytes):
            return cls._str_separator.encode() in positions
        return len(positions) == 2

    @classmethod
    def _parse_range(cls, start, stop):
        if stop is not None:
            raise ValueError(
                "if start is a Sequence (e.g (2, 5), then stop has to be None\n"
                " - Thise is ok : SequenceRange(2, 5)\n"
                " -             : SequenceRange((2, 5))\n"
                " - but not this: SequenceRange((2, 5), 5)\n")
        if isinstance(start, bytes):
            start = start.decode('utf8')
        if isinstance(start, str) and cls._str_separator in start:
            return map(int, start.split(cls._str_separator))
        elif isinstance(start, Sequence) and len(start) == 2:
            return start
        raise ValueError("{} cannot be understood by the constructor".format(start))

    # alternate constructors
    @classmethod
    #  def from_index(cls, start_index: range_types, stop_index: point_types=None, **kwargs):
    def from_index(cls, start_index: Union[int, Sequence], stop_index: Union[int, None]=None,
                   **kwargs):
        """
        Alternative Constructure, using python indexes

        :param start_index: python index of start position
        :param stop_index: python index of stop position
        """

        if isinstance(start_index, cls.mro()[1]):  # isinstance of parent
            raise ValueError("start_index, cannot be of type {}".format(repr(start_index)))
        if isinstance(stop_index, cls.mro()[1]):  # isinstance of parent
            raise ValueError("stop_index, cannot be of type {}".format(repr(stop_index)))
        return cls(start_index, stop_index, _special='index', **kwargs)

    @classmethod
    def from_center_and_window(cls, center: Union[int, SequencePoint], window: int,
                               max_length: point_types=math.inf, **kwargs):
        r"""
        Alternative Constructur, :math:`center\pm{}window`

        :param center: central position
        :param window: extension to the left and right of :code:`center`
        :param max_length: `stop` cannot be extended pass this number

        Example:

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
    def from_slice(cls, start_slice: Union[int, Sequence, slice],
                   stop_slice: Union[int, None]=None, **kwargs):
        """
        Alternative Constructor, from python slice, or slice coordinates

        :param start_slice: slice object or where to slice from
        "param stop_slice: where to slice to or :code:`None` of :code:`start_slice` is a
        :code:`slice`

        .. code-block:: python

            >>> seq = "ELVISLIVES"
            >>> seq[5:9]
            'LIVE'
            >>> peptide = SequenceRange.from_slice(5, 9)
            >>> peptide.slice
            slice(5, 9, None)
            >>> seq[peptide.slice]
            'LIVE'
        """

        if isinstance(start_slice, slice):
            if start_slice.step not in (1, None):
                raise ValueError("Slice has to have a step of 1 to be a valid SequenceRange")
            return cls(start_slice.start + 1, start_slice.stop, **kwargs)
        return cls(start_slice, stop_slice, _special='slice', **kwargs)

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

    # implementation of abstract methods
    def validate(self):
        if self.start.pos < 1:
            raise ValueError("start < 1")
        if self.stop.pos < self.start.pos:
            raise ValueError("stop({}) < start({})".format(self.stop.pos, self.start.pos))

    def _comparison_cast(self, other):
        if super()._comparison_cast(other):
            return True
        elif isinstance(other, Sequence):
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
            except (ValueError, TypeError):
                pass
        return False

    def _contains(self, sequence_point):
        "Helper method that checks if a SequencePoint is in self"
        return self.start.pos <= sequence_point.pos <= self.stop.pos

    # properties, to make it read-only
    @property
    def seq(self) -> Union[str, None]:
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
            #  try:
            other = self.__class__(other, validate=False)
            return self.pos == other.pos and (not compare_seq or self.seq == other.seq)
            #  except (ValueError, TypeError):
            #      pass
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
