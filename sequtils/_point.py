# core imports
from typing import Union
from collections.abc import Sequence

# local imports
from ._base import BaseSequenceLocation


point_types = Union[str, int, float, 'SequencePoint']


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

    def __new__(cls, position, *args, **kwargs):
        if isinstance(position, cls.mro()[1]):  # isinstance of parent
            if isinstance(position, cls):
                # if arg is already of the correct type, then keep it because subclasses are
                # mutable just like other immutables: x = 213124512421312; x is int(x)
                return position
            elif hasattr(position, 'start') and hasattr(position, 'stop'):
                if len(position) == 1:
                    return position.start
                raise TypeError("can only Convert {} to {} if len({}) = 1".format(
                    type(position), cls, type(position)))
        return super().__new__(cls)

    def __init__(self, position: point_types, *, validate=True):
        if isinstance(position, self.__class__.mro()[1]):  # isinstance of parent
            return

        self._pos = int(position)
        if validate:
            self.validate()
        self._index = self._pos - 1
        self._slice = slice(self.index, self.index + 1)

    # alternative constructors
    @classmethod
    def from_index(cls, index: Union[int, Sequence], *, validate=True):
        """
        Alternative Constructure, using python indexes

        :param index: python index of position
        """

        if isinstance(index, cls.mro()[1]):  # isinstance of parent
            raise ValueError("index, cannot be of type {}".format(repr(index)))
        return cls(index + 1, validate=validate)

    # implementation of abstract methods
    def validate(self):
        if self.pos < 1:
            raise ValueError("position({}) < 1".format(self.pos))

    def _join(self, other, operator):
        return self.__class__.from_index(operator(self.index, other.index), validate=False)

    # implementation of abstract properties
    @property
    def pos(self):
        return self._pos

    @property
    def index(self):
        return self._index

    @property
    def slice(self):
        return self._slice

    # dunders
    def __str__(self):
        return str(self.pos)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.pos)

    # needed for pickle protocol 3+
    def __getnewargs__(self):
        return (self.pos,)


