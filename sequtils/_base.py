# core imports
import abc
import functools
import operator


@functools.total_ordering
class BaseSequenceLocation(metaclass=abc.ABCMeta):
    #  __metaclass__ = _ABCMeta

    # read only attributes
    @property
    @abc.abstractmethod
    def pos(self):
        pass
        #  raise NotImplementedError("Please Implement this method")

    @property
    @abc.abstractmethod
    def index(self):
        pass
        #  raise NotImplementedError("Please Implement this method")

    @property
    @abc.abstractmethod
    def slice(self):
        pass
        #  raise NotImplementedError("Please Implement this method")

    # abstract methods
    @abc.abstractmethod
    def _join(self, other, operator):
        "helper methood, needed to make __add__ and __sub__ work"
        pass
        #  raise NotImplementedError("Please Implement this method")

    @abc.abstractmethod
    def validate(self):
        pass
        #  raise NotImplementedError("Please Implement this method")

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


        # if it is a subclass of this class, then try to cast
        cls = type(self)
        if isinstance(other, BaseSequenceLocation):
            if not isinstance(other, cls):
                try:
                    other = cls(other, validate=False)
                except TypeError:
                    return NotImplemented
        # if not, then try to create one assuming that 'other' is an index
        else:
            try:
                other = cls.from_index(other, validate=False)
            except TypeError:
                return NotImplemented
        return self._join(other, operator)

        #  if not isinstance(other, self.__class__):
        #      try:
        #          #  other = getattr(other, 'index', other)
        #          other = self.__class__.from_index(other, validate=False)
        #      except TypeError:
        #          return NotImplemented
        #  return self._join(other, operator)

    #  math dunders
    def __add__(self, other):
        return self._arithmetic(other, operator.add)

    def __sub__(self, other):
        return self._arithmetic(other, operator.sub)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        # other - self
        return self.__class__.from_index(other, validate=False) - self

    # comparison dunders
    def __lt__(self, other):
        if self._comparison_cast(other):
            #  try:
            return self.pos < self.__class__(other, validate=False).pos
            #  except ValueError:
            #      raise TypeError("cannot compare type {} and {}".format(type(self), type(other)))
        return NotImplemented

    def __eq__(self, other):
        if self._comparison_cast(other):
            #  try:
            return self.pos == self.__class__(other, validate=False).pos
            #  except (ValueError, TypeError):
            #      pass
        return NotImplemented

    # other dunders
    def __hash__(self):
        return hash(self.pos)
