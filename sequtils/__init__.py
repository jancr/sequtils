"""
#. :code:`SequencePoint`, useful for emulating Mutations, SNPs, PTM's etc., it's two
   most important attributes are:

   - :code:`SequencePoint.pos`, the human readable number, counting from 1
   - :code:`SequencePoint.index`, the python readable number counting from 0

#. :code:`SequenedRange`, useful for emulating Proteins, domains, secondary structure etc.

   - Its 3 most important attributes are:

       - :code:`SequenceRange.start` is a :code:`SequencePoint` pointing to the first amino acid
       - :code:`SequenceRange.stop` is a :code:`SequencePoint` pointing to the last amino acid
       - :code:`SequenceRange.slice.[start, stop]`: The python slice object, to index strings

   - It also has the following two properties for easy conversion to tuple

       - :code:`SequencePoint.pos.[start, stop]`: tuple containing
         (:code:`self.start.pos`, :code:`self.stop.pos`)
       - :code:`SequencePoint.index[start, stop]`: tuple containing
         (:code:`self.start.index`, :code:`self.stop.index`)

For Developers:

:code:`SequencePoint` and :code:`SequenceRange` both subclass The base class
:code:`BaseSequenceLocation`, which defines most the dunders and :code:`_arithmetic` and
:code:`_comparison_cast` which takes care of most of the math and comparason for the subclasses
"""


from ._point import SequencePoint
from ._range import SequenceRange


__slots__ = ("SequencePoint", "SequenceRange")
__all__ = ("SequencePoint", "SequenceRange")
