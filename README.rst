=================
Sequtils Tutorial
=================

Collection of Classes and functions for working with biological sequences

Overview
========

There are two Public classes

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

       - :code:`SequencePoint.pos.[start, stop]`: tuple containing
         (:code:`self.start.pos`, :code:`self.stop.pos`)
       - :code:`SequencePoint.index.[start, stop]`: tuple containing
         (:code:`self.start.index`, :code:`self.stop.index`)

Example Usage
=============

Example code, lets make glucagon

.. code-block:: python

    >>> from sequtils import SequenceRange, SequencePoint
    >>> glucagon_sequence = ("MKTIYFVAGLLIMLVQGSWQHALQDTEENPRSFPASQTEAHEDPDEMNEDKRHSQGTFTS"
    ...                      "DYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFI"
    ...                      "AWLVKGRGRRDFPEEVAIAEELGRRHADGSFSDEMSTILDNLATRDFINWLIQTKITDKK")
    >>> glucagon = SequenceRange(1, seq=glucagon_sequence)
    >>> glucagon
    SequenceRange(1, 180, seq="MKTIY..ITDKK")

So we now have a protein object, where the stop was inferred from the sequence, :code:`glp1` is a peptide

.. code-block:: python

    >>> glp1 = SequenceRange(98, 127, full_sequence=glucagon_sequence)
    >>> glp1
    SequenceRange(98, 127, seq="HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR")

A :code:`SequenceRange` from 98 to 127 is created, with the peptide sequence inferred
from the protein sequence

Lets see the :code:`start` and :code:`stop` attributes of the peptide:

.. code-block:: python

    >>> glp1.start
    SequencePoint(98)

    >>> glucagon_sequence[glp1.start.index] == glp1.seq[0]
    True

    >>> glp1.stop
    SequencePoint(127)

    >>> glucagon_sequence[glp1.stop.index] == glp1.seq[-1]
    True

Lets try to use the slice object to cut the peptide sequence out of the protein

.. code-block:: python

    >>> glp1.slice
    slice(97, 127, None)

    >>> glucagon_sequence[glp1.slice]
    'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'

    >>> glp1.seq == glucagon.seq[glp1.slice]
    True

GLP-1 is famous for having a canonical G\[KR\]\[KR\] motif, this motif is the 3
N-terminal flaking amino acids, let's find it

.. code-block:: python

    >>> motif = SequenceRange(1 + glp1.stop.pos, 3 + glp1.stop.pos)
    >>> glucagon.seq[motif.slice]
    'GRR'

Math Examples
=============

The objects also supports math... So lets try to do the above with math, but first an explanation.

All math on these objects are performed based on the Indexes, thus

.. code-block:: python

    >>> SequencePoint(1) + SequencePoint(1)
    SequencePoint(1)

    >>> SequenceRange(1, 1) + SequenceRange(1, 1)
    SequenceRange(1, 1, seq=None)

Because :code:`SequencePoint(1).index` is :math:`0` and :math:`0 + 0 = 0`

The above code is equivalent to the following:

.. code-block:: python

    >>> SequencePoint.from_index((SequencePoint(1).index + SequencePoint(1).index))
    SequencePoint(1)

The math is super intuitive for scalars

.. code-block:: python

    >>> SequenceRange(2, 5) + 2
    SequenceRange(4, 7, seq=None)

    >>> SequenceRange(2, 5, seq="EVIL") + 2
    SequenceRange(4, 7, seq="EVIL")

It also works for non scalars, but then seq becomes :code:`None` because the length has changed

.. code-block:: python

    >>> SequenceRange(2, 5, seq="EVIL") + SequenceRange(3, 6)
    SequenceRange(4, 10, seq=None)

If you add numbers or tuples, the code will assume that those are indexes,
thus the following 3 all gives the GRR motif by moving :code:`glp1.stop` by :code:`(1, 3)`

Create new object moving :code:`glp1.stop`

.. code-block:: python

    >>> SequenceRange(glp1.stop + 1, glp1.stop + 3)
    SequenceRange(128, 130, seq=None)

Create new object via math, here we perform :code:`SequenceRange` + :code:`SequencePoint`

.. code-block:: python

    >>> glp1.stop + SequenceRange.from_index(1, 3)
    SequenceRange(128, 130, seq=None)

    >>> glp1.stop + SequenceRange(2, 4)
    SequenceRange(128, 130, seq=None)

Convert :code:`SequencePoint` to :code:`SequenceRange` and then add an offset tuple, **note**
that :code:`SequencePoint` only knows 'scalar' math, so we have to ether convert it
to a :code:`SequenceRange` as here, or convert the :code:`(1, 3)` tuple to a :code:`SequnceRange`
as we did above

.. code-block:: python

    >>> SequenceRange(glp1.stop) + (1, 3)
    SequenceRange(128, 130, seq=None)
    
