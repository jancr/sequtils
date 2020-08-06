========
sequtils
========

Collection of Classes and functions for working with biological sequences

There are two classes
=====================

1. `SequencePoint`, useful for emulating Mutations, SNPs, PTM's etc., it's two
   most important attributes are:
   - `SequencePoint.pos`, the human readable number, counting from 1
   - `SequencePoint.index`, the python readable number counting from 0
	
2. `SequenedRange`, useful for emulating Proteins, domains, secondary structure etc.
   - Its 3 most important attributes are:
       - `SequenceRange.start` is a `SequencePoint` pointing to the first amino acid
       - `SequenceRange.stop` is a `SequencePoint` pointing to the last amino acid
       - `SequenceRange.slice[start, stop]`: The python slice object, to index strings
   - It also has the following two properties for easy conversion to tuple
       - `SequencePoint.pos.[start, stop]`: tuple with (`self.start.pos`, `self.stop.pos`)
       - `SequencePoint.index[start, stop]`: tuple with (`self.start.index`, `self.stop.index`)

   - SequencePoint.slice[start, stop]: The python slice object, to index strings

Example code, lets make glucagon

```
In [1]: from sequtils import SequenceRange, SequencePoint

In [2]: glucagon_sequence = ("MKTIYFVAGLLIMLVQGSWQHALQDTEENPRSFPASQTEAHEDPDEMNEDKRHSQGTFTS"
   ...:                      "DYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFI"
   ...:                      "AWLVKGRGRRDFPEEVAIAEELGRRHADGSFSDEMSTILDNLATRDFINWLIQTKITDKK")

In [3]: glucagon = SequenceRange(1, seq=glucagon_sequence)

In [3]: glucagon
Out[3]: SequenceRange(1, 180, seq="KM..KK")
```

So we now have a protein object, where the stop was inferred from the sequence, `glp1` is a peptide

```
In [4]: glp1 = SequenceRange(98, 127, full_sequence=glucagon_sequence)
In [5]: glp1
Out[5]: SequenceRange(98, 127, seq="HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR")
```

A `SequenceRange` from 98 to 127 is created, with the peptide sequence inferred
from the protein sequence

Lets see the `start` and `stop` attributes of the peptide:

```
In [6]: glp1.start
Out[6]: SequencePoint(98)

In [7]: glucagon_sequence[glp1.start.index] == glp1.seq[0]
Out[7]: True

In [7]: glp1.stop
Out[7]: SequencePoint(127)

In [8]: glucagon_sequence[glp1.stop.index] == glp1.seq[-1]
Out[8]: True

```

Lets try to use the slice object to cut the peptide sequence out of the protein

```
In [9]: glp1.slice
Out[9]: slice(97, 127, None)

In [10]: glucagon_sequence[glp1.slice]
Out[10]: 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'

In [11]: glp1.seq == glucagon.seq[glp1.slice]
Out[11]: True
```

GLP-1 is famous for having a canonical G\[KR\]\[KR\] motif, this motif is the 3
N-terminal flaking amino acids, let's find it

```
In [12]: motif = SequenceRange(1 + glp1.pos.stop, 3 + glp1.pos.stop)

In [13]: glucagon.seq[motif.slice]
Out[13]: 'GRR'
```

## Math API examples

The objects also supports math... So lets try to do the above with math, but first an explanation.

All math on these objects are performed based on the Indexes, thus

```
In [14]: SequencePoint(1) + SequencePoint(1)
Out[14]: SequencePoint(1)

In [15]: SequenceRange(1, 1) + SequenceRange(1, 1)
Out[15]: SequenceRange(1, 1, seq=None)
```

Because `SequencePoint(1).index` is 0 and 0 + 0 = 0

The above code is equivalent to the following:

```
In [16]: SequencePoint.from_index((SequencePoint(1).index + SequencePoint(1).index))
Out[16]: SequencePoint(1)
```

The math is super intuitive for scalars

```
In [17]: SequenceRange(2, 5) + 2
Out[17]: SequenceRange(4, 7, seq=None)

In [18]: SequenceRange(2, 5, seq="EVIL") + 2
Out[18]: SequenceRange(4, 7, seq=EVIL)
```

It also works for non scalars, but then seq becomes `None` because the length has changed

```
In [19]: SequenceRange(2, 5, seq="EVIL") + SequenceRange(3, 6)
Out[19]: SequenceRange(4, 10, seq=None)
```

If you add numbers or tuples, the code will assume that those are indexes,
thus the following 3 all gives the GRR motif by moving `glp1.stop` by `(1, 3)`

Create new object moving `glp1.stop`

```
In [20]: SequenceRange(1 + glp1.stop, 3 + glp1.stop)
Out[20]: SequenceRange(128, 130, seq=None)
```

Create new object via math, here we perform `SequenceRange` + `SequencePoint`

```
In [20]: glp1.stop + SequenceRange.from_index(1, 3)
Out[20]: SequenceRange(128, 130, seq=None)

In [21]: glp1.stop + SequenceRange(2, 4)
Out[21]: SequenceRange(128, 130, seq=None)
```

Convert `SequencePoint` to `SequenceRange` and then add an offset tuple, **note**
that `SequencePoint` only knows 'scalar' math, so we have to ether convert it
to a `SequenceRange` as here, or convert the `(1, 3)` tuple to a `SequnceRange`
as we did above

```
In [22]: SequenceRange(glp1.stop) + (1, 3)
Out[22]: SequenceRange(128, 130, seq=None)
```
