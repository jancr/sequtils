# sequtils
Collection of Classes and functions for working with biological sequences

## There are two classes
1. SequencePoint, usefull for emulation Mutations, SNPs, PTM's etc, it's two
   most important attributes are:
 - SequencePoint.pos, the human readable number, counting from 1
 - SequencePoint.index, the python readable number counting from 0
2. SequenedRange, usefull for emulating Proteins, domains, secondary structure etc.
 - SequencePoint.pos.[start, stop]: The human readable numbers, counting from 1 with inclusive end
 - SequencePoint.index[start, stop]: The python readable numbers counting from 0
 - SequencePoint.slice[start, stop]: The python slice object, to index strings
   
 - Its 3 most important attributes are:

example code, lets make glucagon
```
In [1]: from sequtils import SequenceRange, SequencePoint

In [2]: glucagon_sequence = ("MKTIYFVAGLLIMLVQGSWQHALQDTEENPRSFPASQTEAHEDPDEMNEDKRHSQGTFTS"
   ...:                      "DYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFI"
   ...:                      "AWLVKGRGRRDFPEEVAIAEELGRRHADGSFSDEMSTILDNLATRDFINWLIQTKITDKK")

In [3]: glucagon = SequenceRange(1, seq=glucagon_sequence)

In [3]: glucagon
Out[3]: SequenceRange(1, 180)
```

So we now have a protein object, where the stop was inferred from the sequence, glp1 is a peptide

```
In [7]: glp1 = SequenceRange(98, 127, full_sequence=glucagon_sequence)
In [8]: glp1
Out[8]: SequenceRange(98, 127)

In [9]: glp1.seq
Out[9]: 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'
```

A SequenceRange from 98 to 127 is created, with the peptide sequence inferred
from the protein sequence

Lets see the `pos`, `index` and `slice` attributes of the peptide:

```
In [10]: glp1.pos
Out[10]: Position(start=98, stop=127)

In [11]: glp1.index
Out[11]: Index(start=97, stop=126)

In [12]: glp1.slice
Out[12]: slice(97, 127, None)
```

Lets try to use the slice object to cut the peptide sequence out of the protein

```
In [13]: glucagon_sequence[glp1.slice]
Out[13]: 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'

In [14]: glp1.seq == glucagon.seq[glp1.slice]
Out[14]: True
```

GLP-1 is famos for having a canonical G\[KR\]\[KR\] motif, this motif is the 3
N-terminal flaking amino acids, let's find it

```
motif = S
In [15]: motif = SequenceRange(1 + glp1.pos.stop, 3 + glp1.pos.stop)

In [16]: glucagon.seq[motif.slice]
Out[16]: 'GRR'
```

TODO: explain math :D






