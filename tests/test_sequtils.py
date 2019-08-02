
# core imports
import os

# 3rd party imports
import pytest

# local imports
from sequtils import SequencePoint, SequenceRange

TEST_FOLDER = os.path.abspath(os.path.dirname(__file__))
TEST_FILES_FOLDER = os.path.abspath(os.path.join(TEST_FOLDER, 'test_files'))

########################################
# fixtures
# move to conftest when project becomes larger!
########################################
@pytest.fixture
def glucagon_seq():
    with open(os.path.join(TEST_FILES_FOLDER, 'glucagon.fasta')) as f:
        return "".join(line.strip() for line in f.readlines()[1:])


@pytest.fixture
def glucagon_peptides():
    """ all peptides from glucagon as ((start, stop, seq), ..) """

    def get_peptide(f):
        start, stop = map(int, f.readline().split('|')[-1].split('-'))
        seq = f.readline().strip()
        return start, stop, seq

    peptides = []
    with open(os.path.join(TEST_FILES_FOLDER, 'glucagon_peptides.fasta')) as f:
        while f.tell() != os.fstat(f.fileno()).st_size:
            peptides.append(get_peptide(f))
    return tuple(peptides)


def test_fixtures(glucagon_peptides, glucagon_seq):
    assert len(glucagon_peptides) == 11
    assert len(glucagon_seq) == 60 * 3  # full std fasta lines


########################################
# Base class for testing
########################################
class BaseTestSequence:
    #              1234567890
    protein_seq = "ELVISLIVES"
    #                   ____
    pep_seq = "LIVE"
    pep_start = 6
    pep_stop = 9

    def _assert_hash(self, my_set, item, items_before, items_after):
        assert len(my_set) == items_before
        my_set.add(item)
        assert len(my_set) == items_after

    # math should be index based, ie the Indexes should be added together, and any integer
    # should be assumed to be an "integer offset"
    def test__add__(self):
        # cls(1) has index = 0, pos = 1
        cls = self.test_class
        assert cls(1) + cls(1) == cls(1)  # because index 0 + 0 = 0
        assert cls(1) + 1 == cls(2)  # because index 0 + 1 = 1
        assert 1 + cls(1) == cls(2)  # because index 1 + 0 = 1
        #  assert cls(2) + (-1) == cls(1)  # because index -1 + 1 = 0
        #  assert -1 + cls(2) == cls(1)  # because index -1 + 1 = 0
        assert cls(2) + cls(2) == cls(3)  # because 1 + 1 = 2
        assert 10 + cls(2) == cls(12)  # because 10 + 1 = 11
        assert cls(10) + cls(1) == cls(10)  # because 9 + 0 = 9

    def test__sub__(self):
        # cls(1) has index = 0, pos = 1
        cls = self.test_class
        assert cls(1) - cls(1) == cls(1)  # because index 0 - 0 = 0
        assert cls(2) - 1 == cls(1)  # because index 1 - 1 = 0
        #  assert -1 + cls(1) == cls(2)  # because index 1 + 0 = 1
        assert cls(2) - cls(2) == cls(1)  # because 1 - 1 = 0
        assert cls(10) - cls(10) == cls(1)  # because 9 - 9 = 0
        assert 10 - cls(2) == cls(10)  # because 10 - 1 = 9
        assert cls(10) - cls(1) == cls(10)  # because 9 - 0 = 9

        # because 1 - 2 = -1 -> index = -1, pos = 0 -> invalid pos
        with pytest.raises(ValueError):
            cls(2) - 2
        with pytest.raises(ValueError):
            SequencePoint(2) - SequencePoint(3)
        with pytest.raises(ValueError):
            1 - SequencePoint(3)


########################################
# Tests for SequencePoint
########################################
class TestSequencePoint(BaseTestSequence):
    test_class = SequencePoint

    def _assert(self, start, stop, pep_seq, protein_seq):
        assert protein_seq[start.index] == pep_seq[0]
        assert protein_seq[start.index + 1] == pep_seq[1]

        assert protein_seq[stop.index] == pep_seq[-1]
        assert protein_seq[stop.index - 1] == pep_seq[-2]

        assert protein_seq[start.slice] == pep_seq[0]
        assert protein_seq[stop.slice] == pep_seq[-1]

        assert protein_seq[start.slice.start:stop.slice.stop] == pep_seq

    def test_init(self, glucagon_peptides, glucagon_seq):
        # simple tests
        start = SequencePoint(self.pep_start)
        stop = SequencePoint(self.pep_stop)
        self._assert(start, stop, self.pep_seq, self.protein_seq)

        with pytest.raises(ValueError):
            SequencePoint(0)
        with pytest.raises(ValueError):
            SequencePoint(-10)

        # ensure it can cast itself, just like int(int(1)) works
        assert SequencePoint(SequencePoint(1)) == SequencePoint(1)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            self._assert(SequencePoint(start), SequencePoint(stop), seq, glucagon_seq)

    def test_from_index(self, glucagon_peptides, glucagon_seq):
        # simple tests
        index = self.protein_seq.index(self.pep_seq)
        start = SequencePoint.from_index(index)
        stop = SequencePoint(index + len(self.pep_seq))
        self._assert(start, stop, self.pep_seq, self.protein_seq)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            index = glucagon_seq.index(seq)
            start = SequencePoint.from_index(index)
            # if len(seq) == 1, then start = stop, thus "+ len(seq) - 1"
            stop = SequencePoint.from_index(index + len(seq) - 1)
            self._assert(start, stop, seq, glucagon_seq)

    def test_comparisons(self):
        sl = SequencePoint(1)
        sl1 = SequencePoint(1)
        sl5 = SequencePoint(5)
        sl10 = SequencePoint(10)

        # equal / unequal
        assert sl == sl1
        assert sl == 1
        assert sl == "1"
        assert 1 == sl
        assert sl != sl5
        assert sl != 5
        assert 5 != sl

        assert not 5 < sl5
        assert not 5 > sl5
        assert 5 <= sl5
        assert 5 >= sl5

        with pytest.raises(TypeError):
            sl < "Wrong type!!"
        with pytest.raises(TypeError):
            "Wrong type!!" < sl
        with pytest.raises(TypeError):
            sl > "Wrong type!!"

        assert sl1 < sl5 < sl10 and sl1 <= sl5 <= sl10
        assert sl10 > sl5 > sl1 and sl10 >= sl5 >= sl1
        assert sl1 < 5 < sl10 and sl1 <= 5 <= sl10
        assert sl10 > 5 > sl1 and sl10 >= 5 >= sl1
        assert 1 < sl5 < 10 and 1 <= sl5 <= 10
        assert 10 > sl5 > 1 and 10 >= sl5 >= 1

        assert -1 < sl5 and -1 <= sl5

    def test__str__(self):
        assert str(SequencePoint(10)) == str(10)

    def test__repr__(self):
        assert repr(SequencePoint(10)) == "SequencePoint(10)"

    def test__hash__(self):
        hash(SequencePoint(1))
        my_set = set()
        self._assert_hash(my_set, SequencePoint(1), 0, 1)
        self._assert_hash(my_set, SequencePoint(1), 1, 1)
        self._assert_hash(my_set, SequencePoint(2), 1, 2)
        self._assert_hash(my_set, SequencePoint(2), 2, 2)
        self._assert_hash(my_set, SequencePoint(3), 2, 3)

    def test_immutability(self):
        s = SequencePoint(1)
        with pytest.raises(AttributeError):
            s.pos = 2
        with pytest.raises(AttributeError):
            s.index = 2

        assert s is SequencePoint(s)
        assert s is not SequencePoint(1)


########################################
# Tests for SequenceRange
########################################
class TestSequenceRange(BaseTestSequence):
    test_class = SequenceRange

    def _assert(self, p, pep_seq, protein_seq):
        assert protein_seq[p.index.start] == pep_seq[0]
        assert protein_seq[p.index.start + 1] == pep_seq[1]

        assert protein_seq[p.index.stop] == pep_seq[-1]
        assert protein_seq[p.index.stop - 1] == pep_seq[-2]

        assert protein_seq[p.slice.start:p.slice.stop] == pep_seq
        assert protein_seq[p.slice] == pep_seq

        assert len(pep_seq) == p.length == len(p)

    def test_init(self, glucagon_peptides, glucagon_seq):
        # simple tests
        p = SequenceRange(self.pep_start, self.pep_stop, protein_sequence=self.protein_seq)
        assert p.pos == (self.pep_start, self.pep_stop)
        assert self.pep_seq == p.seq
        self._assert(p, self.pep_seq, self.protein_seq)

        # pep_stop infered from len(pep_seq)
        assert p == SequenceRange(self.pep_start, seq=self.pep_seq)

        # has to be valid numbers!
        with pytest.raises(ValueError):
            SequenceRange(0, 0)
        with pytest.raises(ValueError):
            SequenceRange(-10, 10)
        with pytest.raises(ValueError):
            SequenceRange(10, -10)
        with pytest.raises(ValueError):
            SequenceRange(0, 10)
        with pytest.raises(ValueError):
            SequenceRange(10, 0)
        with pytest.raises(ValueError):
            SequenceRange(15, 10)

        # ensure it can cast itself, just like int(int(1)) == int(1)
        assert SequenceRange(SequenceRange(1)) == SequenceRange(1)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            p = SequenceRange(start, stop, seq)
            self._assert(p, seq, glucagon_seq)

        # if no stop, then the range is 1 amino acid
        p = SequenceRange(self.pep_start)
        assert p.pos.start == p.pos.stop
        assert self.protein_seq[p.slice] == 'L'

    def test_from_index_and_length(self, glucagon_peptides, glucagon_seq):
        # simple tests
        index = self.protein_seq.index(self.pep_seq)
        p = SequenceRange.from_index_and_length(index, len(self.pep_seq))
        self._assert(p, self.pep_seq, self.protein_seq)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            p = SequenceRange.from_index_and_length(glucagon_seq.index(seq), len(seq))
            self._assert(p, seq, glucagon_seq)

    def test_from_index(self, glucagon_peptides, glucagon_seq):
        pep_start_index = 5
        pep_stop_index = 8
        p_index = SequenceRange.from_index(pep_start_index, pep_stop_index)
        p_index2 = SequenceRange.from_index((pep_start_index, pep_stop_index))
        assert p_index == p_index2

        p = SequenceRange(self.pep_start, self.pep_stop)
        assert p == p_index
        assert self.pep_seq[0] == self.protein_seq[p_index.index.start]
        assert self.pep_seq[-1] == self.protein_seq[p_index.index.stop]

    def test_from_slices(self, glucagon_peptides, glucagon_seq):
        pep_start_slice = 5
        pep_stop_slice = 9
        p_slice = SequenceRange.from_slice(pep_start_slice, pep_stop_slice)
        p_slice2 = SequenceRange.from_slice(slice(pep_start_slice, pep_stop_slice))
        p = SequenceRange(self.pep_start, self.pep_stop)
        assert p == p_slice == p_slice2
        assert self.pep_seq == self.protein_seq[p_slice.slice.start:p_slice.slice.stop]
        assert self.pep_seq == self.protein_seq[p_slice.slice]

    def test_from_sequence(self, glucagon_peptides, glucagon_seq):
        expected = SequenceRange(self.pep_start, self.pep_stop)
        observed = SequenceRange.from_sequence(self.protein_seq, self.pep_seq)
        assert observed == expected
        with pytest.raises(IndexError):
            SequenceRange.from_sequence('PROTEINSEQ', "PEPTIDESEQ")

        for (start, stop, seq) in glucagon_peptides:
            p = SequenceRange.from_sequence(glucagon_seq, seq)
            self._assert(p, seq, glucagon_seq)

    def test___len__(self, glucagon_peptides, glucagon_seq):
        assert len(self.pep_seq) == len(SequenceRange(self.pep_start, self.pep_stop))
        for (start, stop, seq) in glucagon_peptides:
            assert len(SequenceRange(start, stop)) == len(seq)

    def test_comparisons(self):
        p = p00 = SequenceRange(self.pep_start, self.pep_stop)
        p2 = SequenceRange(self.pep_start, self.pep_stop)
        p01 = SequenceRange(self.pep_start, self.pep_stop + 1)
        p10 = SequenceRange(self.pep_start + 1, self.pep_stop)
        p_tuple = (self.pep_start, self.pep_stop)

        # equal / unequal
        assert p == p2
        assert p == p_tuple
        assert p_tuple == p
        assert p_tuple != p01
        assert p_tuple != 123
        with pytest.raises(TypeError):
            p < "Wrong type!!"
        with pytest.raises(TypeError):
            "Wrong type!!" < p
        with pytest.raises(TypeError):
            p > "Wrong type!!"

        assert not p < p and not p > p
        assert p <= p and p >= p and p == p

        # less/greater or equal
        # (x, y) < (x, y + 1), < (x+1, y)
        assert p00 < p01 < p10
        assert p10 > p01 > p00

        assert p_tuple < p01 < p10 and p_tuple <= p01 <= p10
        assert p10 > p01 > p_tuple and p10 >= p01 >= p_tuple

    def test__str__(self):
        assert str(SequenceRange(10, 20)) == "10:20"

    def test__repr__(self):
        assert repr(SequenceRange(10, 20)) == "SequenceRange(10, 20)"

    def test__add__extra(self):
        # SequenceRange(1,1) has index = (0, 0), pos = (1, 1)
        assert SequenceRange(1, 1) + SequenceRange(1) == SequenceRange(1, 1)
        assert SequenceRange(1, 1) + SequenceRange(1, 1) == SequenceRange(1, 1)
        assert SequenceRange(1, 5) + SequenceRange(1, 1) == SequenceRange(1, 5)
        assert SequenceRange(1, 5) + 1 == SequenceRange(2, 6)
        assert 1 + SequenceRange(1, 5) == SequenceRange(2, 6)
        assert -1 + SequenceRange(2, 6) == SequenceRange(1, 5)

    def test__sub__extra(self):
        # SequenceRange(1) has index = 0, pos = 1
        SequenceRange = self.test_class
        assert SequenceRange(1, 1) - SequenceRange(1) == SequenceRange(1, 1)
        assert SequenceRange(1, 5) - SequenceRange(1) == SequenceRange(1, 5)
        assert SequenceRange(2, 3) - 1 == SequenceRange(1, 2)
        assert -1 + SequenceRange(3, 4) == SequenceRange(2, 3)
        assert SequenceRange(2, 5) - SequenceRange(2, 5) == SequenceRange(1)
        assert SequenceRange(2, 7) - SequenceRange(2, 5) == SequenceRange(1, 3)
        assert SequenceRange(10, 20) - SequenceRange(10, 20) == SequenceRange(1)

        # because 1 - 2 = -1 -> index = -1, pos = 0 -> invalid pos
        with pytest.raises(ValueError):
            SequenceRange(2, 20) - 2
        with pytest.raises(ValueError):
            SequenceRange(2, 20) - SequenceRange(3)

    def test__iter__(self):
        sr = SequenceRange(5, 10)
        sr_points = list(sr)  # should be equivalent to list(sr.__iter__())
        assert sr.length == len(sr_points)
        assert sr_points[0].pos == sr.pos.start
        assert sr_points[-1].pos == sr.pos.stop

        assert list(SequenceRange(5, 5))[0] == SequencePoint(5)

    def test__getitem__(self):
        p = SequenceRange(5, 10)
        assert p[0] == p[:1][0] == p[:-1][0] == 5
        assert p[1] == p[1:][0] == p[-1:][0] == 10
        assert p[:] == p
        assert type(p[:]) != type(p)

    def test__hash__(self):
        hash(SequenceRange(1, 2))
        my_set = set()
        self._assert_hash(my_set, SequenceRange(1, 1), 0, 1)
        self._assert_hash(my_set, SequenceRange(1, 1), 1, 1)
        self._assert_hash(my_set, SequenceRange(1, 2), 1, 2)
        self._assert_hash(my_set, SequenceRange(2, 2), 2, 3)
        self._assert_hash(my_set, SequenceRange(1, 2), 3, 3)

    def test_immutability(self):
        s = SequenceRange(1, 2)
        with pytest.raises(AttributeError):
            s.pos = (1, 2)
        with pytest.raises(AttributeError):
            s.index = (1, 2)
        with pytest.raises(AttributeError):
            s.slice = (1, 2)

        assert s is SequenceRange(s)
        assert s is not SequenceRange(1, 2)

    def _in(self, cls, arg, peptide):
        assert arg in peptide
        assert cls(arg) in peptide

    def _not_in(self, cls, arg, peptide):
        assert arg not in peptide
        assert cls(arg) not in peptide

    def test__contains__(self):
        peptide = SequenceRange(5, 20)

        self._in(SequencePoint, 5, peptide)
        self._in(SequencePoint, 10, peptide)
        self._in(SequencePoint, 20, peptide)
        self._not_in(SequencePoint, 4, peptide)
        self._not_in(SequencePoint, 21, peptide)

        self._in(SequenceRange, (5, 10), peptide)
        self._in(SequenceRange, (10, 15), peptide)
        self._in(SequenceRange, (15, 20), peptide)
        self._not_in(SequenceRange, (1, 5), peptide)
        self._not_in(SequenceRange, (4, 11), peptide)
        self._not_in(SequenceRange, (10, 21), peptide)

    def test_contains(self):
        #          0        10        20
        #          0123456789012345678901 #
        #  self  = -----ELVISLIVES        #
        #  item1 = ----------L----        <--- part=all -> True,  part=any -> True
        #  item2 = ----------LIVE-        <--- part=all -> True,  part=any -> True
        #  item3 = ----------LIVESANDDIES <--- part=all -> False, part=any -> True
        #  item4 = ELVENELVISLIVESANDDIES <--- part=all -> False, part=any -> True
        #  item5 = ------------------D--- <--- part=all -> False, part=any -> False
        #  item6 = ------------------DIES <--- part=all -> False, part=any -> False

        _self = SequenceRange(4, 14)
        item1 = SequencePoint(10)
        item2 = SequenceRange(10, 13)
        item3 = SequenceRange(10, 21)
        item4 = SequenceRange(1, 21)
        item5 = SequenceRange(18, 21)
        item6 = SequenceRange(18, 21)

        assert _self.contains(item1, part=all)
        assert _self.contains(item2, part=all)
        assert not _self.contains(item3, part=all)
        assert not _self.contains(item4, part=all)
        assert not _self.contains(item5, part=all)
        assert not _self.contains(item6, part=all)

        assert _self.contains(item1, part=any)
        assert _self.contains(item2, part=any)
        assert _self.contains(item3, part=any)
        assert _self.contains(item4, part=any)
        assert not _self.contains(item5, part=any)
        assert not _self.contains(item6, part=any)


class TestInteroperability:
    def test_conversion(self):
        assert SequenceRange(1, 2) == SequenceRange(SequencePoint(1), SequencePoint(2))
        assert SequenceRange(1, 2) == SequenceRange(SequencePoint(1), 2)
        assert SequenceRange(1, 2) == SequenceRange(1, SequencePoint(2))

    def test__add__(self):
        assert SequenceRange(2, 3) + SequencePoint(2) == SequenceRange(3, 4)
        assert SequencePoint(2) + SequenceRange(2, 3) == SequenceRange(3, 4)
        assert SequenceRange(2, 3) + SequencePoint(2) + 5 == SequenceRange(8, 9)
        assert 5 + SequenceRange(2, 3) + SequencePoint(2) == SequenceRange(8, 9)
        assert SequenceRange(2, 3) + 5 + SequencePoint(2) == SequenceRange(8, 9)

    def test__sub__(self):
        assert SequenceRange(5, 20) - SequencePoint(3) == SequenceRange(3, 18)
        with pytest.raises(ValueError):
            # 20 - 5, 20 - 10 -> 15, 10 = makes no sense!!
            SequencePoint(20) - SequenceRange(5, 10)

        assert SequenceRange(10, 15) - SequencePoint(5) == SequenceRange(6, 11)
