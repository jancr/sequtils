
# core imports
import os

# 3rd party imports
import pytest

# local imports
from sequtils import ProteinLocation

TEST_FOLDER = os.path.abspath(os.path.dirname(__file__))
TEST_FILES_FOLDER = os.path.abspath(os.path.join(TEST_FOLDER, 'test_files'))

########################################
# fixtures
# move to conftest when project becomes larger!
########################################
@pytest.fixture
def glucagon_seq():
    with open(os.path.join(TEST_FILES_FOLDER, 'glucagon.fasta')) as f:
        return "".join(f.readlines()[1:])


@pytest.fixture
def glucagon_peptides():
    """ all peptides from glucagon as ((start, stop, seq), ..) """
    peptides = []
    with open(os.path.join(TEST_FILES_FOLDER, 'glucagon_peptides.fasta')) as f:
        while f.tell() == os.fstat(f.fileno()).st_size:
            start, stop = map(int, f.readline().split('|')[-1].split('-'))
            seq = f.readline().strip()
            peptides.append((start, stop, seq))
    return tuple(peptides)

########################################
# Tests for ProteinLoaction
########################################
class TestProteinLocaton():
    #              1234567890
    protein_seq = "ELVISLIVES"
    pep_seq     =   "VISL"
    pep_start   =    3
    pep_stop     =       6

    def _assert(self, p, pep_seq, protein_seq):
        assert protein_seq[p.index.start] == pep_seq[0]
        assert protein_seq[p.index.start+1] == pep_seq[1]

        assert protein_seq[p.index.stop] == pep_seq[-1]
        assert protein_seq[p.index.stop-1] == pep_seq[-2]

        assert protein_seq[p.slice.start:p.slice.stop] == pep_seq
        assert protein_seq[p.slice] == pep_seq

        assert len(pep_seq) == p.length == len(p)

    def test_init(self, glucagon_peptides, glucagon_seq):
        # simple tests
        p = ProteinLocation(self.pep_start, self.pep_stop)
        assert (p.pos.start, p.pos.stop) == (self.pep_start, self.pep_stop)
        self._assert(p, self.pep_seq, self.protein_seq)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            p = ProteinLocation(start, stop, seq)
            self._assert(p, seq, glucagon_seq)

        # if no stop, then the range is 1 amino acid
        p = ProteinLocation(self.pep_start)
        assert p.pos.start == p.pos.stop
        assert self.protein_seq[p.slice] == 'V'

    def test_from_index_and_length(self, glucagon_peptides, glucagon_seq):
        # simple tests
        index = self.protein_seq.index(self.pep_seq)
        p = ProteinLocation.from_index_and_length(index, len(self.pep_seq))
        self._assert(p, self.pep_seq, self.protein_seq)

        # all peptides
        for (start, stop, seq) in glucagon_peptides:
            p = ProteinLocation.from_index_and_length(glucagon_seq.index(seq), len(seq))
            self._assert(p, seq, glucagon_seq)

    def test_from_indexes(self, glucagon_peptides, glucagon_seq):
        pep_start_index = 2
        pep_stop_index = 5
        p_index = ProteinLocation.from_indexes(pep_start_index, pep_stop_index)
        p = ProteinLocation(self.pep_start, self.pep_stop)
        assert p == p_index
        assert self.pep_seq[0] == self.protein_seq[p_index.index.start]
        assert self.pep_seq[-1] == self.protein_seq[p_index.index.stop]

    def test_from_slices(self, glucagon_peptides, glucagon_seq):
        pep_start_slice = 2
        pep_stop_slice = 6
        p_slice = ProteinLocation.from_slice(pep_start_slice, pep_stop_slice)
        p_slice2 = ProteinLocation.from_slice(slice(pep_start_slice, pep_stop_slice))
        p = ProteinLocation(self.pep_start, self.pep_stop)
        assert p == p_slice == p_slice2
        assert self.pep_seq == self.protein_seq[p_slice.slice.start:p_slice.slice.stop]
        assert self.pep_seq == self.protein_seq[p_slice.slice]

    def test_from_sequence(self, glucagon_peptides, glucagon_seq):
        expected = ProteinLocation(self.pep_start, self.pep_stop)
        observed = ProteinLocation.from_sequence(self.protein_seq, self.pep_seq)
        assert observed == expected
        with pytest.raises(IndexError):
            ProteinLocation.from_sequence('PROTEINSEQ', "PEPTIDESEQ")

        for (start, stop, seq) in glucagon_peptides:
            p = ProteinLocation.from_sequence(seq, glucagon_seq)
            self._assert(p, seq, glucagon_seq)

    def test___len__(self, glucagon_peptides, glucagon_seq):
        assert len(self.pep_seq) == len(ProteinLocation(self.pep_start, self.pep_stop))
        for (start, stop, seq) in glucagon_peptides:
            assert len(ProteinLocation(start, stop)) == len(seq)

    def test_comparisons(self):
        p = p00 = ProteinLocation(self.pep_start, self.pep_stop)
        p2 = ProteinLocation(self.pep_start, self.pep_stop)
        p01 = ProteinLocation(self.pep_start, self.pep_stop + 1)
        p10 = ProteinLocation(self.pep_start + 1, self.pep_stop)
        p_tuple = (self.pep_start, self.pep_stop)

        # equal / unequal
        assert p == p2
        assert p == p_tuple
        assert p_tuple == p
        assert p_tuple != p01

        # less/greater or equal
        # (x, y) < (x, y + 1), < (x+1, y)
        assert p00 < p01 < p10
        assert p10 > p01 > p00

        assert p_tuple < p01 < p10
        assert p10 > p01 > p_tuple

    def test__str__(self):
        assert str(ProteinLocation(10, 20)) == "10:20"

    def test__repr__(self):
        assert repr(ProteinLocation(10, 20)) == "ProteinLocation(10, 20)"

    def test__add__(self):
        assert ProteinLocation(2, 2) + ProteinLocation(3, 3) == ProteinLocation(5, 5)
        assert ProteinLocation(2, 2) + 3 == ProteinLocation(5, 5)
        assert 3 + ProteinLocation(2, 2) == ProteinLocation(5, 5)
        assert 3 + ProteinLocation(2, 5) == ProteinLocation(5, 8)

    def test__sub__(self):
        assert ProteinLocation(2, 2) - 1 == ProteinLocation(1, 1)
        assert ProteinLocation(2, 5) - 1 == ProteinLocation(1, 4)
        assert 2 - ProteinLocation(1, 1) == ProteinLocation(1, 1)

        with pytest.raises(ValueError):
            ProteinLocation(2, 2) - ProteinLocation(3, 3)
        with pytest.raises(ValueError):
            ProteinLocation(2, 2) - ProteinLocation(1, 3)
        with pytest.raises(ValueError):
            ProteinLocation(2, 2) - 3
        with pytest.raises(ValueError):
            2 - ProteinLocation(3, 3)


