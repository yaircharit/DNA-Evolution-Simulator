import random

from modules.config_loader import aa_config
from modules.dna_sequence import dna_sequence


def test_random_initialization_generates_valid_sequence():
    seq = dna_sequence(length=9)

    assert len(seq.sequence) == 9
    assert seq.length == 9
    assert seq.original_sequence == "".join(seq.sequence)
    assert seq.original_protein == "".join(seq.protein)

    codons = ["".join(seq.sequence[index : index + 3]) for index in range(0, len(seq.sequence), 3)]
    assert codons[0] == "ATG"
    assert codons[-1] in aa_config.stop_codons
    assert all(codon not in aa_config.stop_codons for codon in codons[:-1])


def test_sequence_initialization_tracks_protein_and_ranges():
    seq = dna_sequence("ATGAAATAG")

    assert seq.original_sequence == "ATGAAATAG"
    assert "".join(seq.sequence) == "ATGAAATAG"
    assert "".join(seq.protein) == "MK*"
    assert seq.mutation_start_index == 3
    assert seq.mutation_end_index == 5


def test_mutate_changes_one_nucleotide_and_returns_summary(monkeypatch):
    seq = dna_sequence("ATGAAATAG")

    monkeypatch.setattr(random, "random", lambda: 0.0)
    monkeypatch.setattr(random, "randint", lambda start, end: start)
    monkeypatch.setattr(
        dna_sequence,
        "get_random_nucliotide",
        staticmethod(lambda start_nucleotide=None, k=1: "C"),
    )

    result = seq.mutate()

    assert result == "3;A;C"
    assert seq.sequence[3] == "C"
    assert "".join(seq.sequence) == "ATGCAATAG"


def test_internal_stop_codon_protection_prevents_new_internal_stop(monkeypatch):
    seq = dna_sequence("ATGTGGGGG")
    seq.internal_stop_codons = False

    values = iter(["A", "C"])
    monkeypatch.setattr(
        dna_sequence,
        "get_random_nucliotide",
        staticmethod(lambda start_nucleotide=None, k=1: next(values)),
    )

    seq.mutate_nucleotide(5)

    codons = ["".join(seq.sequence[index : index + 3]) for index in range(0, len(seq.sequence), 3)]
    assert codons[0] == "ATG"
    assert codons[1] not in aa_config.stop_codons
    assert codons[2] not in aa_config.stop_codons


def test_external_codons_protection_skips_start_and_stop_codons():
    protected = dna_sequence("ATGAAATAG", external_codon_protection=True)
    unprotected = dna_sequence("ATGAAATAG", external_codon_protection=False)

    assert protected.mutation_start_index == 3
    assert protected.mutation_end_index == protected.length - 4
    assert unprotected.mutation_start_index == 0
    assert unprotected.mutation_end_index == unprotected.length - 1


def test_calculate_diff_counts_nucleotide_and_amino_acid_changes():
    seq = dna_sequence("ATGAAATAG")
    seq.sequence = list("ATGCAATAG")
    seq.calculate_diff()

    assert seq.nucleotides_diff == 1
    assert seq.amino_acids_diff == 1
    assert "".join(seq.protein) == "MQ*"
