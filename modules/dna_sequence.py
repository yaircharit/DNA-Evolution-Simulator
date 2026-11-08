import random
from modules.config_loader import nuc_config, aa_config

seq_counter = 0  # Global counter for sequence IDs
class dna_sequence:
    id: int
    original_sequence: str
    sequence: list[str]
    length: int
    
    amino_acids_diff: int
    nucleotides_diff: int
    mutations_count: int
    rejected_mutations_count: int

    mutation_start_index: int = 0
    mutation_end_index: int
    internal_stop_codons: bool
    external_codons_protection: bool

    mutation_rate = 0.01  # Default mutation rate, can be adjusted as needed


    @property
    def protein(self) -> list[str]:
        """Convert the current nucleotide sequence to an amino acid sequence."""
        return self.to_amino_acids("".join(self.sequence))

    def __str__(self) -> str:
        """Return the current nucleotide sequence as a string."""
        return "".join(self.sequence)

    def __repr__(self) -> str:
        """Return a compact representation of the sequence and its translated protein."""
        return f"dna_sequence(id={self.id}, nucleotides='{self}', protein='{''.join(self.protein)}')"

    def __init__(self, sequence: str | None = None, 
                 length: int = nuc_config.sequence_length,
                 internal_stop_codons: bool | None = None,
                 external_codon_protection: bool | None = None) -> None:

        self.internal_stop_codons = aa_config.internal_stop_codons if internal_stop_codons is None else internal_stop_codons
        self.external_codons_protection = aa_config.external_codons_protection if external_codon_protection is None else external_codon_protection

        if sequence is None:
            self.initialize_sequence(self.generate_random(length))
        else:
            self.initialize_sequence(sequence)

    def initialize_sequence(self, sequence: str) -> None:
        """Initialize the sequence with a given string."""
        global seq_counter

        self.id = seq_counter
        seq_counter += 1
        self.original_sequence = sequence
        self.original_protein = ''.join(self.to_amino_acids(sequence))
        self.length = len(sequence)
        self.nucleotides_diff = 0
        self.amino_acids_diff = 0
        self.mutations_count = 0
        self.rejected_mutations_count = 0
        self.sequence = list(sequence)  # Convert the string sequence to a list for mutability

        if self.external_codons_protection:
            self.mutation_start_index = 3  # Skip the start codon
            self.mutation_end_index = self.length - 4  # Skip the stop codon (3+1)
        else:
            self.mutation_start_index = 0
            self.mutation_end_index = self.length - 1

    @staticmethod
    def generate_random(length: int) -> str:
        """Generate a nucleotide sequence using the configured transition matrix."""
        nucleotides = nuc_config.nucleotides
        nucs = list(nucleotides.keys())
        nucs_weights = list(nucleotides.values())

        seq = list(aa_config.start_codons[0])

        # Generate a nucleotide based on the transition matrix
        mid = random.choices(
            population=nucs,
            weights=nucs_weights,
            k=length - 6,
        )

        dna_sequence.remove_stop_codons(mid)  # Ensure no internal stop codons are present

        seq.extend(mid)
        stop_codon = aa_config.stop_codons[random.randint(0, len(aa_config.stop_codons) - 1)]
        seq.extend(list(stop_codon))

        return "".join(seq)
    
    def mutate(self, mutation_rate: float = 1.0) -> str | None:

        if random.random() < mutation_rate:
            rand_pos = random.randint(self.mutation_start_index, self.mutation_end_index)
            return self.mutate_nucleotide(rand_pos)

        return None  # No mutation occurred


    def mutate_nucleotide(self, position: int) -> str | None:
        """Mutate a specific nucleotide at the given position."""
        if 0 <= position < len(self.sequence):
            selected_nuc = self.sequence[position]
            candidate = self.get_random_nucliotide(selected_nuc)
            self.sequence[position] = candidate

            if not self.internal_stop_codons:
                # Ensure no internal stop codons are present after mutation.
                in_codon_position = position % 3
                codon_start_index = position - in_codon_position
                codon = self.sequence[codon_start_index:codon_start_index + 3]

                while "".join(codon) in aa_config.stop_codons:
                    replacement = self.get_random_nucliotide(selected_nuc)
                    codon[in_codon_position] = replacement
                    self.sequence[codon_start_index:codon_start_index + 3] = codon
                    self.rejected_mutations_count += 1

            self.mutations_count += 1
            return f"{position};{selected_nuc};{self.sequence[position]}"

        return None

    @staticmethod
    def get_random_nucliotide(start_nucleotide: str = None, k: int = 1) -> str:
        src = nuc_config.nucleotides if start_nucleotide is None else nuc_config.transition_matrix[start_nucleotide]
        return random.choices(
                population=list(src.keys()),
                weights=list(src.values()),
                k=1
            )[0]

    @staticmethod
    def remove_stop_codons(sequence: list) -> None:
        """Remove internal stop codons from the sequence."""
        for i in range(0, len(sequence), 3):
            codon = ''.join(sequence[i:i + 3])
            if codon in aa_config.stop_codons:
                print(f"Internal stop codon found at position {i}: {codon}")
                replacement = dna_sequence.get_random_nucliotide(sequence[i])
                sequence[i] = replacement
                # Re-check if the codon remains a stop codon.
                while ''.join(sequence[i:i + 3]) in aa_config.stop_codons:
                    sequence[i] = dna_sequence.get_random_nucliotide(sequence[i])

    def calculate_diff(self) -> None:
        """Calculate the number of nucleotides that have changed from the original sequence."""
        self.nucleotides_diff = sum(
            1 for original, current in zip(self.original_sequence, self.sequence) if original != current
        )
        self.amino_acids_diff = len(
            [True for orig, curr in zip(self.original_protein, self.protein) if orig != curr]
        )

    @staticmethod
    def to_amino_acids(sequence: str) -> list[str]:
        """Convert a nucleotide sequence to an amino acid sequence using the codon table."""
        codon_table = aa_config.codon_table
        amino_acids = []

        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i + 3]
            amino_acid = codon_table.get(codon, '?')  # Use '?' for unknown codons
            amino_acids.append(amino_acid)

        return amino_acids


    def print_original_sequence(self) -> None:
        """Print the original nucleotide sequence and its corresponding amino acid sequence."""
        print(f"Original #{self.id} length: {self.length}")
        print(f"> {self.original_sequence}")
        print(f"> {self.original_protein}")


    def print_sequence_stats(self) -> None:
        """Print the current nucleotide sequence and its corresponding amino acid sequence."""
        print(f"Current #{self.id} length: {self.length} nucs diff: {self.nucleotides_diff}, aa diff: {self.amino_acids_diff}")
        print(f"> {''.join(self.sequence)}")
        print(f"> {''.join(self.protein)}")

