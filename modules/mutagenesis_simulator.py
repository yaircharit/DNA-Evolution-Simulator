from pathlib import Path
import random

from modules.config_loader import sim_config
from modules.dna_sequence import dna_sequence
from modules.simulation_logger import SimulationLogger, log_sequence_run, log_simulation_run


class mutagenesis_simulator:
    sequences: list[dna_sequence]

    def __init__(self, iterations: int = sim_config.iterations,
                 sequences_num: int = sim_config.sequences_num,
                 random_seed: int = sim_config.random_seed,
                 generate_sequences: bool = sim_config.generate_sequences,
                 input_file: str = sim_config.input_file) -> None:

        self.iterations = iterations
        self.sequences_num = sequences_num
        self.random_seed = random_seed if random_seed is not None else random.randint(0, 2**32 - 1)
        self.generate_sequences = generate_sequences
        self.input_file = input_file

        self.output_dir = Path(sim_config.output_dir)
        if not self.output_dir.is_absolute():
            self.output_dir = Path(__file__).resolve().parent.parent / self.output_dir
        self.logger = SimulationLogger(self.output_dir)
        self.output_dir = self.logger.output_dir

        if self.generate_sequences:
            random.seed(self.random_seed)
            self.sequences = [dna_sequence() for _ in range(sequences_num)]
        else:
            with open(self.input_file, "r", encoding="utf-8") as f:
                self.sequences = [dna_sequence(line.strip()) for line in f]

    @log_simulation_run
    def run_simulation(self) -> list[list[dict]]:
        sequence_history = []
        for seq in self.sequences:
            sequence_history.append(self.run_sequence(seq))

        return sequence_history

    @log_sequence_run
    def run_sequence(self, seq: dna_sequence) -> list[dict]:
        history: list[dict] = [{
            "iteration": -1,
            "mutation": None,
            "sequence": str(seq),
            "protein": ''.join(seq.protein),
            "stats": {
                "nucleotides_diff": seq.nucleotides_diff,
                "amino_acids_diff": seq.amino_acids_diff,
                "mutations_count": seq.mutations_count,
                "rejected_mutations_count": seq.rejected_mutations_count,
            },
        }]

        seq.print_original_sequence()
        for iteration in range(self.iterations):
            mutation = seq.mutate()
            seq.calculate_diff()

            stats = {
                "nucleotides_diff": seq.nucleotides_diff,
                "amino_acids_diff": seq.amino_acids_diff,
                "mutations_count": seq.mutations_count,
                "rejected_mutations_count": seq.rejected_mutations_count,
            }

            self.logger.log_iteration(
                seq_id=seq.id,
                iteration=iteration,
                mutation=mutation,
                sequence=str(seq),
                protein=''.join(seq.protein),
                stats=stats,
            )
            history.append({
                "iteration": iteration,
                "mutation": mutation,
                "sequence": str(seq),
                "protein": ''.join(seq.protein),
                "stats": stats,
            })

        seq.print_sequence_stats()
        print()
        return history
