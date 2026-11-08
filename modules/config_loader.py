import json
from pathlib import Path
from typing import Any


class base_config:
    _config_name = "config.json"

    def __init__(self, configs_dir: str | Path | None = None) -> None:
        if configs_dir is None:
            self._configs_dir = Path(__file__).resolve().parent.parent / "configs"
        else:
            self._configs_dir = Path(configs_dir)

        self._config_path = self._configs_dir / self._config_name
        self.load()  # Load the configuration upon initialization

    def load(self) -> None:
        """Load the configuration from the file."""
        if not self._config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self._config_path}")
        with self._config_path.open("r", encoding="utf-8") as file:
            data = json.load(file)
            for key, value in data.items():
                setattr(self, key, value)


class nucleotide_config(base_config):
    _config_name = "nucleotide_config.json"

    sequence_length: int = 9
    nucleotides: dict[str, float] = {}
    transition_matrix: dict[str, dict[str, float]] = {}


class amino_acid_config(base_config):
    _config_name = "amino_acids_config.json"

    start_codons: list[str] = []
    stop_codons: list[str] = []
    amino_acids: dict[str, str] = {}
    codon_table: dict[str, str] = {}
    internal_stop_codons: bool = False
    external_codons_protection: bool = True

class simulation_config(base_config):
    _config_name = "simulator_config.json"

    iterations: int 
    sequences_num: int
    output_dir: str = "results"
    random_seed: int
    generate_sequences: bool = True
    input_file: str 

def load_configs(config_dir: str | Path | None = None) -> None:
    global nuc_config, aa_config, sim_config
    nuc_config = nucleotide_config(config_dir)
    aa_config = amino_acid_config(config_dir)
    sim_config = simulation_config(config_dir)

loaded_configs = load_configs()
