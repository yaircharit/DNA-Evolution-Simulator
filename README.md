# DNA Evolution Simulator

A Python simulation for evolving DNA sequences over many generations, tracking how nucleotide mutations affect the translated protein sequence.

This project generates random coding sequences, mutates them according to configured nucleotide transition probabilities, and records how each mutation changes the nucleotide string and the resulting amino acid sequence. It also saves summary CSV files and chart images for analysis.

## What the project does

- Creates DNA sequences that begin with a start codon and end with a stop codon.
- Uses a codon table to translate DNA into amino acids.
- Applies mutation events while protecting start/stop codons by default.
- Tracks:
  - nucleotide differences from the original sequence,
  - amino acid differences from the original protein,
  - mutation counts,
  - rejected mutation counts.
- Logs per-sequence and per-iteration results to the `results/` folder.
- Produces plots for sequence evolution and simulation summaries using Matplotlib.

## Project structure

- `main.py` – entry point for running the simulator.
- `modules/` – implementation of the simulator logic.
  - `dna_sequence.py` – DNA sequence generation, mutation, translation, and stats.
  - `mutagenesis_simulator.py` – orchestration of simulation runs.
  - `config_loader.py` – loads configuration from the `configs/` folder.
  - `data_analyzer.py` – generates summary plots from CSV data.
  - `simulation_logger.py` – writes simulation logs and CSV outputs.
- `configs/` – simulator configuration files.
- `results/` – generated output from simulations.
- `tests/` – unit tests for sequence mutation, logging, and analysis.

## Configuration

The simulator reads its settings from JSON files in the `configs/` directory:

- `simulator_config.json` – high-level run settings such as iteration count, number of sequences, output path, and whether to generate sequences automatically.
- `nucleotide_config.json` – nucleotide composition and mutation transition probabilities.
- `amino_acids_config.json` – codon table, start/stop codons, and translation-related rules.

You can change the following values to tune behavior:

- number of generations to simulate,
- number of sequences per run,
- random seed,
- sequence length,
- mutation probabilities,
- whether to protect start/stop codons,
- whether to generate random sequences or load from an input file.

## Running the simulator

From the project root:

```bash
python main.py
```

If you are using the included virtual environment on Windows:

```powershell
.\.venv\Scripts\python.exe main.py
```

The simulator will run the configured simulation and write output files under the `results/` directory.

## Output files

Each simulation run creates output such as:

- sequence CSV logs,
- summary CSV files,
- PNG plots of mutation trajectories,
- summary-analysis charts.

These files are useful for understanding how mutation pressure changes DNA and protein composition over time.

## Example use cases

- Testing how likely mutations are to change the protein sequence.
- Studying the effect of different nucleotide transition matrices.
- Simulating evolutionary drift under constrained codon protection.
- Comparing mutation outcomes across many random sequences.

## Requirements

This project uses Python 3 and depends on packages such as:

- `matplotlib`
- `numpy`

The repository includes a local virtual environment (`.venv`) that already contains the required dependencies for this project.

## License

This project does not currently include a formal license file. Check the repository status and project owner for licensing terms before publishing or redistributing the code.
