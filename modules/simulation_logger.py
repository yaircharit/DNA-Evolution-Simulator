import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import matplotlib

from modules.config_loader import nuc_config
from modules.data_analyzer import SimulationDataAnalyzer
from modules.dna_sequence import dna_sequence

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class SimulationLogger:
    """Persist and visualize simulation data for each sequence and the whole run."""

    def __init__(self, output_dir: str | Path = "results") -> None:
        self.base_output_dir = Path(output_dir)
        if not self.base_output_dir.is_absolute():
            self.base_output_dir = Path(__file__).resolve().parent.parent / self.base_output_dir
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        run_name = f"simulation_{timestamp}"
        self.output_dir = self.base_output_dir / run_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.log_path = self.output_dir / "simulation_log.jsonl"
        self.settings_path = self.output_dir / "simulation_settings.json"
        self.sequence_csv_path = self.output_dir / "sequence_summary.csv"
        self.simulation_summary_csv_path = self.output_dir / "simulation_summary.csv"
        self.sequences_dir = self.output_dir / "sequences"
        self.sequences_dir.mkdir(parents=True, exist_ok=True)
        self.sequence_histories: list[dict[str, Any]] = []

    def save_settings(self, params: dict[str, Any]) -> None:
        payload = {
            **params,
            "sequence_length": nuc_config.sequence_length,
            "nucleotide_initial_distribution": nuc_config.nucleotides,
            "nucleotide_transition_matrix": nuc_config.transition_matrix,
        }
        self.settings_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    def log_start(self, params: dict[str, Any]) -> None:
        self.save_settings(params)
        self._write_event({"event": "simulation_start", "parameters": params})

    def log_sequence_start(self, seq_id: int, seq: str, protein: str, length: int) -> None:
        self._write_event(
            {
                "event": "sequence_start",
                "sequence_id": seq_id,
                "sequence": seq,
                "protein": protein,
                "length": length,
            },
            path=self.output_dir / "simulation_log.jsonl",
        )

    def log_iteration(self, seq_id: int, iteration: int, mutation: str | None, sequence: str, protein: str, stats: dict[str, Any]) -> None:
        self._write_event(
            {
                "event": "mutation",
                "sequence_id": seq_id,
                "iteration": iteration,
                "mutation": mutation,
                "sequence": sequence,
                "protein": protein,
                "stats": stats,
            },
            path=self.output_dir / "simulation_log.jsonl",
        )

    def log_sequence_end(self, seq_id: int, history: list[dict[str, Any]], final_stats: dict[str, Any]) -> None:
        payload = {
            "sequence_id": seq_id,
            "history": history,
            "final_stats": final_stats,
        }
        self.sequence_histories.append(payload)
        self.save_sequence_history_csv(seq_id, history)
        SimulationDataAnalyzer.plot_sequence_history(seq_id=seq_id, history=history, output_dir=self.sequences_dir)
        self._write_event(
            {
                "event": "sequence_end",
                "sequence_id": seq_id,
                "final_stats": final_stats,
            },
            path=self.output_dir / "simulation_log.jsonl",
        )

    def save_sequence_history_csv(self, seq_id: int, history: list[dict[str, Any]]) -> Path:
        self.sequences_dir.mkdir(parents=True, exist_ok=True)

        fieldnames = [
            "ID",
            "mutation number",
            "mutation",
            "nucleotides diff",
            "amino acids diff",
            "rejected mutations count",
            "sequence",
            "amino acids sequence",
        ]
        rows: list[dict[str, Any]] = []

        for step in history:
            sequence = step.get("sequence") or step.get("sequence_after") or step.get("sequence_before") or ""
            stats = step.get("stats", {}) if isinstance(step.get("stats"), dict) else {}
            mutation = step.get("mutation")
            mutation_number = stats.get("mutations_count", step.get("mutation_count", 0))
            rejected_mutations = stats.get("rejected_mutations_count", step.get("rejected_mutations", 0))
            rows.append(
                {
                    "ID": seq_id,
                    "mutation number": mutation_number,
                    "mutation": mutation or "",
                    "nucleotides diff": stats.get("nucleotides_diff", 0),
                    "amino acids diff": stats.get("amino_acids_diff", 0),
                    "rejected mutations count": rejected_mutations,
                    "sequence": sequence,
                    "amino acids sequence": "".join(dna_sequence.to_amino_acids(sequence)),
                }
            )

        output_path = self.sequences_dir / f"sequence_{seq_id}.csv"
        with output_path.open("w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        SimulationDataAnalyzer.plot_sequence_history(seq_id=seq_id, history=history, output_dir=self.sequences_dir)
        return output_path


    def build_sequence_row(self, entry: dict[str, Any]) -> dict[str, Any]:
        history = entry.get("history", [])
        final_stats = entry.get("final_stats", {})

        if not history:
            initial_sequence = ""
            final_sequence = ""
        else:
            first_entry = history[0]
            last_entry = history[-1]
            initial_sequence = first_entry.get("sequence_before") or first_entry.get("sequence") or first_entry.get("original_sequence") or ""
            final_sequence = last_entry.get("sequence_after") or last_entry.get("sequence") or last_entry.get("final_sequence") or initial_sequence

        initial_protein = "".join(dna_sequence.to_amino_acids(initial_sequence))
        final_protein = "".join(dna_sequence.to_amino_acids(final_sequence))
        nucleotide_diff = final_stats.get("nucleotides_diff")
        amino_acid_diff = final_stats.get("amino_acids_diff")
        if nucleotide_diff is None:
            nucleotide_diff = sum(1 for orig, curr in zip(initial_sequence, final_sequence) if orig != curr)
        if amino_acid_diff is None:
            amino_acid_diff = sum(1 for orig, curr in zip(initial_protein, final_protein) if orig != curr)

        return {
            "ID": entry.get("sequence_id", 0),
            "initial_sequence": initial_sequence,
            "final_sequence": final_sequence,
            "initial_protein": initial_protein,
            "final_protein": final_protein,
            "nucleotide_difference_count": nucleotide_diff,
            "amino_acids_difference_count": amino_acid_diff,
            "differences_ratio": round((amino_acid_diff / nucleotide_diff), 4) if nucleotide_diff else 0.0,
            "rejected_mutation_count": final_stats.get("rejected_mutations_count", 0),
        }

    @staticmethod
    def _nucleotide_change_ratios(initial_sequence: str, final_sequence: str) -> dict[str, float]:
        change_counts = {base: 0 for base in "ATCG"}
        total_changes = 0

        for original, current in zip(initial_sequence, final_sequence):
            if original != current:
                total_changes += 1
                if original in change_counts:
                    change_counts[original] += 1

        if total_changes == 0:
            return {base: 0.0 for base in "ATCG"}

        return {base: round(change_counts[base] / total_changes, 4) for base in "ATCG"}

    def save_sequence_csv(self) -> Path:
        rows: list[dict[str, Any]] = []
        fieldnames = [
            "sequence ID",
            "mutations number",
            "nucleotides difference",
            "amino acids difference",
            "amino acids difference to nucleotides difference ratio",
            "initial sequence",
            "final sequence",
            "initial protein",
            "final protein",
            "rejected mutations count",
        ]

        for entry in self.sequence_histories:
            history = entry.get("history", [])
            final_stats = entry.get("final_stats", {})

            if history:
                first_entry = history[0]
                last_entry = history[-1]
                initial_sequence = first_entry.get("sequence_before") or first_entry.get("sequence") or first_entry.get("original_sequence") or ""
                final_sequence = last_entry.get("sequence_after") or last_entry.get("sequence") or last_entry.get("final_sequence") or initial_sequence
            else:
                initial_sequence = ""
                final_sequence = ""

            initial_protein = "".join(dna_sequence.to_amino_acids(initial_sequence)) if initial_sequence else ""
            final_protein = "".join(dna_sequence.to_amino_acids(final_sequence)) if final_sequence else ""
            nucleotide_diff = final_stats.get("nucleotides_diff")
            amino_acid_diff = final_stats.get("amino_acids_diff")
            if nucleotide_diff is None:
                nucleotide_diff = sum(1 for orig, curr in zip(initial_sequence, final_sequence) if orig != curr)
            if amino_acid_diff is None:
                amino_acid_diff = sum(1 for orig, curr in zip(initial_protein, final_protein) if orig != curr)

            aa_to_nt_ratio = round((amino_acid_diff / nucleotide_diff), 4) if nucleotide_diff else 0.0
            rows.append(
                {
                    "sequence ID": entry.get("sequence_id", 0),
                    "mutations number": final_stats.get("mutations_count", 0),
                    "nucleotides difference": nucleotide_diff,
                    "amino acids difference": amino_acid_diff,
                    "amino acids difference to nucleotides difference ratio": aa_to_nt_ratio,
                    "initial sequence": initial_sequence,
                    "final sequence": final_sequence,
                    "initial protein": initial_protein,
                    "final protein": final_protein,
                    "rejected mutations count": final_stats.get("rejected_mutations_count", 0),
                }
            )

        with self.sequence_csv_path.open("w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return self.sequence_csv_path

    def save_summary_csv(self) -> Path:
        rows: list[dict[str, Any]] = []
        fieldnames = [
            "sequence ID",
            "mutations number",
            "nucleotides difference",
            "amino acids difference",
            "amino acids difference to nucleotides difference ratio",
            "initial sequence",
            "final sequence",
            "initial protein",
            "final protein",
            "rejected mutations count",
        ]

        for entry in self.sequence_histories:
            history = entry.get("history", [])
            final_stats = entry.get("final_stats", {})

            if history:
                first_entry = history[0]
                last_entry = history[-1]
                initial_sequence = first_entry.get("sequence_before") or first_entry.get("sequence") or first_entry.get("original_sequence") or ""
                final_sequence = last_entry.get("sequence_after") or last_entry.get("sequence") or last_entry.get("final_sequence") or initial_sequence
            else:
                initial_sequence = ""
                final_sequence = ""

            initial_protein = "".join(dna_sequence.to_amino_acids(initial_sequence)) if initial_sequence else ""
            final_protein = "".join(dna_sequence.to_amino_acids(final_sequence)) if final_sequence else ""
            nucleotide_diff = final_stats.get("nucleotides_diff")
            amino_acid_diff = final_stats.get("amino_acids_diff")
            if nucleotide_diff is None:
                nucleotide_diff = sum(1 for orig, curr in zip(initial_sequence, final_sequence) if orig != curr)
            if amino_acid_diff is None:
                amino_acid_diff = sum(1 for orig, curr in zip(initial_protein, final_protein) if orig != curr)

            aa_to_nt_ratio = round((amino_acid_diff / nucleotide_diff), 4) if nucleotide_diff else 0.0
            rows.append(
                {
                    "sequence ID": entry.get("sequence_id", 0),
                    "mutations number": final_stats.get("mutations_count", 0),
                    "nucleotides difference": nucleotide_diff,
                    "amino acids difference": amino_acid_diff,
                    "amino acids difference to nucleotides difference ratio": aa_to_nt_ratio,

                    "initial sequence": initial_sequence,
                    "final sequence": final_sequence,
                    "initial protein": initial_protein,
                    "final protein": final_protein,
                    "rejected mutations count": final_stats.get("rejected_mutations_count", 0),
                }
            )

        with self.simulation_summary_csv_path.open("w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return self.simulation_summary_csv_path

    def summarize_all_sequences(self) -> None:
        if not self.sequence_histories:
            return

        self.save_summary_csv()

        summary_csv = self.simulation_summary_csv_path if self.simulation_summary_csv_path.exists() else self.sequence_csv_path
        if summary_csv.exists():
            try:
                SimulationDataAnalyzer(summary_csv).analyze()
            except Exception as exc:  # pragma: no cover - defensive fallback for plot generation
                print(f"Warning: failed to generate summary plots from {summary_csv}: {exc}")


    def _write_event(self, payload: dict[str, Any], path: Path | None = None) -> None:
        target = path or self.log_path
        with target.open("a", encoding="utf-8") as file:
            file.write(json.dumps(payload, default=str))
            file.write("\n")


def log_simulation_run(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to log the overall simulation parameters and finalize the summary plots."""

    def wrapper(self, *args: Any, **kwargs: Any) -> Any:
        if not hasattr(self, "logger"):
            self.logger = SimulationLogger(getattr(self, "output_dir", "results"))

        parameters = {
            "mutations_num": self.iterations,
            "sequences_num": self.sequences_num,
            "random_seed": self.random_seed,
            "generate_sequences": self.generate_sequences,
            "input_file": self.input_file,
        }
        self.logger.log_start(parameters)
        result = func(self, *args, **kwargs)
        self.logger.summarize_all_sequences()
        self.output_dir = self.logger.output_dir
        return result

    return wrapper


def log_sequence_run(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to log sequence-level start/end state and capture detailed per-iteration data."""

    def wrapper(self, *args: Any, **kwargs: Any) -> Any:
        if not hasattr(self, "logger"):
            self.logger = SimulationLogger(getattr(self, "output_dir", "results"))

        sequence = args[0] if args else kwargs.get("sequence")
        if sequence is None:
            return func(self, *args, **kwargs)

        self.logger.log_sequence_start(sequence.id, str(sequence), sequence.protein, sequence.length)
        history = func(self, *args, **kwargs)

        final_stats = {
            "nucleotides_diff": sequence.nucleotides_diff,
            "amino_acids_diff": sequence.amino_acids_diff,
            "mutations_count": sequence.mutations_count,
            "rejected_mutations_count": sequence.rejected_mutations_count,
        }
        self.logger.log_sequence_end(sequence.id, history, final_stats)
        return history

    return wrapper
