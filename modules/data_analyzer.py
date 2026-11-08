import csv
import sys
from pathlib import Path
from typing import Any
from re import fullmatch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class SimulationDataAnalyzer:
    """Read a simulation summary CSV and produce summary plots."""

    @staticmethod
    def _add_statistic_lines(
        ax: Any,
        values: list[float],
        label_prefix: str,
        color: str,
    ) -> None:
        if not values:
            return

        statistics = (
            ("Max", max(values), "--"),
            ("Min", min(values), "-."),
            ("Avg", sum(values) / len(values), ":"),
        )
        for label, value, linestyle in statistics:
            ax.axhline(
                value,
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
                label="_nolegend_",
            )
            ax.annotate(
                f"{label_prefix} {label}: {value:g}",
                xy=(1, value),
                xycoords=("axes fraction", "data"),
                xytext=(-4, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                color=color,
                fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none", "pad": 1},
            )

    @staticmethod
    def _normalize_path(raw_path: str | Path) -> Path:
        text = str(raw_path).strip()
        if not text:
            raise ValueError("CSV path is empty.")

        while len(text) >= 2 and text[0] == text[-1] and text[0] in {'"', "'"}:
            text = text[1:-1].strip()

        if text.startswith("&"):
            text = text[1:].strip()

        shell_alias = " & "
        if shell_alias in text:
            left, right = text.split(shell_alias, 1)
            if left.strip():
                text = left.strip()
            else:
                text = right.strip()

        while len(text) >= 2 and text[0] == text[-1] and text[0] in {'"', "'"}:
            text = text[1:-1].strip()

        path = Path(text).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        return path

    def __init__(self, csv_path: str | Path) -> None:
        self.csv_path = self._normalize_path(csv_path)

    @staticmethod
    def plot_sequence_history(
        seq_id: int,
        history: list[dict[str, Any]] | None = None,
        csv_path: str | Path | None = None,
        output_dir: str | Path | None = None,
    ) -> Path:
        if csv_path is not None:
            csv_file = Path(csv_path)
            if not csv_file.exists():
                raise FileNotFoundError(f"Sequence CSV not found: {csv_file}")
            with csv_file.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                rows = list(reader)
            if not rows:
                raise ValueError(f"Sequence CSV is empty: {csv_file}")
            seq_id = seq_id or int(float(rows[0].get("ID", 0) or 0))
            sequence = rows[1].get("sequence", "???")
            iterations = [int(float(row.get("mutation number", row.get("mutation_count", 0)) or 0)) for row in rows]
            nucleotide_diffs = [float(row.get("nucleotides diff", row.get("nucleotides_difference", 0)) or 0) for row in rows]
            amino_acid_diffs = [float(row.get("amino acids diff", row.get("amino_acids_difference", 0)) or 0) for row in rows]
            target_dir = Path(output_dir) if output_dir is not None else csv_file.parent
        else:
            if not history:
                raise ValueError("No history data was provided for the sequence plot.")
            sequence = history[0].get("sequence", "???")
            iterations = [entry.get("iteration", index) for index, entry in enumerate(history)]
            nucleotide_diffs = [float(entry.get("stats", {}).get("nucleotides_diff", 0) or 0) for entry in history]
            amino_acid_diffs = [float(entry.get("stats", {}).get("amino_acids_diff", 0) or 0) for entry in history]
            target_dir = Path(output_dir) if output_dir is not None else Path.cwd()

        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(iterations, nucleotide_diffs, label="Nucleotide diff", marker="o", linewidth=2, color="tab:blue")
        ax.plot(iterations, amino_acid_diffs, label="Amino acid diff", marker="s", linewidth=2, color="tab:orange")
        SimulationDataAnalyzer._add_statistic_lines(ax, nucleotide_diffs, "Nucleotide", "tab:blue")
        SimulationDataAnalyzer._add_statistic_lines(ax, amino_acid_diffs, "Amino acid", "tab:orange")
        ax.set_title(f"Sequence {seq_id} mutation trajectory (length={len(sequence)})")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Difference count")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        output_path = target_dir / f"sequence_{seq_id}.png"
        fig.savefig(output_path, dpi=200)
        plt.close(fig)
        return output_path

    def _read_csv_rows(self) -> list[dict[str, str]]:
        if not self.csv_path.exists():
            raise FileNotFoundError(f"Simulation summary CSV not found: {self.csv_path}")

        with self.csv_path.open("r", newline="", encoding="utf-8") as csv_file:
            reader = csv.DictReader(csv_file)
            rows = list(reader)

        if not rows:
            raise ValueError(f"Simulation summary CSV is empty: {self.csv_path}")

        return rows

    def _value(self, row: dict[str, str], *names: str) -> float:
        for name in names:
            if name in row and row.get(name, "") not in (None, ""):
                value = row.get(name, "").strip()
                try:
                    return float(value)
                except ValueError:
                    return 0.0
        return 0.0

    def _sequence_id(self, row: dict[str, str]) -> int:
        for name in ("sequence ID", "sequence_id", "ID"):
            if name in row and row.get(name, "").strip() not in ("", "None"):
                try:
                    return int(float(row[name]))
                except ValueError:
                    return 0
        return 0

    def _build_arrays(self) -> tuple[list[int], list[float], list[float], list[float]]:
        rows = self._read_csv_rows()
        sequence_ids: list[int] = []
        nucleotide_diffs: list[int] = []
        amino_acid_diffs: list[int] = []
        ratios: list[float] = []

        for row in rows:
            sequence_ids.append(self._sequence_id(row))
            nucleotide_diffs.append(self._value(row, "nucleotides difference", "nucleotide_difference_count", "nucleotides_diff"))
            amino_acid_diffs.append(self._value(row, "amino acids difference", "amino_acids_difference_count", "amino_acids_diff"))
            ratios.append(self._value(row, "amino acids difference to nucleotides difference ratio", "differences_ratio", "aa_to_nt_ratio"))

        return sequence_ids, nucleotide_diffs, amino_acid_diffs, ratios

    def analyze(self) -> Path:
        sequence_ids, nucleotide_diffs, amino_acid_diffs, ratios = self._build_arrays()
        output_path = self.csv_path.with_name("simulation_summary_analysis.png")

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"Simulation summary analysis: {self.csv_path.name}, length={len(sequence_ids)}", fontsize=14)

        if sequence_ids:
            axes[0].bar(sequence_ids, nucleotide_diffs, label="Nucleotide diff", color="tab:blue")
            axes[0].bar(sequence_ids, amino_acid_diffs, label="Amino acid diff", color="tab:orange")
            self._add_statistic_lines(axes[0], nucleotide_diffs, "Nucleotide", "tab:blue")
            self._add_statistic_lines(axes[0], amino_acid_diffs, "Amino acid", "tab:orange")
            axes[0].set_title("Per-sequence difference trends")
            axes[0].set_xlabel("Sequence ID")
            axes[0].set_ylabel("Difference count")
            axes[0].legend()
            axes[0].grid(alpha=0.3)
        else:
            axes[0].text(0.5, 0.5, "No sequence data", ha="center", va="center", transform=axes[0].transAxes)
            axes[0].set_title("Per-sequence difference trends")

        if ratios:
            axes[1].bar(sequence_ids, ratios, align="center", color="tab:green", edgecolor="black")
            self._add_statistic_lines(axes[1], ratios, "Ratio", "black")
            axes[1].set_title("Amino acid / nucleotide ratio distribution")
            axes[1].set_xlabel("Sequence ID")
            axes[1].set_ylabel("Ratio")
            axes[1].text(0.99, 0.01, f"Total sequences: {len(ratios)}", transform=axes[1].transAxes, ha="right", va="bottom", fontsize=9)
        else:
            axes[1].text(0.5, 0.5, "No ratio data", ha="center", va="center", transform=axes[1].transAxes)
            axes[1].set_title("Amino acid / nucleotide ratio distribution")

        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(output_path, dpi=200)
        plt.close(fig)
        return output_path


def _prompt_for_csv_path() -> Path:
    raw_path = input("Enter the path to the simulation_summary.csv or sequence CSV file: ").strip().strip('"')
    if not raw_path:
        raise ValueError("No CSV path was entered.")
    return Path(raw_path).expanduser()


def _is_sequence_csv(csv_path: Path) -> bool:
    return bool(fullmatch(r"sequence_\d+\.csv", csv_path.name, flags=0))

def main() -> None:
    raw_path = sys.argv[1] if len(sys.argv) > 1 else _prompt_for_csv_path()
    csv_path = SimulationDataAnalyzer._normalize_path(raw_path)
    analyzer = SimulationDataAnalyzer(csv_path)
    if _is_sequence_csv(csv_path):
        output_path = analyzer.plot_sequence_history(seq_id=0, csv_path=csv_path)
        print(f"Sequence trajectory plot saved to: {output_path}")
    else:
        output_path = analyzer.analyze()
        print(f"Analysis complete. Plots saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)
