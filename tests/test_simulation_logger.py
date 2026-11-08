from pathlib import Path

from modules.simulation_logger import SimulationLogger


def test_save_sequence_csv_exports_required_columns(tmp_path):
    logger = SimulationLogger(tmp_path)
    logger.sequence_histories = [
        {
            "sequence_id": 7,
            "history": [
                {"iteration": 0, "sequence_before": "ATGAAA", "sequence_after": "ATGCAA"},
                {"iteration": 1, "sequence_before": "ATGCAA", "sequence_after": "ATGCAA"},
            ],
            "final_stats": {
                "nucleotides_diff": 1,
                "amino_acids_diff": 1,
                "mutations_count": 2,
                "rejected_mutations_count": 2,
            },
        }
    ]

    csv_path = logger.save_sequence_csv()

    assert csv_path.exists()
    assert csv_path.name == "sequence_summary.csv"
    text = csv_path.read_text(encoding="utf-8")
    assert "sequence ID" in text
    assert "mutations number" in text
    assert "nucleotides difference" in text
    assert "amino acids difference" in text
    assert "amino acids difference to nucleotides difference ratio" in text
    assert "initial sequence" in text
    assert "final sequence" in text
    assert "initial protein" in text
    assert "final protein" in text
    assert "rejected mutations count" in text
    assert "A change ratio" in text
    assert "T change ratio" in text
    assert "C change ratio" in text
    assert "G change ratio" in text
    assert "ATGAAA" in text
    assert "ATGCAA" in text


def test_save_sequence_history_csv_exports_step_rows(tmp_path):
    logger = SimulationLogger(tmp_path)
    history = [
        {"iteration": -1, "mutation": None, "sequence": "ATGAAA", "stats": {"nucleotides_diff": 0, "amino_acids_diff": 0, "mutations_count": 0, "rejected_mutations_count": 0}},
        {"iteration": 0, "mutation": "2;A;C", "sequence": "ATGCAA", "stats": {"nucleotides_diff": 1, "amino_acids_diff": 1, "mutations_count": 1, "rejected_mutations_count": 0}},
        {"iteration": 1, "mutation": "5;G;T", "sequence": "ATGCTA", "stats": {"nucleotides_diff": 2, "amino_acids_diff": 1, "mutations_count": 2, "rejected_mutations_count": 1}},
    ]

    csv_path = logger.save_sequence_history_csv(9, history)

    assert csv_path.exists()
    assert csv_path.name == "sequence_9.csv"
    assert csv_path.parent == logger.sequences_dir
    assert (logger.sequences_dir / "sequence_9.png").exists()
    text = csv_path.read_text(encoding="utf-8")
    assert "ID" in text
    assert "mutation number" in text
    assert "mutation" in text
    assert "nucleotides diff" in text
    assert "amino acids diff" in text
    assert "rejected mutations count" in text
    assert "sequence" in text
    assert "amino acids sequence" in text
    lines = text.strip().splitlines()
    assert lines[1].startswith("9,0,")
    assert lines[-1].startswith("9,2,")
    assert "ATGAAA" in text
    assert "ATGCTA" in text
