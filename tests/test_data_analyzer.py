from pathlib import Path

from modules.data_analyzer import SimulationDataAnalyzer


def test_data_analyzer_creates_summary_plots_from_csv(tmp_path):
    csv_path = tmp_path / "simulation_summary.csv"
    csv_path.write_text(
        "sequence ID,mutations number,nucleotides difference,amino acids difference,amino acids difference to nucleotides difference ratio,initial sequence,final sequence,initial protein,final protein,rejected mutations count\n"
        "0,3,5,2,0.4,ATGAAA,ATGCAA,MXK,MXK,1\n"
        "1,2,2,1,0.5,GGGCCC,GGGTTT,AAA,BBB,0\n",
        encoding="utf-8",
    )

    analyzer = SimulationDataAnalyzer(csv_path)
    plot_path = analyzer.analyze()

    assert plot_path.exists()
    assert plot_path.name == "simulation_summary_analysis.png"
    assert plot_path.parent == csv_path.parent
