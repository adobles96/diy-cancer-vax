import csv
import io
import json
from pathlib import Path

import pytest
from rich.console import Console

from neovax import cli
from neovax.config import Config
from neovax.iedb import Cancelled
from neovax.io import digest
from neovax.report import render_report


class TerminalInput(io.StringIO):
    def isatty(self):
        return True


def menu(monkeypatch, answers):
    output = io.StringIO()
    console = Console(file=output, width=100)
    monkeypatch.setattr(cli, "Console", lambda: console)
    monkeypatch.setattr(cli.sys, "stdin", TerminalInput(answers))
    monkeypatch.setattr(cli.sys, "argv", ["neovax"])
    monkeypatch.setattr(cli, "prepare_inputs", lambda config, *args: config)
    jobs = []
    monkeypatch.setattr(cli, "execute_job", lambda *args: jobs.append(args))
    return output, jobs


def test_menu_patient_defaults(monkeypatch):
    output, jobs = menu(monkeypatch, "1\n3\ny\n")
    cli.main()
    config, command, _, _ = jobs[0]
    assert config == Config(patient=3)
    assert command == "run"
    assert "Use defaults?" in output.getvalue()
    assert "DIY Cancer Vaccine" in output.getvalue()


def test_menu_custom_patient_and_model_setup(monkeypatch):
    # Invalid budget is retried before accepting five targets.
    output, jobs = menu(monkeypatch, "1\n2\nn\n3\n9\n0\n5\n21\n9\n501\n125.5\ny\n")
    setup_calls = []
    monkeypatch.setattr(
        cli, "ensure_deepimmuno", lambda root, **kwargs: setup_calls.append(root) or "model.json"
    )
    cli.main()
    config = jobs[0][0]
    assert config.patient == 2 and config.ranking == "combined"
    assert config.lengths == (9,) and config.budget == 5 and config.long_length == 21
    assert config.threshold == 125.5 and config.rescue_unknown_hla
    assert "Enter a number from 10 to 500" in output.getvalue()
    assert config.deepimmuno_manifest == "model.json" and setup_calls == [Path("data")]
    assert "Enter a number from 1 to 100" in output.getvalue()


def test_blog_menu_announces_destination(monkeypatch):
    output, jobs = menu(monkeypatch, "2\n")
    monkeypatch.setattr(cli, "ensure_deepimmuno", lambda *args, **kwargs: "model.json")
    cli.main()
    config, command, destination, _ = jobs[0]
    assert command == "blog" and destination == Path("data/outputs/blog")
    assert config.deepimmuno_manifest == "model.json"
    assert "Blog results will be saved in" in output.getvalue()
    assert "six patients" in output.getvalue()


def test_missing_python_shows_actionable_error_without_traceback(monkeypatch):
    output, _ = menu(monkeypatch, "2\n")

    def missing(*args, **kwargs):
        raise RuntimeError("Install Python 3.11")

    monkeypatch.setattr(cli, "ensure_deepimmuno", missing)
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2
    assert "Could not finish: Install Python 3.11" in output.getvalue()


def test_cancelled_blog_has_cancel_exit(monkeypatch):
    output, _ = menu(monkeypatch, "2\n")
    monkeypatch.setattr(cli, "ensure_deepimmuno", lambda *args, **kwargs: "model.json")

    def cancel(*args):
        raise Cancelled("Cancelled")

    monkeypatch.setattr(cli, "execute_job", cancel)
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 130
    assert "Stopped" in output.getvalue()


def test_menu_requires_terminal(monkeypatch):
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO())
    monkeypatch.setattr(cli.sys, "argv", ["neovax"])
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2


@pytest.mark.integration
@pytest.mark.parametrize("patient", range(1, 7))
def test_real_patient_cli_exports_and_prints_shortlist(patient, tmp_path, monkeypatch):
    if not Path("data/outputs/cache/binding.sqlite").exists():
        pytest.skip("Patient prediction cache required")
    output = io.StringIO()
    console = Console(file=output, width=100)
    monkeypatch.setattr(cli, "Console", lambda: console)
    config = Config(patient=patient, cached_only=True, output_dir=str(tmp_path / "runs"))
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config.to_dict()))
    monkeypatch.setattr(cli.sys, "argv", ["neovax", "run", "--config", str(config_path)])
    cli.main()
    text = output.getvalue()
    assert f"Patient {patient} · COMPLETE" in text
    assert "Gene /" in text and "Driving" in text and "IC50 (nM)" in text
    runs = list((tmp_path / "runs").iterdir())
    assert len(runs) == 1
    data = json.loads((runs[0] / "run.json").read_text())
    assert len(data["proposed"]) == 20
    assert all(row["protein_change"] for row in data["proposed"])
    assert {path.name for path in runs[0].iterdir()} == {"report.html", "proposed_vaccine.csv", "run.json"}
    with (runs[0] / "proposed_vaccine.csv").open() as handle:
        shortlist = list(csv.DictReader(handle))
    assert len(shortlist) == 20
    assert [row["sequence"] for row in shortlist] == [row["sequence"] for row in data["proposed"]]
    assert [float(row["ic50_nm"]) for row in shortlist] == [row["ic50_nm"] for row in data["proposed"]]
    assert {"config", "manifest", "summary", "pairs", "assays", "reference_filtered"} <= data.keys()
    assert data["manifest"]["exports"] == {
        name: digest(runs[0] / name) for name in ["report.html", "proposed_vaccine.csv"]
    }
    if patient == 1:
        render_report(runs[0])
        saved = json.loads((runs[0] / "run.json").read_text())
        assert saved["manifest"]["exports"]["report.html"] == digest(runs[0] / "report.html")
        assert {path.name for path in runs[0].iterdir()} == {
            "report.html",
            "proposed_vaccine.csv",
            "run.json",
        }


def test_shortlist_never_truncates_sequences_in_narrow_terminal():
    output = io.StringIO()
    sequence = "ACDEFGHIKLMNPQRSTVWYACDEFG"
    cli.print_summary(
        {
            "config": {"ranking": "binding"},
            "summary": {
                "patient": 1,
                "status": "complete",
                "errors": [],
                "run_directory": "/tmp/results",
                "overlap": {
                    "target_overlap": 1,
                    "historical_components": 13,
                    "reconstructable_historical_mutations": 11,
                },
            },
            "proposed": [
                {
                    "gene": "GENE1",
                    "protein_change": "p.A123V",
                    "sequence": sequence,
                    "driving_epitope": "LMNPQRSTV",
                    "ranking_cost": 12.4,
                    "historical_target": True,
                }
            ],
        },
        Console(file=output, width=70),
    )
    text = output.getvalue()
    assert "…" not in text
    # Recover the long-peptide column across wrapped lines; no residues may disappear.
    assert all(residue in text for residue in set(sequence))
    assert "p.A123V" in text


@pytest.mark.parametrize("lengths, expected", [("9", (9,)), ("10", (10,)), ("both", (9, 10))])
def test_clear_lengths_and_numeric_cutoff(monkeypatch, lengths, expected):
    output, jobs = menu(monkeypatch, f"1\n1\nn\n1\n{lengths}\n20\n25\n50\n")
    cli.main()
    assert jobs[0][0].lengths == expected
    assert jobs[0][0].threshold == 50
    assert "[9/10/both]" in output.getvalue()
    assert "[500/50/none]" not in output.getvalue()


def test_custom_cutoff_default_is_visible_and_accepts_enter(monkeypatch):
    output, jobs = menu(monkeypatch, "1\n1\nn\n1\nboth\n20\n25\n\n")
    cli.main()
    assert jobs[0][0].threshold == 500
    assert "pMHC binding affinity cutoff (IC50, nM; 10–500) (500.0)" in output.getvalue()
