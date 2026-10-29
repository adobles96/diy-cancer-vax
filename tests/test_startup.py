import hashlib
import io
import json
from dataclasses import replace

import httpx
import pytest
from rich.console import Console

from neovax import cli, data
from neovax.config import Config


@pytest.fixture
def sources(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    payloads = {"table.xlsx": b"test workbook", "reference.txt.gz": b"test proteins"}
    manifest = [
        {"file": name, "url": f"https://sources.test/{name}", "sha256": hashlib.sha256(body).hexdigest()}
        for name, body in payloads.items()
    ]
    monkeypatch.setattr(cli, "source_manifest", lambda: manifest)
    monkeypatch.setattr(data, "source_manifest", lambda: manifest)
    data.verify_sources.cache_clear()
    requested = []

    def respond(request):
        requested.append(request.url.path)
        return httpx.Response(200, content=payloads[request.url.path[1:]])

    with httpx.Client(transport=httpx.MockTransport(respond)) as client:
        monkeypatch.setattr(data.httpx, "stream", lambda method, url, **kwargs: client.stream(method, url))
        yield requested, manifest
    data.verify_sources.cache_clear()


def test_first_run_downloads_to_chosen_root_and_remembers_it(tmp_path, monkeypatch, sources):
    requested, manifest = sources
    root = tmp_path / "chosen data"
    console = Console(file=io.StringIO())
    prompts = []
    monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(console, "input", lambda prompt, **kwargs: prompts.append(prompt) or str(root))
    config = cli.prepare_inputs(Config(patient=3), None, False, console)
    assert prompts and len(requested) == 2
    assert config.patient == 3
    assert config.data_dir == str(root / "raw")
    assert config.cache_dir == str(root / "outputs/cache")
    assert config.output_dir == str(root / "outputs/runs")
    assert all((root / "raw" / source["file"]).exists() for source in manifest)
    assert json.loads((tmp_path / ".neovax.json").read_text())["data_root"] == str(root)

    monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: False)
    again = cli.prepare_inputs(Config(patient=3), None, False, console)
    assert again == config and len(requested) == 2 and len(prompts) == 1


def test_partial_download_reuses_valid_files_and_repairs_missing_ones(tmp_path, sources):
    requested, _ = sources
    root = tmp_path / "data"
    (root / "raw").mkdir(parents=True)
    workbook = root / "raw/table.xlsx"
    workbook.write_bytes(b"test workbook")
    original_mtime = workbook.stat().st_mtime_ns
    cli.prepare_inputs(Config(), root, False, Console(file=io.StringIO()))
    assert requested == ["/reference.txt.gz"]
    assert workbook.stat().st_mtime_ns == original_mtime


def test_noninteractive_first_run_requires_explicit_destination(monkeypatch, sources):
    requested, _ = sources
    monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: False)
    with pytest.raises(ValueError, match="--data-root"):
        cli.prepare_inputs(Config(), None, False, Console(file=io.StringIO()))
    assert requested == []


def test_cached_only_never_downloads_sources(tmp_path, sources):
    requested, _ = sources
    with pytest.raises(FileNotFoundError):
        cli.prepare_inputs(Config(cached_only=True), tmp_path / "offline", False, Console(file=io.StringIO()))
    assert requested == [] and not (tmp_path / ".neovax.json").exists()


def test_explicit_config_paths_take_priority_over_saved_root(tmp_path, sources):
    (tmp_path / ".neovax.json").write_text(json.dumps({"data_root": str(tmp_path / "saved")}))
    config = replace(
        Config(), data_dir=str(tmp_path / "custom/raw"), cache_dir="custom-cache", output_dir="custom-runs"
    )
    prepared = cli.prepare_inputs(config, None, True, Console(file=io.StringIO()))
    assert prepared == config
    assert (tmp_path / "custom/raw/table.xlsx").exists()
    assert not (tmp_path / "saved").exists()


def test_failed_checksum_does_not_remember_destination(tmp_path, monkeypatch, sources):
    _, manifest = sources
    bad_manifest = [dict(manifest[0], sha256="0" * 64)]
    monkeypatch.setattr(data, "source_manifest", lambda: bad_manifest)
    with pytest.raises(ValueError, match="Hash mismatch"):
        cli.prepare_inputs(Config(), tmp_path / "bad", False, Console(file=io.StringIO()))
    assert not (tmp_path / ".neovax.json").exists()
    assert not (tmp_path / "bad/raw/table.xlsx").exists()
