from subprocess import CompletedProcess

import pytest

from neovax import runtime
from neovax.config import Config


def test_hla_rescue_is_explicit_and_saved(tmp_path):
    assert not Config().rescue_unknown_hla
    config = Config(rescue_unknown_hla=True)
    path = tmp_path / "config.json"
    path.write_text('{"rescue_unknown_hla": true}')
    assert config.to_dict()["rescue_unknown_hla"]
    assert Config.load(path).rescue_unknown_hla


def existing_model(tmp_path):
    manifest = tmp_path / "outputs/models/deepimmuno/manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}")
    return manifest


def test_reuse_valid_model_and_runtime(tmp_path, monkeypatch):
    manifest = existing_model(tmp_path)
    monkeypatch.setattr(runtime, "validate_model", lambda *args: {"python": "/model/python"})
    monkeypatch.setattr(runtime, "_runtime_ready", lambda python: python == "/model/python")
    messages = []
    assert runtime.ensure_deepimmuno(tmp_path, messages.append) == str(manifest)
    assert messages == ["DeepImmuno model ready"]


def test_cached_only_does_not_require_runtime_or_install(tmp_path, monkeypatch):
    manifest = existing_model(tmp_path)
    monkeypatch.setattr(runtime, "validate_model", lambda *args: {"python": "/missing/python"})
    monkeypatch.setattr(runtime.subprocess, "run", lambda *args, **kwargs: pytest.fail("Started a process"))
    monkeypatch.setattr(runtime, "setup_model", lambda *args: pytest.fail("Downloaded model"))
    assert runtime.ensure_deepimmuno(tmp_path, cached_only=True) == str(manifest)


def test_cached_only_missing_model_leaves_directory_untouched(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "setup_model", lambda *args: pytest.fail("Downloaded model"))
    with pytest.raises(FileNotFoundError, match="without --cached-only"):
        runtime.ensure_deepimmuno(tmp_path, cached_only=True)
    assert list(tmp_path.iterdir()) == []


def test_prepare_isolated_runtime_and_pinned_model(tmp_path, monkeypatch):
    interpreter = tmp_path / "python3.11"
    interpreter.touch()
    monkeypatch.setattr(runtime.shutil, "which", lambda name: str(interpreter))
    readiness = iter([False, True])
    monkeypatch.setattr(runtime, "_runtime_ready", lambda python: next(readiness))
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        return CompletedProcess(command, 0, stdout="3.11\n")

    monkeypatch.setattr(runtime.subprocess, "run", run)
    setups = []

    def setup(directory, python):
        setups.append((directory, python))
        return directory / "manifest.json"

    monkeypatch.setattr(runtime, "setup_model", setup)
    result = runtime.ensure_deepimmuno(tmp_path)
    environment = tmp_path / "outputs/deepimmuno-env"
    model = tmp_path / "outputs/models/deepimmuno"
    assert result == str(model / "manifest.json")
    assert commands[1] == [str(interpreter), "-m", "venv", str(environment)]
    assert commands[2][:4] == [str(environment / "bin/python"), "-m", "pip", "install"]
    assert "tensorflow==2.15.1" in commands[2]
    assert setups == [(model, environment / "bin/python")]


def test_missing_python_has_actionable_error(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime.shutil, "which", lambda name: None)
    monkeypatch.setattr(runtime.sys, "version_info", (3, 12))
    with pytest.raises(RuntimeError, match="Install Python 3.11"):
        runtime.ensure_deepimmuno(tmp_path)


def test_runtime_rejects_missing_or_incompatible_interpreter(tmp_path, monkeypatch):
    assert not runtime._runtime_ready(tmp_path / "absent")
    python = tmp_path / "python"
    python.touch()
    monkeypatch.setattr(runtime.subprocess, "run", lambda *args, **kwargs: CompletedProcess([], 1))
    assert not runtime._runtime_ready(python)
