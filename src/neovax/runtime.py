"""Prepare the optional model without changing the main application environment."""

import shutil
import subprocess
import sys
from pathlib import Path

from neovax.deepimmuno import setup_model, validate_model


def _runtime_ready(python):
    if not Path(python).is_file():
        return False
    result = subprocess.run(
        [
            str(python),
            "-c",
            "import sys, tensorflow, numpy, pandas; "
            "assert sys.version_info[:2] == (3, 11); "
            "assert tensorflow.__version__ == '2.15.1'",
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def ensure_deepimmuno(data_root: Path, progress=lambda message: None, cached_only=False) -> str:
    """Return a pinned model manifest, installing its isolated runtime when needed.

    ``data_root`` is the directory containing ``raw`` and ``outputs``. Cached-only
    runs need the validated model files but do not need an inference interpreter.
    """
    outputs = Path(data_root).resolve() / "outputs"
    directory = outputs / "models" / "deepimmuno"
    manifest_path = directory / "manifest.json"
    manifest = None
    if manifest_path.is_file():
        try:
            manifest = validate_model(manifest_path, [], (9, 10))
        except (OSError, ValueError, KeyError):
            pass
    if manifest and (cached_only or _runtime_ready(manifest["python"])):
        progress("DeepImmuno model ready")
        return str(manifest_path)
    if cached_only:
        raise FileNotFoundError(
            "The pinned DeepImmuno model is unavailable. Run once without --cached-only to prepare it."
        )

    environment = outputs / "deepimmuno-env"
    python = environment / "bin" / "python"
    if not _runtime_ready(python):
        candidates = [shutil.which("python3.11")]
        if sys.version_info[:2] == (3, 11):
            candidates.append(sys.executable)
        candidates.extend([manifest.get("python") if manifest else None, str(python)])
        interpreter = None
        for candidate in candidates:
            if not candidate or not Path(candidate).is_file():
                continue
            candidate = Path(candidate).resolve()
            version = subprocess.run(
                [
                    str(candidate),
                    "-c",
                    "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')",
                ],
                capture_output=True,
                text=True,
            )
            if version.returncode == 0 and version.stdout.strip() == "3.11":
                interpreter = candidate
                break
        if interpreter is None:
            raise RuntimeError(
                "DeepImmuno needs Python 3.11. Install Python 3.11, make python3.11 available on PATH, "
                "then rerun; its model environment will be prepared automatically."
            )
        progress("Preparing isolated DeepImmuno environment (Python 3.11 / TensorFlow 2.15.1)")
        subprocess.run([str(interpreter), "-m", "venv", str(environment)], check=True)
        subprocess.run(
            [
                str(python),
                "-m",
                "pip",
                "install",
                "--disable-pip-version-check",
                "--quiet",
                "tensorflow==2.15.1",
                "numpy==1.26.4",
                "pandas>=2,<4",
            ],
            check=True,
        )
        if not _runtime_ready(python):
            raise RuntimeError(f"DeepImmuno runtime verification failed: {python}")
    progress("Downloading missing DeepImmuno model files and verifying pinned hashes")
    return str(setup_model(directory, python))
