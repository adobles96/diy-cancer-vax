"""Set up a repository-local environment using standard Python and pip."""

import hashlib
import os
import subprocess
import sys
import venv
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    environment = root / ".venv"
    python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    stamp = environment / ".neovax-bootstrap"
    requirements = hashlib.sha256((root / "pyproject.toml").read_bytes()).hexdigest()
    if not python.exists():
        print("Creating a local Python environment in .venv/ ...", flush=True)
        venv.EnvBuilder(with_pip=True, symlinks=os.name != "nt").create(environment)
    if not stamp.exists() or stamp.read_text() != requirements:
        print("Installing DIY Cancer Vaccine and its dependencies (first run or changed requirements) ...", flush=True)
        pip_check = subprocess.run(
            [str(python), "-m", "pip", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        if pip_check.returncode:
            subprocess.run([str(python), "-m", "ensurepip"], check=True)
        subprocess.run(
            [str(python), "-m", "pip", "install", "--disable-pip-version-check", "--quiet", "-e", str(root)],
            check=True,
        )
        stamp.write_text(requirements)
    arguments = sys.argv[1:]
    if arguments and arguments[0].startswith("-") and arguments[0] not in {"-h", "--help"}:
        arguments = ["run", *arguments]
    os.chdir(root)
    os.execv(str(python), [str(python), "-m", "neovax.cli", *arguments])


if __name__ == "__main__":
    main()
