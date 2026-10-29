"""Download the pinned, unmodified BigMHC IM source and seven-model ensemble."""

import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import httpx

from neovax.io import digest, write_json

COMMIT = "199e890fa69a63646730bcecea0a520e08ca4f02"
REPOSITORY = "https://github.com/KarchinLab/bigmhc"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="data/outputs/bigmhc-pilot/upstream")
    args = parser.parse_args()
    root = Path(args.output_dir).resolve()
    response = httpx.get(f"https://api.github.com/repos/KarchinLab/bigmhc/git/trees/{COMMIT}?recursive=1")
    response.raise_for_status()
    tree = response.json()
    if tree.get("truncated"):
        raise ValueError("Incomplete upstream file tree")
    records = [
        row
        for row in tree["tree"]
        if row["type"] == "blob"
        and (
            row["path"].startswith(("src/", "models/"))
            or row["path"] in {"LICENSE", "README.md", "data/pseudoseqs.csv", "data/example1.csv"}
        )
    ]

    def fetch(record):
        path = root / record["path"]
        url = f"https://raw.githubusercontent.com/KarchinLab/bigmhc/{COMMIT}/{record['path']}"
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            with httpx.stream("GET", url, follow_redirects=True, timeout=180) as download:
                download.raise_for_status()
                temporary = path.with_suffix(f"{path.suffix}.part")
                with temporary.open("wb") as handle:
                    for chunk in download.iter_bytes():
                        handle.write(chunk)
                temporary.replace(path)
        contents = path.read_bytes()
        blob = hashlib.sha1(f"blob {len(contents)}\0".encode() + contents).hexdigest()
        if blob != record["sha"]:
            raise ValueError(f"Official Git blob checksum mismatch: {path}")
        return {"file": record["path"], "url": url, "git_blob": blob, "sha256": digest(path)}

    with ThreadPoolExecutor(max_workers=8) as pool:
        files = list(pool.map(fetch, records))
    write_json(
        root / "manifest.json",
        {
            "repository": REPOSITORY,
            "commit": COMMIT,
            "root": str(root),
            "model": "BigMHC_IM",
            "ensemble_members": 7,
            "aggregation": "mean of seven sigmoid outputs, as in official src/predict.py",
            "license": "BigMHC Academic License; noncommercial use including personal non-profit research",
            "files": files,
        },
    )
    print(f"Verified {len(files)} official files in {root}")


if __name__ == "__main__":
    main()
