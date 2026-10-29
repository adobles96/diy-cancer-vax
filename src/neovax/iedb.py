import csv
import datetime
import hashlib
import io
import json
import math
import sqlite3
import time
from email.utils import parsedate_to_datetime
from pathlib import Path

import httpx

from neovax import __version__
from neovax.io import digest, write_json
from neovax.peptides import AMINO_ACIDS

ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
METHOD = "netmhcpan_ba-4.1"


class InferenceError(ValueError):
    pass


class Cancelled(InferenceError):
    pass


def validate_response(text, sequences, allele, length):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    expected = {
        (index, start, start + length - 1)
        for index, sequence in enumerate(sequences, 1)
        for start in range(1, len(sequence) - length + 2)
    }
    seen, validated, peptide_scores = set(), [], {}
    required = {"seq_num", "start", "end", "length", "allele", "peptide", "ic50", "percentile_rank"}
    for row in rows:
        if not required <= row.keys():
            raise InferenceError("IEDB returned an unexpected schema")
        index, start, end = int(row["seq_num"]), int(row["start"]), int(row["end"])
        key = (index, start, end)
        if key not in expected or key in seen or row["allele"] != allele or int(row["length"]) != length:
            raise InferenceError(
                "IEDB returned inconsistent sequence coordinates, allele, length or duplicate rows"
            )
        if sequences[index - 1][start - 1 : end] != row["peptide"]:
            raise InferenceError("IEDB sequence mapping mismatch")
        ic50, percentile = float(row["ic50"]), float(row["percentile_rank"])
        if (
            not math.isfinite(ic50)
            or ic50 <= 0
            or not math.isfinite(percentile)
            or not 0 <= percentile <= 100
        ):
            raise InferenceError("IEDB returned non-finite or invalid scores")
        previous = peptide_scores.get(row["peptide"])
        if previous and (
            not math.isclose(previous[0], ic50, abs_tol=0.02)
            or not math.isclose(previous[1], percentile, abs_tol=0.02)
        ):
            raise InferenceError("IEDB returned inconsistent scores for an identical peptide–HLA pair")
        peptide_scores[row["peptide"]] = (ic50, percentile)
        seen.add(key)
        validated.append(
            {
                "peptide": row["peptide"],
                "allele": allele,
                "ic50": ic50,
                "percentile_rank": percentile,
                "seq_num": index,
                "start": start,
                "end": end,
            }
        )
    if seen != expected:
        raise InferenceError(
            f"Incomplete IEDB response: {len(seen)} of {len(expected)} expected windows; response begins {text[:160]!r}"
        )
    return validated


class BindingCache:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.directory / "binding.sqlite")
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS predictions (endpoint TEXT, method TEXT, allele TEXT, peptide TEXT, ic50 REAL, percentile REAL, request_id TEXT, PRIMARY KEY(endpoint,method,allele,peptide))"
        )

    def close(self):
        self.db.close()

    def lookup(self, alleles):
        result = {}
        for allele in alleles:
            for peptide, ic50, percentile, request_id in self.db.execute(
                "SELECT peptide,ic50,percentile,request_id FROM predictions WHERE endpoint=? AND method=? AND allele=?",
                (ENDPOINT, METHOD, allele),
            ):
                result[(peptide, allele)] = {
                    "peptide": peptide,
                    "allele": allele,
                    "ic50": ic50,
                    "percentile_rank": percentile,
                    "request_id": request_id,
                }
        return result

    def save(self, rows, request_id):
        with self.db:
            self.db.executemany(
                "INSERT OR REPLACE INTO predictions VALUES (?,?,?,?,?,?,?)",
                [
                    (ENDPOINT, METHOD, r["allele"], r["peptide"], r["ic50"], r["percentile_rank"], request_id)
                    for r in rows
                ],
            )

    def import_legacy(self, directory):
        directory = Path(directory)
        fasta = (directory / "prediction_input.fasta").read_text()
        sequences = parse_fasta(fasta)
        count = 0
        for path in sorted(directory.glob("iedb_*.json")):
            metadata = json.loads(path.read_text())
            response = path.with_suffix(".tsv")
            if metadata["method"] != METHOD or metadata["endpoint"] != ENDPOINT:
                raise InferenceError("Legacy cache method/endpoint mismatch")
            if (
                hashlib.sha256(fasta.encode()).hexdigest() != metadata["input_sha256"]
                or digest(response) != metadata["response_sha256"]
            ):
                raise InferenceError("Legacy cache hash mismatch")
            rows = validate_response(response.read_text(), sequences, metadata["allele"], metadata["length"])
            request_id = f"legacy-{metadata['response_sha256']}"
            target = self.directory / "requests" / request_id
            target.mkdir(parents=True, exist_ok=True)
            write_json(
                target / "metadata.json",
                {**metadata, "software_version": __version__, "legacy_import": True, "sequence_text": fasta},
            )
            (target / "response.tsv").write_text(response.read_text())
            self.save(rows, request_id)
            count += len(rows)
        return count


def parse_fasta(text):
    sequences = []
    for line in text.strip().splitlines():
        if line.startswith(">"):
            sequences.append("")
        else:
            sequences[-1] += line.strip()
    return sequences


def contexts_for(proteins, lengths):
    sequences = set()
    flank = max(lengths) - 1
    for protein in proteins:
        p = protein["position"]
        for state in ("mutant_protein", "reference_protein"):
            sequence = protein[state][max(0, p - flank) : p + flank + 1]
            if not set(sequence) - AMINO_ACIDS:
                sequences.add(sequence)
            else:
                # Preserve valid windows when an ambiguous residue occurs elsewhere in context.
                for length in lengths:
                    for start in range(max(0, p - length + 1), min(p, len(protein[state]) - length) + 1):
                        window = protein[state][start : start + length]
                        if not set(window) - AMINO_ACIDS:
                            sequences.add(window)
    return sorted(sequences)


def plan_batches(sequences, alleles, lengths, lookup, batch_size):
    batches = []
    for allele in alleles:
        for length in lengths:
            missing = [
                s
                for s in sequences
                if len(s) >= length
                and any((s[i : i + length], allele) not in lookup for i in range(len(s) - length + 1))
            ]
            batches.extend(
                (allele, length, missing[i : i + batch_size]) for i in range(0, len(missing), batch_size)
            )
    return batches


def request_batch(client, payload, cancel, attempts=3, sleep=time.sleep):
    for attempt in range(attempts):
        if cancel():
            raise Cancelled("Cancelled between inference requests")
        try:
            response = client.post(ENDPOINT, data=payload)
        except httpx.TransportError as error:
            if attempt == attempts - 1:
                raise InferenceError(f"IEDB request failed: {error}") from error
            delay = 2**attempt
        else:
            if response.status_code not in {429, 500, 502, 503, 504}:
                if response.is_error:
                    raise InferenceError(f"IEDB HTTP {response.status_code}: {response.text[:200]}")
                return response.text
            if attempt == attempts - 1:
                raise InferenceError(f"IEDB temporarily unavailable (HTTP {response.status_code})")
            retry = response.headers.get("Retry-After", "")
            if retry.isdigit():
                delay = float(retry)
            elif retry:
                try:
                    deadline = parsedate_to_datetime(retry)
                    delay = max(0, (deadline - datetime.datetime.now(datetime.UTC)).total_seconds())
                except (ValueError, TypeError):
                    delay = 2**attempt
            else:
                delay = 2**attempt
            # If a server asks for a long pause, leave a resumable run instead of retrying too early.
            if delay > 60:
                raise InferenceError(f"IEDB rate limited; retry after {delay:g} seconds")
        for _ in range(math.ceil(delay * 4)):
            if cancel():
                raise Cancelled("Cancelled during retry wait")
            sleep(0.25)
    raise InferenceError("Retry limit exhausted")


def predict(
    cache,
    sequences,
    alleles,
    lengths,
    batch_size=128,
    cached_only=False,
    cancel=lambda: False,
    progress=lambda message: None,
):
    lookup = cache.lookup(alleles)
    batches = plan_batches(sequences, alleles, lengths, lookup, batch_size)
    statistics = {
        "planned_requests": len(batches),
        "completed_requests": 0,
        "cached_predictions": len(lookup),
        "remaining_requests": len(batches),
    }
    progress(f"Binding: {len(batches)} requests needed; {len(lookup)} cached peptide–HLA scores available")
    if cached_only:
        return lookup, statistics, "partial" if batches else "complete", []
    errors, status = [], "complete"
    with httpx.Client(timeout=120, follow_redirects=True) as client:
        for number, (allele, length, batch) in enumerate(batches, 1):
            fasta = "".join(f">sequence_{i}\n{s}\n" for i, s in enumerate(batch, 1))
            payload = {"method": METHOD, "sequence_text": fasta, "allele": allele, "length": length}
            request_id = hashlib.sha256(
                json.dumps({"endpoint": ENDPOINT, **payload}, sort_keys=True).encode()
            ).hexdigest()
            path = cache.directory / "requests" / request_id
            try:
                if cancel():
                    raise Cancelled("Cancelled between batches")
                if (path / "metadata.json").exists():
                    metadata = json.loads((path / "metadata.json").read_text())
                    if metadata.get("payload") != payload or metadata.get("endpoint") != ENDPOINT:
                        raise InferenceError("Cached request provenance mismatch")
                    response = (path / "response.tsv").read_text()
                    if hashlib.sha256(response.encode()).hexdigest() != metadata["response_sha256"]:
                        raise InferenceError("Cached response hash mismatch")
                else:
                    response = request_batch(client, payload, cancel)
                rows = validate_response(response, batch, allele, length)
            except Cancelled as error:
                status = "cancelled"
                errors.append(str(error))
                break
            except (InferenceError, ValueError, KeyError) as error:
                status = "partial"
                errors.append(str(error))
                failure_path = cache.directory / "failures" / request_id
                failure_path.mkdir(parents=True, exist_ok=True)
                write_json(
                    failure_path / "request.json",
                    {"endpoint": ENDPOINT, "payload": payload, "error": str(error)},
                )
                break
            path.mkdir(parents=True, exist_ok=True)
            (path / "response.tsv").write_text(response)
            write_json(
                path / "metadata.json",
                {
                    "endpoint": ENDPOINT,
                    "payload": payload,
                    "method": METHOD,
                    "software_version": __version__,
                    "retrieved_at_utc": datetime.datetime.now(datetime.UTC).isoformat(),
                    "response_sha256": hashlib.sha256(response.encode()).hexdigest(),
                    "sequence_mapping": {str(i): s for i, s in enumerate(batch, 1)},
                },
            )
            cache.save(rows, request_id)
            for row in rows:
                lookup[(row["peptide"], allele)] = {**row, "request_id": request_id}
            statistics.update(completed_requests=number, remaining_requests=len(batches) - number)
            progress(f"Binding {number}/{len(batches)} complete · {allele} · {length}-mers")
            if not cancel():
                time.sleep(0.25)
    return lookup, statistics, status, errors
