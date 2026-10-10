import httpx
import pytest

from neovax import iedb
from neovax.iedb import BindingCache, Cancelled, InferenceError, predict, request_batch, validate_response

ALLELE = "HLA-A*02:01"


def response_for(sequences, length=9):
    rows = ["allele\tseq_num\tstart\tend\tlength\tpeptide\tic50\tpercentile_rank"]
    for index, sequence in enumerate(sequences, 1):
        for start in range(len(sequence) - length + 1):
            rows.append(
                f"{ALLELE}\t{index}\t{start + 1}\t{start + length}\t{length}\t{sequence[start : start + length]}\t12.3\t0.2"
            )
    return "\n".join(rows) + "\n"


def test_response_coordinates_completeness_and_scores():
    sequences = ["ACDEFGHIKL"]
    text = response_for(sequences)
    assert len(validate_response(text, sequences, ALLELE, 9)) == 2
    for invalid in [
        "\n".join(text.splitlines()[:-1]),
        text.replace("ACDEFGHIK", "AAAAAAAAA"),
        text.replace("12.3", "nan"),
        text + text.splitlines()[1] + "\n",
    ]:
        with pytest.raises(InferenceError):
            validate_response(invalid, sequences, ALLELE, 9)


def test_temporary_retry_and_rate_limit():
    codes = iter([503, 429, 200])
    client = httpx.Client(
        transport=httpx.MockTransport(
            lambda request: httpx.Response(next(codes), text="success", headers={"Retry-After": "0"})
        )
    )
    assert request_batch(client, {}, lambda: False, sleep=lambda _: None) == "success"
    client.close()
    client = httpx.Client(
        transport=httpx.MockTransport(lambda request: httpx.Response(429, headers={"Retry-After": "120"}))
    )
    with pytest.raises(InferenceError, match="retry after 120"):
        request_batch(client, {}, lambda: False, sleep=lambda _: None)
    client.close()
    with pytest.raises(Cancelled):
        request_batch(None, {}, lambda: True)


def test_cancel_resume_and_cache_reuse(tmp_path, monkeypatch):
    cache = BindingCache(tmp_path)
    calls = []
    cancelled = [False]

    def fake(client, payload, cancel):
        sequences = iedb.parse_fasta(payload["sequence_text"])
        calls.append(sequences)
        return response_for(sequences)

    monkeypatch.setattr(iedb, "request_batch", fake)
    monkeypatch.setattr(iedb.time, "sleep", lambda _: None)

    def progress(message):
        if "1/2 complete" in message:
            cancelled[0] = True

    sequences = ["ACDEFGHIK", "LMNPQRSTV"]
    _, stats, status, _ = predict(
        cache, sequences, [ALLELE], [9], batch_size=1, cancel=lambda: cancelled[0], progress=progress
    )
    assert status == "cancelled" and stats["completed_requests"] == 1
    cancelled[0] = False
    lookup, stats, status, _ = predict(cache, sequences, [ALLELE], [9], batch_size=1)
    assert status == "complete" and len(lookup) == 2 and len(calls) == 2
    _, stats, status, _ = predict(cache, sequences, [ALLELE], [9], cached_only=True)
    assert status == "complete" and stats["planned_requests"] == 0 and len(calls) == 2
    cache.close()


def test_truncated_response_is_not_cached(tmp_path, monkeypatch):
    cache = BindingCache(tmp_path)
    monkeypatch.setattr(iedb, "request_batch", lambda *args: "allele\tseq_num\n")
    lookup, _, status, errors = predict(cache, ["ACDEFGHIK"], [ALLELE], [9])
    assert status == "partial" and errors and not lookup and not cache.lookup([ALLELE])
    assert list((tmp_path / "failures").glob("*/request.json"))
    cache.close()


def test_cached_only_missing_is_partial(tmp_path):
    cache = BindingCache(tmp_path)
    _, stats, status, _ = predict(cache, ["ACDEFGHIK"], [ALLELE], [9], cached_only=True)
    assert status == "partial" and stats["remaining_requests"] == 1
    cache.close()


def test_duplicate_peptide_scores_must_agree():
    sequences = ["ACDEFGHIK", "ACDEFGHIK"]
    text = response_for(sequences)
    lines = text.splitlines()
    lines[-1] = lines[-1].replace("12.3", "99.0")
    with pytest.raises(InferenceError, match="inconsistent scores"):
        validate_response("\n".join(lines), sequences, ALLELE, 9)
